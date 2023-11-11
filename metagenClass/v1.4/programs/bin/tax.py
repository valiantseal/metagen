import os
from ete3 import NCBITaxa
import pandas as pd
import itertools
from concurrent.futures import ProcessPoolExecutor
from collections import defaultdict

# Initialize NCBITaxa
ncbi = NCBITaxa()

# Function to get the rank of a taxid
def get_rank(taxid):
    return ncbi.get_rank([taxid]).get(taxid, None)

# Function to check if a taxid is classified
def is_classified(taxid):
    lineage = ncbi.get_lineage(taxid)
    names = ncbi.get_taxid_translator(lineage)
    return "unclassified" not in " ".join(names.values()).lower()

# Function to get the LCA for a list of taxids
def get_lca(taxids):
    try:
        combinations = list(itertools.combinations(taxids, 2))
        lcas = [(pair, ncbi.get_topology(pair)) for pair in combinations if pair[0] and pair[1]]
        if not lcas:
            return None, None, None
        lcas.sort(key=lambda x: ncbi.get_rank([x[1].taxid])[x[1].taxid])
        return lcas[0][0][0], lcas[0][0][1], lcas[0][1].taxid
    except KeyError as e:
        print(f"KeyError encountered: {e}")
        return None, None, None

# Function to process each row in parallel
def process_row(row):
    bid_used, kid_used, lca_taxid = determine_final_lca(row)
    return row['Main_Read'], bid_used, kid_used, lca_taxid

# Pre-process dataframe to merge paired-end reads into single rows
def preprocess_dataframe(df):
    df['Main_Read'] = df['Read'].str.extract(r'(.*?)(?:/\d+)?$')
    df_1 = df[df['Read'].str.endswith('/1')].add_suffix('_1')
    df_2 = df[df['Read'].str.endswith('/2')].add_suffix('_2')
    merged_df = pd.merge(df_1, df_2, left_on='Main_Read_1', right_on='Main_Read_2', how='outer')
    return merged_df.rename(columns={'Main_Read_1': 'Main_Read'})

# Function to determine the final LCA with a preference for classified species
def determine_final_lca(row):
    all_bids = [row[f"BID{i}_1"] for i in range(1, 4) if pd.notnull(row[f"BID{i}_1"])] + \
               [row[f"BID{i}_2"] for i in range(1, 4) if pd.notnull(row[f"BID{i}_2"])]
    all_kids = [row[f"KID{i}_1"] for i in range(1, 4) if pd.notnull(row[f"KID{i}_1"])] + \
               [row[f"KID{i}_2"] for i in range(1, 4) if pd.notnull(row[f"KID{i}_2"])]
    all_taxids = all_bids + all_kids

    lcas = []
    for taxid1, taxid2 in itertools.combinations(all_taxids, 2):
        if taxid1 and taxid2:
            try:
                lca = ncbi.get_topology([taxid1, taxid2])
                lcas.append((taxid1, taxid2, lca.taxid))
            except KeyError:
                pass

    species_lcas = [(lca, get_rank(lca[2])) for lca in lcas if get_rank(lca[2]) == 'species']
    species_lcas.sort(key=lambda x: (x[1], not is_classified(x[0][2])))

    if species_lcas:
        return species_lcas[0][0]
    elif lcas:
        lcas.sort(key=lambda x: ncbi.get_rank([x[2]])[x[2]])
        return lcas[0]
    else:
        return None, None, None

# Function to fetch the lineage of a given taxid
def fetch_lineage(taxid):
    try:
        lineage = ncbi.get_lineage(taxid)
        return ncbi.get_rank(lineage), ncbi.get_taxid_translator(lineage)
    except ValueError:
        return {}, {}

def build_hierarchy(lca_counts):
    """Build a hierarchy based on parent-child relationships."""
    parent_of = {}
    taxon_count = defaultdict(int)
    for lca, count in lca_counts.items():
        try:
            lineage = ncbi.get_lineage(lca)
            prev_taxid = None
            for taxid in lineage:
                if prev_taxid is not None:
                    parent_of[taxid] = prev_taxid
                prev_taxid = taxid
                taxon_count[taxid] += count
        except ValueError:
            print(f"Error processing LCA: {lca}")
    return parent_of, taxon_count

def print_hierarchy_to_file(parent_of, taxon_count, taxid, file, indent=0):
    """Recursively print the hierarchy to a file, excluding taxa with a count of 0."""
    children = [child for child, parent in parent_of.items() if parent == taxid]
    
    # Calculate the taxon's own count by subtracting the children's counts
    own_count = taxon_count[taxid] - sum(taxon_count[child] for child in children)

    if own_count > 0 and taxid != 1:  # Skip taxa with 0 count and the root
        taxon_name = ncbi.get_taxid_translator([taxid]).get(taxid, "Unknown")
        file.write(" " * indent + f"{taxon_name} (TaxID: {taxid}) {own_count}\n")

    for child in children:
        print_hierarchy_to_file(parent_of, taxon_count, child, file, indent + 4)

def translate_taxids_to_names(taxids):
    # Filter out NaN values and convert remaining IDs to integers
    valid_taxids = [int(taxid) for taxid in taxids if pd.notnull(taxid)]
    names_dict = ncbi.get_taxid_translator(valid_taxids)

    # Create a list of names, using "Unknown" for NaN values
    return [names_dict.get(int(taxid), "Unknown") if pd.notnull(taxid) else "Unknown" for taxid in taxids]

def main(dataframe):
    preprocessed_df = preprocess_dataframe(dataframe)
    with ProcessPoolExecutor(max_workers=8) as executor:
        results = list(executor.map(process_row, [row for _, row in preprocessed_df.iterrows()]))

    # Create DataFrame with numbers (taxonomic IDs) and save it
    results_df = pd.DataFrame(results, columns=['Main_Read', 'BID_Used', 'KID_Used', 'LCA'])
    results_df.to_csv('output/lca_numbers.csv', sep=',', index=False)

    # Generate summary statistics for unique LCAs and build hierarchy
    lca_counts = results_df['LCA'].value_counts().to_dict()
    parent_of, taxon_count = build_hierarchy(lca_counts)

    # Write hierarchy to file in output directory
    with open('output/taxonomy_hierarchy.txt', 'w') as file:
        roots = set(parent_of.values()) - set(parent_of.keys())
        for root in roots:
            root_name = ncbi.get_taxid_translator([root]).get(root, "Unknown")
            file.write(f"Root: {root_name} (TaxID: {root})\n")
            print_hierarchy_to_file(parent_of, taxon_count, root, file)

    # Translate taxonomic IDs to names for results_df
    results_df['BID_Used'] = translate_taxids_to_names(results_df['BID_Used'])
    results_df['KID_Used'] = translate_taxids_to_names(results_df['KID_Used'])
    results_df['LCA'] = translate_taxids_to_names(results_df['LCA'])

    # Save the DataFrame with taxonomic names
    results_df.to_csv('output/lca_results.csv', sep=',', index=False)

    return results_df

if __name__ == "__main__":
    dataframe = pd.read_csv('output/krakBlastConfReads.csv', sep=',', quotechar='"')
    final_combined_df = main(dataframe)
