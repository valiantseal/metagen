import os
from ete3 import NCBITaxa
import pandas as pd
import itertools
from concurrent.futures import ProcessPoolExecutor
from collections import defaultdict

ncbi = NCBITaxa()

os.makedirs('output', exist_ok=True)

def get_rank(taxid):
    return ncbi.get_rank([taxid]).get(taxid, None)

def is_classified(taxid):
    lineage = ncbi.get_lineage(taxid)
    names = ncbi.get_taxid_translator(lineage)
    return "unclassified" not in " ".join(names.values()).lower()

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

def process_row(row):
    bid_used, kid_used, lca_taxid = determine_final_lca(row)
    return (row['Main_Read'], row['BID1_1'], row['BID2_1'], row['BID3_1'], row['BID1_2'], row['BID2_2'], 
            row['BID3_2'], row['KID1_1'], row['KID2_1'], row['KID3_1'], row['KID1_2'], row['KID2_2'], 
            row['KID3_2'], bid_used, kid_used, lca_taxid, row['Sample'])

def preprocess_dataframe(df):
    df['Main_Read'] = df['Read'].str.extract(r'(.*?)(?:/\d+)?$')
    df_1 = df[df['Read'].str.endswith('/1')].add_suffix('_1')
    df_2 = df[df['Read'].str.endswith('/2')].add_suffix('_2')
    merged_df = pd.merge(df_1, df_2, left_on='Main_Read_1', right_on='Main_Read_2', how='outer')
    merged_df = merged_df.rename(columns={'Main_Read_1': 'Main_Read'})
    merged_df['Sample'] = merged_df['Sample_1'].combine_first(merged_df['Sample_2'])
    return merged_df.drop(columns=['Sample_1', 'Sample_2'])

standard_ranks = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom', 'superkingdom']

def elevate_taxid(taxid):
    rank = get_rank(taxid)
    if rank in ['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom', 'superkingdom']:
        return taxid
    else:
        lineage = ncbi.get_lineage(taxid)
        ranks = ncbi.get_rank(lineage)
        for tid in reversed(lineage):
            current_rank = ranks[tid]
            if current_rank in ['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom', 'superkingdom']:
                return tid
    return taxid  # Return original if no elevation is possible

def determine_final_lca(row):
    all_taxids = set()
    for i in range(1, 4):
        all_taxids.update(filter(pd.notnull, [row.get(f'BID{i}_1'), row.get(f'BID{i}_2'), row.get(f'KID{i}_1'), row.get(f'KID{i}_2')]))

    # Function to get rank or lineage-based rank position
    def get_rank_position(taxid):
        rank = get_rank(taxid)
        if rank in standard_ranks:
            return standard_ranks.index(rank)
        else:
            lineage = ncbi.get_lineage(taxid)
            ranks = ncbi.get_rank(lineage)
            for tid in reversed(lineage):
                current_rank = ranks.get(tid, None)
                if current_rank in standard_ranks:
                    return standard_ranks.index(current_rank)
            # Penalize non-standard ranks without a standard rank in their lineage
            return len(standard_ranks) + 1

    if len(all_taxids) == 1:
        taxid = all_taxids.pop()
        elevated_taxid = elevate_taxid(taxid)
        return elevated_taxid, elevated_taxid, elevated_taxid

    lcas = []
    for taxid1, taxid2 in itertools.combinations(all_taxids, 2):
        lca = ncbi.get_topology([taxid1, taxid2])
        lcas.append((taxid1, taxid2, lca.taxid))

    valid_lcas = []
    for lca in lcas:
        taxid = elevate_taxid(lca[2])
        lineage = ncbi.get_lineage(taxid)
        ranks = ncbi.get_rank(lineage)
        if any(rank in standard_ranks for rank in ranks.values()):
            valid_lcas.append((lca[0], lca[1], taxid))

    if valid_lcas:
        valid_lcas.sort(key=lambda x: get_rank_position(x[2]))
        return valid_lcas[0]
    else:
        return None, None, None  # Return this if no valid LCA is found

def fetch_lineage(taxid):
    try:
        lineage = ncbi.get_lineage(taxid)
        return ncbi.get_rank(lineage), ncbi.get_taxid_translator(lineage)
    except ValueError:
        return {}, {}

def translate_taxids_to_names(taxids):
    valid_taxids = [int(taxid) for taxid in taxids if pd.notnull(taxid)]
    names_dict = ncbi.get_taxid_translator(valid_taxids)
    return [names_dict.get(int(taxid), "Unknown") if pd.notnull(taxid) else "Unknown" for taxid in taxids]

def main():
    process_dir = './process'
    dfs = []
    for sample_dir in os.listdir(process_dir):
        sample_path = os.path.join(process_dir, sample_dir)
        blast_file = os.path.join(sample_path, 'blastResTop_3.tsv')
        kraken_file = os.path.join(sample_path, 'krakenSelVirReads.tsv')
        blast_df = pd.read_csv(blast_file, sep='\t', quotechar='"')
        kraken_df = pd.read_csv(kraken_file, sep=',', quotechar='"')
        merged_df = pd.merge(blast_df, kraken_df, on='Read', how='outer')
        merged_df['Sample'] = sample_dir
        dfs.append(merged_df)
    dataframe = pd.concat(dfs, ignore_index=True)

    preprocessed_df = preprocess_dataframe(dataframe)
    with ProcessPoolExecutor(max_workers=os.cpu_count()) as executor:
        results = list(executor.map(process_row, [row for _, row in preprocessed_df.iterrows()]))
    results_df = pd.DataFrame(results, columns=['Main_Read', 'BID1_1', 'BID2_1', 'BID3_1', 'BID1_2', 'BID2_2', 'BID3_2', 'KID1_1', 'KID2_1', 'KID3_1', 'KID1_2', 'KID2_2', 'KID3_2', 'BID_Used', 'KID_Used', 'LCA', 'Sample'])
    results_df['LCA_Name'] = translate_taxids_to_names(results_df['LCA'])
    lca_counts_df = results_df.groupby(['Sample', 'LCA']).size().reset_index(name='Count')
    lca_names_df = results_df[['LCA', 'LCA_Name']].drop_duplicates()
    lca_summary_df = pd.merge(lca_counts_df, lca_names_df, on='LCA', how='left')
    lca_summary_df = lca_summary_df[['Sample', 'LCA_Name', 'LCA', 'Count']]
    lca_summary_df.to_csv('output/lca_summary.csv', sep=',', index=False)
    results_df.to_csv('output/fullresults.csv', sep=',', index=False)
    return results_df, lca_summary_df

if __name__ == '__main__':
    main()
