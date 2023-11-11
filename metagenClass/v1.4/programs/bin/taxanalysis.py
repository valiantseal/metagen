import pandas as pd
from ete3 import NCBITaxa

# Initialize NCBITaxa
ncbi = NCBITaxa()

# Load data
krak_blast_conf_reads_df = pd.read_csv('output/krakBlastConfReads.csv', sep=',', quotechar='"')
lca_numbers_df = pd.read_csv('output/lca_numbers.csv')
read_label_library_df = pd.read_csv('read_label_library.csv')

# Convert 'Main_Read' in lca_numbers_df to string and process it
lca_numbers_df['Main_Read'] = lca_numbers_df['Main_Read'].astype(str)
lca_numbers_df['Read'] = lca_numbers_df['Main_Read'].apply(lambda x: x.split('/')[0])

# Duplicate the entries to represent both ends of the paired-end reads
lca_numbers_df_1 = lca_numbers_df.copy()
lca_numbers_df_2 = lca_numbers_df.copy()
lca_numbers_df_1['Read'] = lca_numbers_df_1['Read'] + '/1'
lca_numbers_df_2['Read'] = lca_numbers_df_2['Read'] + '/2'
lca_numbers_expanded = pd.concat([lca_numbers_df_1, lca_numbers_df_2])

# Merge DataFrames on 'Read'
merged_df = read_label_library_df.merge(krak_blast_conf_reads_df, left_on='ReadID', right_on='Read', how='left')
merged_df = merged_df.merge(lca_numbers_expanded, on='Read', how='left')

# Select relevant columns and create a new DataFrame
final_df = merged_df[['Read', 'Label', 'BID1', 'KID1', 'LCA']].copy()
final_df.rename(columns={'BID1': 'BLAST', 'KID1': 'Kraken'}, inplace=True)

def is_match_at_level(taxid, label_taxid, level):
    """Check if a taxid matches or is nested within the label_taxid at a given level."""
    if pd.isna(taxid) or pd.isna(label_taxid):
        return False
    if taxid == label_taxid:
        return True
    lineage = ncbi.get_lineage(taxid)
    label_lineage = ncbi.get_lineage(label_taxid)
    label_rank = ncbi.get_rank([label_taxid]).get(label_taxid, None)
    if label_rank != level:
        label_lineage = [ancestor for ancestor in label_lineage if ncbi.get_rank([ancestor])[ancestor] == level]
        label_taxid = label_lineage[0] if label_lineage else None
    return label_taxid in lineage

def calculate_precision(df, label_col, pred_col, level):
    """Calculate precision at a given taxonomic level."""
    matches = df.apply(lambda x: is_match_at_level(x[pred_col], x[label_col], level), axis=1)
    TP = sum(matches)
    FN = sum(~matches & ~pd.isna(df[label_col]))
    return TP / (TP + FN) if TP + FN > 0 else 0

# Levels to check
levels = ['species', 'genus', 'family']

def calculate_precision(df, label_col, pred_col, level):
    """Calculate precision at a given taxonomic level, accounting for NAs."""
    matches = df.apply(lambda x: is_match_at_level(x[pred_col], x[label_col], level) if pd.notna(x[pred_col]) else False, axis=1)
    TP = sum(matches)
    FN = sum(~matches & ~pd.isna(df[label_col]))
    return TP / (TP + FN) if TP + FN > 0 else 0

# Levels to check
levels = ['species', 'genus', 'family']

# Calculate precision for each level
precision_data = []
for level in levels:
    kraken_precision = calculate_precision(final_df, 'Label', 'Kraken', level)
    blast_precision = calculate_precision(final_df, 'Label', 'BLAST', level)
    pipeline_precision = calculate_precision(final_df, 'Label', 'LCA', level)
    precision_data.append([level.capitalize(), kraken_precision, blast_precision, pipeline_precision])

# Create a DataFrame for the precision data
metrics_df = pd.DataFrame(precision_data, columns=['Precision Level', 'Kraken', 'BLAST', 'Pipeline'])

# Write the DataFrame to a CSV file
metrics_df.to_csv('output/pipeline_metrics.csv', index=False)
