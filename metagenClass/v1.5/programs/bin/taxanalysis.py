import pandas as pd
from ete3 import NCBITaxa

ncbi = NCBITaxa()

full_results_df = pd.read_csv('output/fullresults.csv')
read_label_library_df = pd.read_csv('read_label_library.csv')
read_label_library_df['Main_Read'] = read_label_library_df['ReadID'].str.extract(r'(.*?)(?:/\d+)?$')
merged_df = pd.merge(full_results_df, read_label_library_df[['Main_Read', 'Label']], on='Main_Read', how='left')

def get_lineage_info(taxid):
    if pd.isna(taxid):
        return {}
    lineage = ncbi.get_lineage(taxid)
    ranks = ncbi.get_rank(lineage)
    names = ncbi.get_taxid_translator(lineage)
    return {ranks[taxid]: names[taxid] for taxid in lineage}

def is_correct_species_or_below(classifier_lineage, label_lineage):
    return 'species' in label_lineage and label_lineage['species'] in classifier_lineage.values()

def is_correct_family_or_below(classifier_lineage, label_lineage):
    return 'family' in label_lineage and label_lineage['family'] in classifier_lineage.values()

def calculate_stats(row):
    stats = {
        'Kraken_Species_PID': 0, 'Kraken_Species_NID': 0, 'Kraken_Family_PID': 0,
        'Pipeline_Species_PID': 0, 'Pipeline_Species_NID': 0, 'Pipeline_Family_PID': 0
    }

    label_lineage = get_lineage_info(row['Label'])

    for classifier_key in ['KID1_1', 'KID1_2', 'LCA']:
        classifier_taxid = row[classifier_key]
        classifier_lineage = get_lineage_info(classifier_taxid)

        if is_correct_species_or_below(classifier_lineage, label_lineage):
            if 'KID' in classifier_key:
                stats['Kraken_Species_PID'] += 1
            else:
                stats['Pipeline_Species_PID'] += 1
        elif 'species' in classifier_lineage:
            if 'KID' in classifier_key:
                stats['Kraken_Species_NID'] += 1
            else:
                stats['Pipeline_Species_NID'] += 1

        if is_correct_family_or_below(classifier_lineage, label_lineage):
            if 'KID' in classifier_key:
                stats['Kraken_Family_PID'] += 1
            else:
                stats['Pipeline_Family_PID'] += 1

    return pd.Series(stats)

stats_df = merged_df.apply(calculate_stats, axis=1)
summary_stats = stats_df.groupby(merged_df['Label']).sum()
summary_stats.reset_index(inplace=True)
summary_stats.rename(columns={'Label': 'Virus'}, inplace=True)
summary_stats.to_csv('summary_stats.csv', index=False)
