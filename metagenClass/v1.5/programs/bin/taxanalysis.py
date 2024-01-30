import pandas as pd

# Load the full results and read label library data
full_results_df = pd.read_csv('output/fullresults.csv')
read_label_library_df = pd.read_csv('read_label_library.csv')

# Extract the parent read ID from the ReadID column in the read label library
read_label_library_df['Main_Read'] = read_label_library_df['ReadID'].str.extract(r'(.*?)(?:/\d+)?$')

# Merge the dataframes on the Main_Read column
merged_df = pd.merge(full_results_df, read_label_library_df[['Main_Read', 'Label']], on='Main_Read', how='left')

# Save the merged dataframe to a new CSV file
merged_df.to_csv('output/fullresults_with_labels.csv', index=False)
