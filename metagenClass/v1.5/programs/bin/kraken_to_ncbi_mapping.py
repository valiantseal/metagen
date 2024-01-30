import csv
import os

taxdb_file_path = '~/extraVol/krakenUniq/viral/taxDB'
taxdb_file_path = os.path.expanduser(taxdb_file_path)

# Function to create a dictionary mapping KrakenDB IDs to their NCBI parent IDs
def create_kraken_to_ncbi_mapping(file_path):
    kraken_to_ncbi = {}
    parent_mapping = {}

    # Read the file and build the mappings
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split()
            if len(parts) >= 2:
                kraken_id, parent_id = int(parts[0]), int(parts[1])
                # Store the parent mapping for all IDs
                parent_mapping[kraken_id] = parent_id
                # If the Kraken ID is over 1,000,000,000, start the mapping process
                if kraken_id > 1000000000:
                    # Loop to find the first parent ID under 1,000,000,000
                    while parent_id > 1000000000:
                        parent_id = parent_mapping.get(parent_id, parent_id)
                    kraken_to_ncbi[kraken_id] = parent_id

    return kraken_to_ncbi

# Create the mapping dictionary
kraken_to_ncbi_mapping = create_kraken_to_ncbi_mapping(taxdb_file_path)

# Save the dictionary to a CSV file
with open('kraken_to_ncbi_mapping.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['KrakenDB', 'NCBI_Parent'])  # Writing header
    for key, value in kraken_to_ncbi_mapping.items():
        writer.writerow([key, value])
