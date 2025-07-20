import pandas as pd
import os
import zipfile

# Read the metadata
metadata_path = '/Users/vinayak/Documents/gemini-cli/BixBench_data/BixBench_metadata.csv'
df = pd.read_csv(metadata_path)

# Target IDs
target_ids = ['bix-5', 'bix-10']

for short_id in target_ids:
    # Find the row for the target ID
    row = df[df['short_id'] == short_id]
    if row.empty:
        print(f"Could not find metadata for {short_id}")
        continue

    # Get the data folder name (which is also the zip file name)
    data_folder = row['data_folder'].iloc[0]
    zip_filename = f"{data_folder}.zip" # Corrected: removed the extra .zip
    zip_path = os.path.join('/Users/vinayak/Documents/gemini-cli/BixBench_data', zip_filename)

    # Define the destination directory for unzipping
    dest_dir = f'/Users/vinayak/Documents/gemini-cli/bixbench_analysis/{short_id}'
    
    print(f"Processing {short_id}:")
    print(f"  - Zip file: {zip_path}")
    print(f"  - Destination: {dest_dir}")

    # Unzip the file
    if os.path.exists(zip_path):
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_ref.extractall(dest_dir)
        print(f"  - Successfully unzipped.")
        
        # Remove any .ipynb files
        for root, dirs, files in os.walk(dest_dir):
            for file in files:
                if file.endswith('.ipynb'):
                    ipynb_path = os.path.join(root, file)
                    os.remove(ipynb_path)
                    print(f"  - Removed notebook: {ipynb_path}")
    else:
        print(f"  - ERROR: Zip file not found at {zip_path}")