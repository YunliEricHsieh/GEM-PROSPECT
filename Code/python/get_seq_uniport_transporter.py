import pandas as pd
import requests
import time
import re

# Read the Excel file
file_path = 'Results/SPOT/potential_transporter.xlsx'
df = pd.read_excel(file_path)

# Extract all unique UniProt IDs, even from "ID1 or ID2" strings
def get_primary_ids(enzyme_col):
    primary_ids = set()
    for entry in enzyme_col.dropna():
        # Split by ' or ', '/', ';', or ',' and take the first part
        parts = re.split(r'\s+or\s+|\s*[/,;]\s*', str(entry), flags=re.IGNORECASE)
        first_id = parts[0].strip()
        if first_id:
            primary_ids.add(first_id)
    return list(primary_ids)

unique_ids = get_primary_ids(df['Enzymes'])

# 3. Fetch sequences from UniProt API
sequence_map = {}
print(f"Fetching sequences for {len(unique_ids)} unique primary enzymes...")

for uniID in unique_ids:
    try:
        response = requests.get(f'https://rest.uniprot.org/uniprotkb/{uniID}', timeout=15)
        if response.status_code == 200:
            data = response.json()
            seq = data.get('sequence', {}).get('value', None)
            if seq:
                sequence_map[uniID] = seq
                print(f"✅ {uniID}")
            else:
                sequence_map[uniID] = "Sequence field missing in API"
        else:
            sequence_map[uniID] = f"Error: Status {response.status_code}"
        
        # Rate limiting to stay within UniProt guidelines
        time.sleep(0.2) 
    except Exception as e:
        sequence_map[uniID] = f"Connection failed: {e}"

# 4. Map the sequence to the table using only the first ID in the cell
def fetch_first_sequence(enzyme_str):
    if pd.isna(enzyme_str):
        return None
    
    # Split the string and take the first ID
    parts = re.split(r'\s+or\s+|\s*[/,;]\s*', str(enzyme_str), flags=re.IGNORECASE)
    first_id = parts[0].strip()
    
    return sequence_map.get(first_id, "ID not found/processed")

# Apply and Save
df['Protein'] = df['Enzymes'].apply(fetch_first_sequence)
df.to_excel(file_path, index=False)

print(f"\nDone! Updated table saved to {file_path}")