import requests
import pandas as pd
import time

# --- CONFIG ---
csv_file = 'Results/CataPro/enzymes_and_substrates.csv'

# 1. Read Data
df = pd.read_csv(csv_file)

# 2. Fetch Protein Sequences (UniProt)
unique_enzymes = set(str(e).split(' or ')[0].strip() for e in df['Enzyme'].dropna())
protein_map = {}
print(f"🧬 Fetching sequences for {len(unique_enzymes)} unique enzymes...")

for uniID in unique_enzymes:
    try:
        res = requests.get(f'https://rest.uniprot.org/uniprotkb/{uniID}', timeout=10)
        if res.status_code == 200:
            protein_map[uniID] = res.json()['sequence']['value']
            print(f"✅ Protein: {uniID}")
        time.sleep(0.1)
    except:
        print(f"🔴 Failed: {uniID}")

# 3. Fetch SMILES (PubChem Two-Step)
unique_metabolites = df['Substrate'].dropna().unique().tolist()
smiles_map = {}
print(f"\n🧪 Fetching SMILES for {len(unique_metabolites)} metabolites...")

for name in unique_metabolites:
    try:
        cid_res = requests.get(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/cids/JSON", timeout=10)
        if cid_res.status_code == 200:
            cid = cid_res.json()['IdentifierList']['CID'][0]
            prop_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/IsomericSMILES,CanonicalSMILES,SMILES/JSON"
            prop_res = requests.get(prop_url, timeout=10)
            if prop_res.status_code == 200:
                props = prop_res.json()['PropertyTable']['Properties'][0]
                smiles = props.get('IsomericSMILES') or props.get('SMILES') or props.get('CanonicalSMILES')
                smiles_map[name] = smiles
                print(f"✅ SMILES: {name}")
        time.sleep(0.2)
    except:
        print(f"❓ SMILES Not Found: {name}")

# 4. Map back and Save
def get_clean_id(val):
    return str(val).split(' or ')[0].strip() if pd.notna(val) else None

df['Sequence'] = df['Enzyme'].apply(get_clean_id).map(protein_map)
df['smiles'] = df['Substrate'].map(smiles_map)

df.to_csv(csv_file, index=False)
print(f"\n✅ File ready for manual curation: {csv_file}")