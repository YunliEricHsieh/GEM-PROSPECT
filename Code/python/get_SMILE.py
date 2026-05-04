import pandas as pd
import requests
import time

# Read your Excel file
file_path = 'Results/SPOT/potential_transporter.xlsx'
df = pd.read_excel(file_path)

# Identify unique metabolite names to minimize API calls
unique_metabolites = df['MetaboliteID'].dropna().unique().tolist()
smiles_map = {}

print(f"Fetching SMILES for {len(unique_metabolites)} metabolites...\n")

for name in unique_metabolites:
    try:
        # STEP 1: Get the CID
        cid_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/cids/JSON"
        cid_res = requests.get(cid_url, timeout=10)
        
        if cid_res.status_code == 200:
            cid = cid_res.json()['IdentifierList']['CID'][0]
            
            # STEP 2: Request ALL likely structural properties
            # Note: We add 'SMILES' and 'InChI' as broad backups
            prop_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/IsomericSMILES,CanonicalSMILES,SMILES,InChI/JSON"
            prop_res = requests.get(prop_url, timeout=10)
            
            if prop_res.status_code == 200:
                props = prop_res.json()['PropertyTable']['Properties'][0]
                
                # FALLBACK LOGIC: Try keys in order of preference
                # PubChem is moving some IsomericSMILES to just 'SMILES' in their RDF/JSON responses
                smiles = (props.get('IsomericSMILES') or 
                          props.get('SMILES') or 
                          props.get('CanonicalSMILES'))
                
                if smiles:
                    smiles_map[name] = smiles
                    print(f"✅ {name} (CID:{cid}): {smiles[:30]}...")
                elif props.get('InChI'):
                    # Final fallback: If SMILES is missing, record the InChI
                    smiles_map[name] = props.get('InChI')
                    print(f"⚠️ {name} (CID:{cid}): SMILES missing, used InChI instead")
                else:
                    smiles_map[name] = "Structure data missing"
                    print(f"❌ {name} (CID:{cid}): No structural fields found")
            else:
                smiles_map[name] = "Property API Error"
        else:
            smiles_map[name] = "Name not found"
            print(f"❓ {name}: Not found")
            
        time.sleep(0.2)
        
    except Exception as e:
        print(f"🔴 {name}: {str(e)}")
        smiles_map[name] = "Error"

# Save results
df['Metabolite'] = df['MetaboliteID'].map(smiles_map)
df.to_excel(file_path, index=False)
print(f"\nSaved results to {file_path}")