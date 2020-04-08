from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import Descriptors
import pandas as pd

suppl = Chem.SDMolSupplier('/home/pawan/Downloads/pawan-50000.sdf')

def AtomNumber(m):
    M = []
    for atom in m.GetAtoms():
        M.append(atom)
    
    return(len(M))

def readMOl(A):
	data = pd.DataFrame([])
	
	for mol in suppl:
		if mol is not None:
			name = mol.GetProp('chembl_id')
			core= MurckoScaffold.GetScaffoldForMol(mol)
			MW = Descriptors.MolWt(core)
			RoundMW = round(MW,2)
			S = Chem.MolToSmiles(core)
			AN = AtomNumber(core)
			data = data.append(pd.DataFrame({'CompName':name,'Smiles':S,'MolWt':RoundMW,'AtomNumber':AN},index=[0]),ignore_index=True)
	
	
	return(data)
	


df = readMOl(suppl)
df= df.drop_duplicates("Smiles",keep = 'first')
df=df[df.AtomNumber > 10]
print(df.shape)
print(df.head(20))
