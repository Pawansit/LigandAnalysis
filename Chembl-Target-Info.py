from chembl_webresource_client.new_client import new_client
import pandas as pd
import math
from rdkit.Chem import PandasTools


def BioactivityIC50(CHEMBLID):
	bioact = bioactivities.filter(target_chembl_id = CHEMBLID) \
                      .filter(type = 'IC50') \
                      .only('activity_id', 'assay_description', 'assay_type', \
                             'molecule_chembl_id','type', 'units', 'relation', 'value', \
                            'target_chembl_id')
	bioact_df = pd.DataFrame.from_records(bioact)
	return(bioact_df)

def Compoundlist(record):
	#record = bioact_df
	cmpd_id_list = list(record['molecule_chembl_id'])
	compound_list = compounds.filter(molecule_chembl_id__in = cmpd_id_list) \
                         .only('molecule_chembl_id','molecule_structures')
	compound_df = pd.DataFrame.from_records(compound_list)
	compound_df = compound_df.drop_duplicates('molecule_chembl_id', keep = 'first')
	for i, cmpd in compound_df.iterrows():
		if compound_df.loc[i]['molecule_structures'] != None:
			compound_df.loc[i]['molecule_structures'] = cmpd['molecule_structures']['canonical_smiles']

	output_df = pd.merge(record[['target_chembl_id','molecule_chembl_id','units','value']], compound_df, on='molecule_chembl_id')
	output_df = output_df.rename(columns= {'molecule_structures':'smiles', 'value':'IC50'})
	return(output_df)

targets = new_client.target
compounds = new_client.molecule
bioactivities = new_client.activity
#chembl_id = 'CHEMBL2052035'				## Define the input CHEMBL targte ID

listID = ['CHEMBL4222','CHEMBL2052035','CHEMBL4968','CHEMBL3559677','CHEMBL6147','CHEMBL3486','CHEMBL2169724']
#listID = ['CHEMBL4222']
for i in listID:
	print(i)
	bioact = BioactivityIC50(i)
	print(bioact.shape[0])
	if bioact.shape[0]  > 0:
		complist = Compoundlist(bioact)		
		name = i+".csv"
		complist.to_csv("/home/pawan/Downloads/"+name+"")
		print(complist)


