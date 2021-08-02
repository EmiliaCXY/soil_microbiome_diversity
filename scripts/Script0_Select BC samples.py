import pandas as pd

# importing data
metadata = pd.read_csv('MICB_421_Soil_Metadata.tsv',sep='\t')
manifest = pd.read_csv('soil_manifest.txt',sep='\t')

# filtering based on our specs
bc_metadata = metadata[metadata['Region']=='British Columbia']
bc_metadata = bc_metadata[bc_metadata['Horizon'] == 'O horizon']
bc_metadata['pH'] = bc_metadata['pH'].astype(float)
bc_metadata = bc_metadata[bc_metadata['pH']>0]
bc_sample_id = list(bc_metadata['#SampleID'])
manifest = manifest[manifest['sample-id'].isin(bc_sample_id)]

# saving the manifest table
manifest.to_csv('bc_o_ph_soil_manifest.txt',index=False, sep='\t')