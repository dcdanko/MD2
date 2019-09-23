import pandas as pd
import numpy as np

def metasub_process(filename, metadata_filename, feature_name='city', subtext='metasub'):
    metasub_merged = {}
    metadata = pd.read_csv(metadata_filename, index_col=0)
    metadata = metadata.loc[filename.index]
    feature = metadata[feature_name]
    factorized, name_map = pd.factorize(feature)
    for i in range(0, len(name_map)):
        metasub_df = filename[factorized == i]
        species_count = metasub_df.count()
        species_array = (species_count.values/metasub_df.shape[0]) * 100
        if i == 0:
            data = np.reshape(species_array, (len(species_array)),1)
            metasub_merged = pd.DataFrame(data, index = list(filename.columns.values), columns=[str(subtext + '_' + name_map[i])]) 
        else:
            metasub_merged[str(subtext + '_' + name_map[i])] = np.reshape(species_array, (len(species_array)),1)
    metasub_merged[metasub_merged.apply(lambda x: (x>0) & (x<=25))] = 2
    metasub_merged[metasub_merged.apply(lambda x: (x>25) & (x<75))] = 3
    metasub_merged[metasub_merged.apply(lambda x: (x>=75) & (x<100))] = 4
    metasub_merged = metasub_merged.replace([0, 2, 3, 4, 100], ['Never Observed', 'Rarely Observed', 'Fairly Observed', 'Mostly Observed', 'Always Observed'])
    return metasub_merged
	
def convert_taxa_tree(filename, biom_file, metadata_file, feature_name):
    otu_file = pd.read_csv(filename, index_col=0)
    otu_file.replace('__', ' __', inplace=True, regex=True)
    otu_file.replace(r'(.*)\s(__.*)', r'\2', inplace=True, regex=True)
    otu_file.replace('__', '', inplace=True, regex=True)
    otu_to_microbe = {}
    for index, rows in otu_file.iterrows():
        for value in ['Species', 'Genus', 'Family', 'Order', 'Class', 'Phylum', 'Kingdom']:
            if rows[value] != '':
                otu_to_microbe[index] = {'OTU': index, 'Microbe_Name': rows[value]} 
                break
    annotated_otu = pd.DataFrame.from_dict(otu_to_microbe, orient='index')
    biom = pd.read_csv(biom_file, index_col=0)
    biom  = biom.merge(annotated_otu, left_on='#OTU ID' , right_on='OTU')
    biom = biom[biom.columns.drop(list(biom.filter(regex='OTU')))]
    biom = biom.set_index('Microbe_Name')
    biom = biom.groupby('Microbe_Name').sum()
    biom_new = biom.transpose()
    biom_new.index.name = "StudyID"
    biom_new = biom_new.reset_index()
    metadata = pd.read_csv(metadata_file)
    metadata_sliced = metadata[['#SampleID', feature_name]]
    metadata_sliced = metadata_sliced.merge(biom_new, right_on= 'StudyID', left_on= '#SampleID')
    metadata_sliced  = metadata_sliced[metadata_sliced.columns.drop(['StudyID', '#SampleID'])]
    metadata_sliced = metadata_sliced.groupby(feature_name).sum()
    metadata_sliced.index.name = feature_name
    metadata_sliced = metadata_sliced.fillna(0).astype('float')
    metadata_2 = (metadata_sliced.T / metadata_sliced.T.sum()).T
    metadata_2[metadata_2.apply(lambda x: (x>0) & (x<=.25))] = 2
    metadata_2[metadata_2.apply(lambda x: (x>.25) & (x<.75))] = 3
    metadata_2[metadata_2.apply(lambda x: (x>=.75) & (x<1))] = 4
    metadata_2 = metadata_2.replace([0, 2, 3, 4, 1], ['Never Observed', 'Rarely Observed', 'Fairly Observed', 'Mostly Observed', 'Always Observed'])
    return metadata_2.transpose()

	
def taxa_to_organism(otu_file):
    """Convert to the most relevant level of organism for ITS/16S/18S/amplicons studies"""
    taxa_level = ['organism', 'species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom', 'domain']
    regex_column = ['id', 'otu', 'reference']
    otu_file.replace(r'(.*)\s(__.*)', r'\2', inplace=True, regex=True)
    otu_file.replace('__', '', inplace=True, regex=True)
    for levels in taxa_level: 
        if levels not in otu_file.columns:
            otu_file[levels] = ''
    for index, rows in otu_file.iterrows():
        value = ''
        for value in taxa_level:
            if pd.isnull(str(rows[value])) == False and rows[value]!='':
                otu_file.loc[index, 'organism_name'] = rows[value]
                break
    for reg in regex_column:    
        taxa_level.extend(list(otu_file.filter(regex=reg)))
    otu_file = otu_file.drop(columns=taxa_level, axis=1)
    otu_file = otu_file.groupby('organism_name').sum()
    otu_file = otu_file.reset_index()
    otu_file.replace('\'', '', inplace=True, regex=True)
    return otu_file
