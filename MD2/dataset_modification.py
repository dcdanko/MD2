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
	
def convert_taxa_tree(filename, biom_file, metadata_file, feature_name, subtext):
    otu_file = pd.read_csv(filename, index_col=0)
    otu_file.replace('__', ' ', inplace=True, regex=True)
    otu_file.replace(r'(.*)\s(.*)', r'\2', inplace=True, regex=True)
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
    biom_new = metasub_process(biom.transpose(), metadata_file, feature_name, subtext)
    return biom_new
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
	
def convert_taxa_tree(filename, biom_file, metadata_file, feature_name, subtext):
    otu_file = pd.read_csv(filename, index_col=0)
    otu_file.replace('__', ' ', inplace=True, regex=True)
    otu_file.replace(r'(.*)\s(.*)', r'\2', inplace=True, regex=True)
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
    biom_new = metasub_process(biom.transpose(), metadata_file, feature_name, subtext)
    return biom_new
