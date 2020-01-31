import pandas as pd
import numpy as np


def metasub_process(taxa_tbl, metadata, feature_name='city', prefix='metasub_', min_size=16):
    shared_keys = set(taxa_tbl.index) & set(metadata.index)
    metadata, taxa_tbl = metadata.loc[shared_keys], taxa_tbl.loc[shared_keys]
    taxa_tbl = taxa_tbl > 0  # convert to presense/absence matrix
    taxa_tbl[feature_name] = metadata[feature_name]
    taxa_tbl = taxa_tbl.groupby(feature_name).filter(lambda df: df.shape[0] >= min_size)
    taxa_tbl = taxa_tbl.groupby(feature_name).mean()  # gives the prevalence of each taxa per group
    taxa_tbl = taxa_tbl.T  # features in columns, taxa in rows
    taxa_tbl.columns = taxa_tbl.columns.map(lambda el: prefix + el)

    def relabel_ranges(fraction):
        if fraction < 0.000001:
            return 'Never Observed'
        if fraction <= 0.25:
            return 'Rarely Observed'
        if fraction <= 0.75:
            return 'Observed'
        if fraction < 0.98:
            return 'Often Observed'
        if fraction <= 1.0:
            return 'Always Observed'
        assert False, f'Invalid fraction: {fraction}'

    taxa_tbl = taxa_tbl.applymap(relabel_ranges)
    return taxa_tbl


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
