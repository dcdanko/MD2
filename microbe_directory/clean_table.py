import pandas as pd
import numpy as np

REGEX_COLUMN = [
    'specie',
    'genus', 
    'family', 
    'order', 
    'phylum', 
    'domain', 
    'kingdom', 
    'organism', 
    'unnamed',
    'citation',
    'evidence',
    'microbe_id',
    'otu',
    ' id',
    'non_null',
    'pmid',
    'sample type', 
    ]
	
REGEX_COUNT_COL = [
    'count',
    'sols',
    '_id_',
    ]
	

def file_clean(file):
    column_clean = reduce_col(file)
    file_renamed = modify_dataset_value(column_clean)	
    new_file = file_renamed.groupby(['scientific_name', 'taxonomic_id', 'rank'], as_index=True).agg(lambda x: ( ', '.join( repr(e) for e in list(set(x)))))
    return new_file.replace('nan, ', '', regex=True).replace('\'', '').replace('[', '').replace(']', '').applymap(lambda x: x.replace('\'', ''))
	
	
def reduce_col(file):
    """Remove empty columns, ids and taxonomy columns"""
    drop_col = file.dropna(axis='columns', how='all')
    drop_col.columns = map(str.lower, drop_col.columns)
    col_names = list()
    for reg in REGEX_COLUMN:    
        col_names.extend(list(drop_col.filter(regex=reg)))
    drop_col = drop_col.drop(columns=col_names, axis=0)
    drop_col.columns = rename_col(drop_col)	
    final_file = rename_MD1_tables(drop_col)
    return final_file


def rename_col(file):
    """Convert column names as per snakelowercase standards"""
    file.columns = file.columns.str.strip().str.lower().str.replace(' ', '_').str.replace('_x', '_').str.replace('_y', '_').str.replace('.', '_').str.replace(',', '_')
    return file.columns

def rename_MD1_tables(file):
    """Replace numeric entries from MD1"""
    file['gram_stain'] = file['gram_stain'].replace([0, 1, 2], ['Negative', 'Positive', 'Intermediate'])
    file['extreme_environment'] = file['extreme_environment'].replace([0, 1], ['Mesophiles', 'Extremophile'])
    file['antimicrobial_susceptibility'] = file['antimicrobial_susceptibility'].replace([0, 1], ['Maybe Not', 'Sometimes'])
    file['biofilm_forming'] = file['biofilm_forming'].replace([0, 1], ['Never', 'Always'])
    file['animal_pathogen'] = file['animal_pathogen'].replace([0, 1], ['Maybe Not', 'Sometimes'])
    file['plant_pathogen'] = file['plant_pathogen'].replace([0, 1], ['Maybe Not', 'Sometimes'])
    file['microbiome_location'] = file['microbiome_location'].replace([0, 1], ['Maybe', 'Sometimes'])
    file['spore_forming'] = file['spore_forming'].replace([0, 1], ['Never', 'Always'])
    file = file.rename(columns={'microbiome_location': 'human_pathogen'})
    return file
	
def modify_dataset_value(file):
    """Convert datasets with count values to interpretable values"""
    if 'drylands' in file.columns:
        file['drylands'] = file['drylands'].replace([0, 1], ['Not Observed', 'Observed'])
        file['low_productivity'] = file['low_productivity'].replace([0, 1], ['Not Observed', 'Observed'])
        file['low_ph'] = file['low_ph'].replace([0, 1], ['Not Observed', 'Observed'])
        file['high_ph'] = file['high_ph'].replace([0, 1], ['Not Observed', 'Observed'])
    for reg in REGEX_COUNT_COL:
        regex_columns = [cols for cols in file.columns if reg in cols]
        final_file = clean_count_datasets(file, regex_columns)
    return final_file
	
def clean_count_datasets(file, regex_list):
    """Logic for count conversion"""
    for reg in regex_list:
        file[reg] = pd.to_numeric(file[reg], errors='coerce')
        file[reg] = (file[reg] / file[reg].sum()) * 100
        file[reg] = file[reg].mask((file[reg]>0) & (file[reg]<=25), 2)
        file[reg] = file[reg].mask((file[reg]>25) & (file[reg]<75), 3)
        file[reg] = file[reg].mask((file[reg]>=75) & (file[reg]<100), 4)
        file[reg] = file[reg].replace([0, 2, 3, 4, 100], ['Never Observed', 'Rarely Observed', 'Fairly Observed', 'Mostly Observed', 'Always Observed'])
    return file
