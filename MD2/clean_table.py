import pandas as pd
import numpy as np


REGEX_COLUMN = [
    'specie',
    'genus', 
    'family', 
    'order', 
    'class_', 
    'phylum', 
    'domain', 
    'kingdom', 
    'organism', 
    'unnamed',
    ' id',
    'citation',
    'evidence',
    ]
	
REGEX_TAXANOMY = [
    'eubacterium',
    'bacterium ',
    'sp.',
    'phytoplasma',
    'uncultured',
    'unidentified',
    'symbiont',
    'methanotroph',
    'clinical sample',
    'isolate',
    'clone',
    'associated',
    'vouchered',
    'fungal endophyte',
    'cf.',
    's.l.',
    'sect.',
    'aff.',
    'strain',
    '16Sr',
	'str.',
    'culture',
    'enrichment',
    'human',
    'diazotroph',
    ' of ',
    'et al.',
    'planctomycete ',
    ' bacterium',
    'phytoplasma'
    ]

def reduce_col(isvirus, file):
    """Remove empty columns, ids and taxonomy columns"""
    drop_col = file.dropna(axis='columns', how='all')
    drop_col.columns = map(str.lower, drop_col.columns)
    col_names = ['class', 'pmid', 'id']
    for reg in REGEX_COLUMN:
        col_names.extend(list(drop_col.filter(regex=reg)))
    if isvirus == 'True': 
        col_names.remove('id')
    drop_col = drop_col.drop(columns=col_names, axis=1)
    drop_col.columns = rename_col(drop_col)	
    final_file = rename_MD1_tables(drop_col)
    return final_file

def reduce_row(isvirus, file):
    """Remove duplicate rows and subspecies"""
    if isvirus == 'True': 
        REGEX_TAXANOMY.remove('bacterium ')
        REGEX_TAXANOMY.remove('human')
    for reg in REGEX_TAXANOMY:
        filter = file['scientific_name'].str.contains(reg)
        file = file[~filter]
    new_file = file.groupby(['scientific_name', 'taxonomic_id', 'rank'], as_index=True).agg(lambda x: ( ', '.join( repr(e) for e in list(set(x)))))
    return new_file.replace('nan, ', '', regex=True).replace('\'', '').replace('[', '').replace(']', '').applymap(lambda x: x.replace('\'', ''))
	
def rename_col(file):
    """Convert column names as per snakelowercase standards"""
    file.columns = file.columns.str.strip().str.lower().str.replace(' ', '_').str.replace('_x', '_').str.replace('_y', '_').str.replace('.', '_').str.replace(',', '_')
    return file.columns

def rename_MD1_tables(file):
    """Replace numeric entries from MD1"""
    file['gram_stain'] = file['gram_stain'].replace([0, 1, 2], ['Negative', 'Positive', 'Intermediate'])
    file['extreme_environment'] = file['extreme_environment'].replace([0, 1], ['Mesophiles', 'Extremophile'])
    file['antimicrobial_susceptibility'] = file['antimicrobial_susceptibility'].replace([0, 1], ['No', 'Yes'])
    file['biofilm_forming'] = file['biofilm_forming'].replace([0, 1], ['No', 'Yes'])
    file['animal_pathogen'] = file['animal_pathogen'].replace([0, 1], ['No', 'Yes'])
    file['plant_pathogen'] = file['plant_pathogen'].replace([0, 1], ['No', 'Yes'])
    file['microbiome_location'] = file['microbiome_location'].replace([0, 1], ['No', 'Yes'])
    file['spore_forming'] = file['spore_forming'].replace([0, 1], ['No', 'Yes'])
    file = file.rename(columns={'microbiome_location': 'human_disease_causing'})
    return file


