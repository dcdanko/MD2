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
    'id',
    'citation',
    'evidence',
    ]
	
REGEX_TAXANOMY = [
    'eubacterium',
    'bacterium ',
    'sp.',
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
    ]

def reduce_col(file):
    """Remove empty columns, ids and taxonomy columns"""
    drop_col = file.dropna(axis='columns', how='all')
    drop_col.columns = map(str.lower, drop_col.columns)
    col_names = ['class']
    for reg in REGEX_COLUMN:
        col_names.extend(list(drop_col.filter(regex=reg)))
    col_names.remove('taxonomic_id')
    drop_col = drop_col.drop(columns=col_names, axis=1)		 
    return drop_col

def reduce_row(file):
    """Remove duplicate rows and subspecies"""
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
    """"""
    pass


