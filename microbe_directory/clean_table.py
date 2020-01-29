import pandas as pd
import click

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


def file_clean(tbl):
    column_clean = reduce_col(tbl)
    tbl_renamed = modify_dataset_value(column_clean)
    new_tbl = tbl_renamed.drop_duplicates(
        subset=['scientific_name', 'taxonomic_id', 'rank']
    )  # not likely to actually be duplicates
    new_tbl = new_tbl.replace('nan, ', '', regex=True)
    new_tbl = new_tbl.replace('\'', '').replace('[', '').replace(']', '')
    new_tbl = new_tbl.applymap(lambda x: x.replace('\'', '') if isinstance(x, str) else x)
    return new_tbl


def clean_columns(tbl):
    unnamed = [el for el in tbl.columns if 'unnamed' in el.lower()]
    tbl = tbl.drop(columns=unnamed)
    halo = [el for el in tbl.columns if 'halotolerance_classification' in el.lower()]
    if halo:
        tbl['halotolerance'] = tbl[halo[0]].iloc[:, 0].map(
            lambda el: 'Moderate' if 'Moderate' in str(el) else str(el).strip()
        )
        tbl = tbl.drop(columns=halo)
    return tbl


def reduce_col(tbl):
    """Remove empty columns, ids and taxonomy columns"""
    drop_col = tbl.dropna(axis='columns', how='all')
    drop_col.columns = map(str.lower, drop_col.columns)
    col_names = list()
    for reg in REGEX_COLUMN:
        col_names.extend(list(drop_col.filter(regex=reg)))
    drop_col = drop_col.drop(columns=col_names, axis=0)
    drop_col.columns = rename_col(drop_col)
    final_tbl = rename_MD1_tables(drop_col)
    return final_tbl


def rename_col(tbl):
    """Convert column names as per snakelowercase standards"""
    tbl.columns = tbl.columns.str.strip().str.lower()
    pairs = [(' ', '_'), ('_x', '_'), ('_y', '_'), ('.', '_'), (',', '_')]
    for a, b in pairs:
        tbl.columns = tbl.columns.str.replace(a, b)
    return tbl.columns


def rename_MD1_tables(tbl):
    """Replace numeric entries from MD1"""
    md1_cols = [
        ('gram_stain', ['Negative', 'Positive', 'Intermediate']),
        ('extreme_environment', ['Mesophiles', 'Extremophile']),
        ('antimicrobial_susceptibility', ['Maybe Not', 'Sometimes']),
        ('biofilm_forming', ['Never', 'Always']),
        ('animal_pathogen', ['Maybe Not', 'Sometimes']),
        ('plant_pathogen', ['Maybe Not', 'Sometimes']),
        ('microbiome_location', ['Maybe', 'Sometimes']),
        ('spore_forming', ['Never', 'Always'])
    ]
    for col_name, new_vals in md1_cols:
        old_vals = [i for i in range(len(new_vals))]
        tbl[col_name] = tbl[col_name].replace(old_vals, new_vals)
    tbl = tbl.rename(columns={'microbiome_location': 'human_commensal'})
    return tbl


def modify_dataset_value(tbl):
    """Convert datasets with count values to interpretable values"""
    if 'drylands' in tbl.columns:
        for col in ['drylands', 'low_productivity', 'low_ph', 'high_ph']:
            tbl[col] = tbl[col].replace([0, 1], ['Not Observed', 'Observed'])
    for reg in REGEX_COUNT_COL:
        regex_columns = [cols for cols in tbl.columns if reg in cols]
        final_tbl = clean_count_datasets(tbl, regex_columns)
    return final_tbl


def clean_count_datasets(tbl, regex_list):
    """Logic for count conversion"""
    for reg in regex_list:
        tbl[reg] = pd.to_numeric(tbl[reg], errors='coerce')
        tbl[reg] = (tbl[reg] / tbl[reg].sum()) * 100
        tbl[reg] = tbl[reg].mask((tbl[reg] > 0) & (tbl[reg] <= 25), 2)
        tbl[reg] = tbl[reg].mask((tbl[reg] > 25) & (tbl[reg] < 75), 3)
        tbl[reg] = tbl[reg].mask((tbl[reg] >= 75) & (tbl[reg] < 100), 4)
        tbl[reg] = tbl[reg].replace(
            [0, 2, 3, 4, 100],
            ['Never Observed', 'Rarely Observed', 'Sometimes Observed',
             'Usually Observed', 'Always Observed']
        )
    return tbl
