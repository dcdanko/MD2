# Add composite fields to a MD2 Table
import click
import pandas as pd

from microbe_directory.taxa_tree import NCBITaxaTree


COMPOSITES = {
    'composite_polar': [
        'emp_polluted_polar_coastal_sediments_counts',
        'emp_mcmurdo_lakes_counts',
        'emp_antarctic_counts',
    ],
    'composite_desert': [
        'emp_atacama_desert_counts',
        'emp_arabic_counts',
        'emp_thar_counts',
    ],
    'composite_deep_water': [
        'emp_great_lake_microbiome_counts',
        'emp_tara_oceans_counts',
        'count_tara',
    ],
    'composite_soils': [
        'andisols',
        'gelisols',
        'vertisols',
        'mollisols',
        'inceptisols',
        'alfisols',
        'ultisols',
        'sand_rock_ice',
        'entisols',
    ],
}


def rectify_presence(values):

    def count_in(*els):
        count = 0
        for val in values:
            val = str(val).lower()
            for el in els:
                count += 1 if el in val else 0
        return count

    never = count_in('never')
    rarely = count_in('rarely')

    if never == len(values):
        return 'Never Observed'
    if (never + rarely) == len(values):
        return 'Rarely Observed'
    if count_in('often', 'always') >= (len(values) - 1):
        return 'Often Observed'
    return 'Observed in some'


def add_composite_fields(tbl, verbose=True):
    for col_name, composite_cols in COMPOSITES.items():
        tbl[col_name] = tbl[composite_cols].apply(
            lambda row: rectify_presence(list(row)),
            axis=1,
        )
    return tbl
