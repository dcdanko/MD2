# Programatically infill missing fields in a MD2 Table
import click
import pandas as pd

from microbe_directory.taxa_tree import NCBITaxaTree
from microbe_directory.constants import (
    RANK_LIST,
    SPORE_FORMING_GENERA,
    NON_SPORE_FORMING_GENERA,
    PHYLUM_GRAM_STAINS,
    PSYCHRO,
    PSYCHRO_GENUS,
    RADIO,
    RADIO_GENUS,
)




def infill_bacterial_fields(tbl, verbose=True):
    """Return a copy of tbl with certain fields filled in.

    Currently fills in
     - Spore forming status for genus and below
     - Gram stain based on stain of ancestor
     - Presence of relevant genus in the EMP
    """

    def spore_forming(scientific_name):
        try:
            genus = taxa_tree.genus(scientific_name, default=None)
        except KeyError:
            return None
        if genus:
            if genus in SPORE_FORMING_GENERA:
                return 'Always'
            elif genus in NON_SPORE_FORMING_GENERA:
                return 'Never'
        return None

    def gram_stain(scientific_name):
        try:
            rank = taxa_tree.rank(scientific_name)
        except KeyError:
            return None
        try:
            phylum = taxa_tree.phylum(scientific_name, default=None)
            return PHYLUM_GRAM_STAINS[phylum]
        except KeyError:
            pass
        return None

    def infill_emp_col_func(col):
        def infill_emp_col(scientific_name):
            try:
                genus = taxa_tree.genus(scientific_name, default=None)
                if genus:
                    genus_value = col[genus]
                    if genus_value and str(genus_value).lower() not in ['nan', 'na', 'n/a', '']:
                        return f'{genus_value} in Genus'
            except KeyError:
                pass
            return None
        return infill_emp_col

    def psychro(scientific_name):
        try:
            temp = float(tbl.loc[scientific_name, 'optimal_temperature'])
            if temp < 10:
                return 'Psychrophilic'
        except TypeError:
            pass
        except ValueError:
            pass
        if scientific_name in PSYCHRO:
            return 'Psychrophilic'
        elif taxa_tree.genus(scientific_name, default='') in PSYCHRO_GENUS:
            return 'Psychrophilic Genus'
        return None

    def radiophilic(scientific_name):
        if scientific_name in RADIO:
            return 'Radiophilic'
        elif taxa_tree.genus(scientific_name, default='') in RADIO_GENUS:
            return 'Radiophilic Genus'
        return None

    def pif(el):
        if verbose:
            click.echo(el, err=True)

    taxa_tree = NCBITaxaTree.parse_files()
    pif('built taxa tree.')

    tbl['spore_forming'] = tbl['spore_forming'].fillna(tbl.index.to_series().map(spore_forming))
    pif('filled in spore status.')

    tbl['gram_stain'] = tbl['gram_stain'].fillna(tbl.index.to_series().map(gram_stain))
    pif('filled in gram stain')

    tbl['psychrophilic'] = tbl.index.to_series().map(psychro)
    pif('filled in psychrophilia')

    tbl['radiophilic'] = tbl.index.to_series().map(radiophilic)
    pif('filled in radiophilia')

    emp_cols = [col for col in tbl.columns if col.startswith('emp') or col.startswith('count')]
    pif(f'filling in {len(emp_cols)} EMP columns')
    for col in emp_cols:
        filler_func = infill_emp_col_func(tbl[col])
        tbl[col] = tbl[col].fillna(tbl.index.to_series().map(filler_func))
        pif(f'filled in {col}')

    return tbl
