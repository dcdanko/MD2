# Programatically infill missing fields in a MD2 Table
import click
import pandas as pd

from microbe_directory.taxa_tree import NCBITaxaTree
from microbe_directory.constants import (
    RANK_LIST,
    SPORE_FORMING_GENERA,
    NON_SPORE_FORMING_GENERA,
    PHYLUM_GRAM_STAINS,
)


def infill_bacterial_fields(table, verbose=True):
    """Return a copy of table with certain fields filled in.

    Currently fills in
     - Spore forming status for genus and below
     - Gram stain based on stain of ancestor
     - Presence of relevant genus in the EMP
    """
    def pif(el):
        if verbose:
            click.echo(el, err=True)

    taxa_tree = NCBITaxaTree.parse_files()
    pif('built taxa tree.')

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

    table['spore_forming'] = table['spore_forming'].fillna(table.index.to_series().map(spore_forming))
    pif('filled in spore status.')
    table['gram_stain'] = table['gram_stain'].fillna(table.index.to_series().map(gram_stain))
    pif('filled in gram stain')
    emp_cols = [col for col in table.columns if col.startswith('emp') or col.startswith('count')]
    pif(f'filling in {len(emp_cols)} EMP columns')
    for col in emp_cols:
        filler_func = infill_emp_col_func(table[col])
        table[col] = table[col].fillna(table.index.to_series().map(filler_func))
        pif(f'filled in {col}')
        
    return table
