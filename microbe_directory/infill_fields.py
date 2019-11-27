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


def infill_bacterial_fields(table):
    """Return a copy of table with certain fields filled in.

    Currently fills in
     - Spore forming status for genus and below
     - Gram stain based on stain of ancestor
     - Presence of relevant genus in the EMP
    """
    taxa_tree = NCBITaxaTree.parse_files()

    def spore_forming(scientific_name):
        genus = taxa_tree.genus(scientific_name, default=None)
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
        phylum = taxa_tree.phylum(scientific_name, default=None)
        try:
            return PHYLUM_GRAM_STAINS[phylum]
        except KeyError:
            pass
        return None

    def infill_emp_col_func(col):
        def infill_emp_col(scientific_name):
            genus = taxa_tree.genus(scientific_name, default=None)
            if genus:
                try:
                    genus_value = col[genus]
                    if genus_value and str(genus_value).lower() not in ['nan', 'na', 'n/a', '']:
                        return f'{genus_value} in Genus'
                except KeyError:
                    pass
            return None
        return infill_emp_col

    table['spore_forming'] = table['spore_forming'].fillna(table.index.map(spore_forming))
    table['gram_stain'] = table['gram_stain'].fillna(table.index.map(gram_stain))
    for col in [col for col in table.columns if col.startswith('emp') or col.startswith('count')]:
        filler_func = infill_emp_col_func(col)
        table[col] = table[col].fillna(table.index.map(filler_func))

    return table
