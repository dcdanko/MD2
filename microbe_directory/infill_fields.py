# Programatically infill missing fields in a MD2 Table
import click
import pandas as pd

from microbe_directory.taxa_tree import NCBITaxaTree
from microbe_directory.constants import (
    RANK_LIST,
    SPORE_FORMING_GENERA,
    NON_SPORE_FORMING_GENERA,
)


def infill_bacterial_fields(table):
    """Return a copy of table with certain fields filled in.

    Currently fills in
     - Spore forming status for genus and below
     - Gram stain based on stain of ancestor
     - Presence of relevant genus in the EMP
    """
    table = table.copy(deep=True)
    emp_cols = [
        col for col in table.columns
        if col.startswith('emp') or col.startswith('count')
    ]
    taxa_tree = NCBITaxaTree.parse_files()
    
    for index, row in table.iterrows():
        scientific_name = index  # row['scientific_name']
        try:
            rank = taxa_tree.rank(scientific_name)
        except KeyError:
            continue
        genera = taxa_tree.genus(scientific_name, default=None)
        levels_below_genus = RANK_LIST.index('genus') - RANK_LIST.index(rank)

        # Update Spore_Forming in a Top-Down approach based on Genus
        if levels_below_genus >= 0:
            if genera in SPORE_FORMING_GENERA:
                table.loc[index, 'spore_forming'] = 'Always'
            elif genera in NON_SPORE_FORMING_GENERA:
                table.loc[index, 'spore_forming'] = 'Never'

        # Update Gram_Stain in a Mixed approach based on Genus
        ancestor_list = taxa_tree.ancestors_list(scientific_name)
        if ancestor_list:
            for ancestor in ancestor_list[:-1:-1]:
                ancestor_stain = table.loc[ancestor]['gram_stain']
                if 'Positive' == ancestor_stain:
                    table.loc[scientific_name, 'gram_stain'] = 'Positive'
                    break
                elif 'Negative' == ancestor_stain:
                    table.loc[scientific_name, 'gram_stain'] = 'Negative'
                    break

        # Update EMP dataset based on Top-Down approach
        if levels_below_genus >= 1:
            for col in emp_cols:
                if pd.isnull(row[col]):
                    try:
                        genus_value = table.loc[genera, col]
                        if genus_value and str(genus_value).lower() not in ['nan', 'na', 'n/a', '']:
                            table.loc[scientific_name, col] = f'{genus_value} in Genus'
                    except KeyError:
                        pass
    return table
