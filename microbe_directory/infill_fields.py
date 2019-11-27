# Programatically infill missing fields in a MD2 Table

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
    taxa_tree = NCBITaxaTree.parse_files()
    for index, row in table.iterrows():
        scientific_name = row['scientific_name']
        rank = taxa_tree.rank(scientific_name)
        genera = taxa_tree.genus(scientific_name, default=None)
        levels_below_genus = RANK_LIST.index('genus') - RANK_LIST.index(rank)

        # Update Spore_Forming in a Top-Down approach based on Genus
        if levels_below_genus >= 0:
            if genera in SPORE_FORMING_GENERA:
                table.loc[index, 'spore_forming'] = 'Always'
            elif genera in NON_SPORE_FORMING_GENERA:
                table.loc[index, 'spore_forming'] = 'Never'

        # Update Gram_Stain in a Mixed approach based on Genus
        ancestor_list = taxa_tree.ancestors_list('genus', scientific_name, default=None)
        if ancestor_list:
            ancestor_stains = table[table['scientific_name'].isin(ancestor_list)]['gram_stain']
            if 'Positive' in ancestor_stains:
                table.loc[table['scientific_name'].isin(ancestor_list), 'gram_stain'] = 'Positive'
            elif 'Negative' in ancestor_stains:
                table.loc[table['scientific_name'].isin(ancestor_list), 'gram_stain'] = 'Negative'

        # Update EMP dataset based on Top-Down approach
        if levels_below_genus >= 1:
            emp_cols = [
                col for col in table.columns
                if col.startswith('emp') or col.startswith('count')
            ]
            for col in emp_cols:
                if pd.isnull(row[col]):
                    genus_value = table[table['scientific_name'] == genera][col].values
                    if len(genus_value) and not pd.isnull(genus_value[0]):
                        table.loc[table['scientific_name'] == scientific_name, col] = f'{genus_value[0]} in Genus'
    return table
