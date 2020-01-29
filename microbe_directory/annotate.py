

import pandas as pd

from .final_tables import (
    bacteria,
    virus,
    eukaryote,
)


def annotate_taxa(taxa, kind='bacteria'):
    """Return a pandas dataframe with annotations for taxa in the list.

    Drops any columns that is null.
    """
    tbl = bacteria
    if kind[0] == 'e':
        tbl = eukaryote
    elif kind[0] == 'v':
        tbl = virus
    md_annotations = pd.DataFrame.from_dict({
            'taxa': taxa,
        },
        orient='columns'
    ).set_index('taxa').join(tbl(), how='left')
    md_annotations = md_annotations.dropna(axis=1, how='all')
    return md_annotations
