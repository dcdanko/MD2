
import click
import pandas as pd
import numpy as np

from microbe_directory.constants import DOMAINS
from microbe_directory.dataset_stats import (
    verify_column_names,
    column_compare,
)


@click.group('stats')
def stats():
    pass


@stats.command('table')
@click.option('-s', '--sep', default='\t')
@click.option('-o', '--outfile', default='-', type=click.File('w'))
@click.argument('table', type=click.File('r'))
def cli_md2_table_stats(sep, outfile, table):
    """Collate some statistics for a microbe directory table."""
    tbl = pd.read_csv(table, index_col=0)
    pf = lambda key, val: print(f'{key}{sep}{val}', file=outfile)
    pf(f'Number of Taxa', tbl.shape[0])
    pf(f'Number of columns', tbl.shape[1])
    pf(f'Columns', ','.join(tbl.columns))
    for col_name in tbl.columns:
        col = tbl[col_name]
        pf(f'{col_name} Filled', tbl.shape[0] - sum(col.isna()))
        pf(f'{col_name} Unique Vals', len(col.unique()))
