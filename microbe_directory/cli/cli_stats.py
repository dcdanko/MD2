
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


@stats.command('stats-file')
@click.option('--microbe', default='bacteria', type=click.Choice(DOMAINS),
              help='Biological domain of the table')
@click.option('-o', '--outfile', default='-', type=click.File('w'))
@click.argument('table', type=click.File('r'))
def clean_file(microbe, outfile, table):
    """Outputs the statistics of column that needs to be filled."""
    file_name = pd.read_csv(table)
    final_file, stats_file = verify_column_names(file_name, microbe, default=None)
    final_file.to_csv(outfile)
    stats_file.to_csv(str(microbe + '_stats.csv'))


@stats.command('column-compare')
@click.option('-o', '--outfile', default='-', type=click.File('w'))
@click.argument('table_1', type=click.File('r'))
@click.argument('table_2', type=click.File('r'))
def dataset_column_compare(outfile, table_1, table_2):
    """Compare the before and after statistics of a file."""
    table_1 = pd.read_csv(table_1)
    table_2 = pd.read_csv(table_2)
    table_1 = table_1.replace(r'^\s*$', np.nan, regex=True)
    table_2 = table_2.replace(r'^\s*$', np.nan, regex=True)
    stats1, stats2 = column_compare(table_1, table_2)
    stats = pd.concat([stats1, stats2], sort=False, names=['title', 'before', 'after'], axis=1)
    stats.to_csv(outfile)
