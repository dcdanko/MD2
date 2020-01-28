
import pandas as pd

from os.path import join, dirname


def final_table_path(name):
    return join(
        dirname(__file__),
        'stored_final_tables',
        name,
    )


def parse(name):
    return pd.read_csv(
        final_table_path(name),
        index_col=0
    )


def bacteria():
    return parse('bacteria.csv.gz')


def virus():
    return parse('viruses.csv.gz')


def eukaryote():
    return parse('eukaryota.csv.gz')
