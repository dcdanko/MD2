
import pandas as pd

from os.path import join, dirname


def final_table_path(name):
    return join(
        dirname(__file__),
        'stored_final_tables',
        name,
    )


def parse(name, ind=0):
    return pd.read_csv(
        final_table_path(name),
        index_col=ind
    )


def bacteria():
    """Return a pandas dataframe of bacterial & archaeal annotations."""
    return parse('bacteria.csv.gz')


def virus():
    """Return a pandas dataframe of viral annotations."""
    return parse('viruses.csv.gz')


def eukaryote():
    """Return a pandas dataframe of eukaryotic annotations."""
    return parse('eukaryota.csv.gz')


def protein_group(kind):
    fname = final_table_path({
        'b': 'biocideResistanceProts.txt',
        'd': 'drugResistanceProts.txt',
        'm': 'mobilityProts.txt',
        'r': 'repairProts.txt',
        's': 'sporeProts.txt',
    }[kind[0].lower()])
    with open(fname) as f:
        return [line.strip() for line in f if line.strip()]
    assert False


def md1():
    """Return a pandas dataframe of the original microbe directory."""
    return parse('microbe-directory.csv.gz', ind=7)
