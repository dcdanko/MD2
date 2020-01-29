import click
import pandas as pd

from microbe_directory.taxa_tree import NCBITaxaTree, TaxonomicRankError
from microbe_directory.clean_table import file_clean, clean_columns
from microbe_directory.dataset_modification import taxa_to_organism
from microbe_directory.infill_fields import infill_bacterial_fields
from microbe_directory.constants import DOMAINS, FUNGI

from .cli_preprocessing import preprocessing
from .cli_stats import stats


@click.group()
def main():
    pass


main.add_command(stats)
main.add_command(preprocessing)


@main.command('all-taxa-table')
@click.option('-o', '--outfile', default='-', type=click.File('w'))
def filter_microbes(outfile):
    """Write a table of all taxa and their taxonomy to outfile."""
    taxa_tree = NCBITaxaTree.parse_files()
    taxonomy = {}
    for taxon in taxa_tree.all_names():
        try:
            taxonomy[taxon] = taxa_tree.canonical_taxonomy(taxon)
        except TaxonomicRankError:
            pass
    annotated = pd.DataFrame.from_dict(taxonomy, orient='index')
    annotated.to_csv(outfile)


@main.command('taxa-table')
@click.option('-o', '--outfile', default='-', type=click.File('w'))
@click.option('-s', '--superkingdom', default='bacteria', type=click.Choice(DOMAINS),
              help='Biological domain of the table')
def annotate(outfile, superkingdom):
    """Make a master CSV table of taxa."""
    taxa_tree = NCBITaxaTree.parse_files()
    out_table = {}
    for taxon in taxa_tree.all_names():
        try:
            taxon_id, rank, taxon_superkingdom = taxa_tree.place_microbe(taxon)
            if taxon_superkingdom.lower() == superkingdom:
                if superkingdom != FUNGI.lower() or 'Fungi' in taxa_tree.ancestors(taxon):
                    out_table[taxon] = (taxon, taxon_id, rank)
        except TaxonomicRankError:
            pass
    annotated = pd.DataFrame.from_dict(
        out_table, columns=['scientific_name', 'taxonomic_id', 'rank'], orient='index'
    )
    annotated = annotated.dropna()
    annotated.reset_index(drop=True)
    annotated.to_csv(outfile)


def __find_name_column(df, scientific_names):
    """Return a tuple of (name_column, overla_size)."""
    max_overlap, name_col_name = -1, None
    for col_name in df.columns:
        unique_col_vals = set(df[col_name])
        overlap_size = len(unique_col_vals & scientific_names)
        if overlap_size > max_overlap:
            max_overlap = overlap_size
            name_col_name = col_name
    return name_col_name, max_overlap


@main.command('merge-csvs')
@click.option('-o', '--outfile', default='-', type=click.File('w'))
@click.argument('master_table', type=click.File('r'))
@click.argument('csv_files', nargs=-1, type=click.File('r'))
def merge_csv_files(outfile, master_table, csv_files):
    master_table = pd.read_csv(master_table)
    with click.progressbar(csv_files) as pbar:
        for csv_file in pbar:
            df = pd.read_csv(csv_file)
            df.columns = [col.lower() for col in df.columns]
            if ('species' in df.columns or 'genus' in df.columns) and 'class' in df.columns:
                df = taxa_to_organism(df)
            col_name, overlap_size = __find_name_column(df, set(master_table['scientific_name']))
            if overlap_size > 0:
                master_table = master_table.merge(
                    df, left_on='scientific_name', right_on=col_name, how='left'
                )
    master_table.to_csv(outfile)


@main.command('clean-file')
@click.option('-o', '--outfile', default='-', type=click.File('w'))
@click.argument('file', type=click.File('r'))
def clean_file(outfile, file):
    """Clean up data-table to be used for Microbe Directory 2.0 and above"""
    tbl = pd.read_csv(file, index_col=0)
    cleaned_file = file_clean(tbl)
    cleaned_file = cleaned_file.replace('nan', '')
    cleaned_file.reset_index()
    cleaned_file = clean_columns(cleaned_file)
    cleaned_file.to_csv(outfile)


@main.command('infill-bacteria')
@click.option('-o', '--outfile', default='-', type=click.File('w'))
@click.argument('table', type=click.File('r'))
def update_bacteria(outfile, table):
    """Parse NCBI File to fill a column which exhibits certain hierarchical traits"""
    table = pd.read_csv(table, index_col=1, dtype=str)
    table = infill_bacterial_fields(table)
    table.to_csv(outfile)


if __name__ == '__main__':
    main()
