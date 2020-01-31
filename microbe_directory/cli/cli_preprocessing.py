
import click
import pandas as pd

from microbe_directory.dataset_modification import (
    metasub_process,
    convert_taxa_tree,
)


@click.group('preprocess')
def preprocessing():
    pass


@preprocessing.command('metasub')
@click.option('-f', '--feature-name', default='city', help='The feature to condense')
@click.option('-p', '--prefix', default='metasub_', help='Name of Study')
@click.option('-o', '--out', type=click.File('w'), default='-')
@click.argument('taxa_tbl', type=click.File('r'))
@click.argument('metadata_tbl', type=click.File('r'))
def metasub_preprocess(feature_name, prefix, out, taxa_tbl, metadata_tbl):
    """Construct a table to integrate MetaSUB data based on chosen feature"""
    taxa_tbl = pd.read_csv(taxa_tbl, index_col=0)
    metadata = pd.read_csv(metadata_tbl, index_col=0)
    compiled_metasub = metasub_process(taxa_tbl, metadata, feature_name, prefix)
    compiled_metasub.to_csv(out)


@preprocessing.command('dataset')
@click.option('--feature-name', default='city', help='The feature to consider')
@click.option('--subtext', default='MetaSUB', help='Name of Study')
@click.argument('file', type=click.File('r'))
@click.argument('biom-file', type=click.File('r'))
@click.argument('metadata-file', type=click.File('r'))
@click.argument('out', type=click.File('w'))
def dataset_preprocess(feature_name, subtext, file, biom_file, metadata_file, out):
    """Construct a table to integrate other datasets based on chosen features"""
    compiled_dataset = convert_taxa_tree(file, biom_file, metadata_file, feature_name)
    compiled_dataset = compiled_dataset.add_prefix(subtext)
    compiled_dataset.to_csv(out)
