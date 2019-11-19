
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
@click.option('--feature-name', default='city', help='The feature to consider')
@click.option('--subtext', default='metasub', help='Name of Study')
@click.argument('file', type=click.File('r'))
@click.argument('metadata-file', type=click.File('r'))
@click.argument('out', type=click.File('w'))
def metasub_preprocess(feature_name, subtext, file, metadata_file, out):
    """Construct a table to integrate MetaSUB data based on chosen feature"""
    file_name = pd.read_csv(file, index_col=0)
    compiled_metasub = metasub_process(file_name, metadata_file, feature_name, subtext)
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
