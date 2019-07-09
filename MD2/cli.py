import click
import pandas as pd
from .taxa_tree import NCBITaxaTree

@click.group()
def main():
    pass
	
@main.command('data-table')
def data_table():
    pass
	
@main.command('filter-taxa')
def filter_taxa():
    pass
	
@main.command('annotate-taxa')
@click.argument('out', type=click.File('w'))
def annotate(out):
    """Parse NCBI File to return the taxonomic classification of a species"""
    taxa_tree, sci_name = NCBITaxaTree.parse_files()
    phyla = list()
    for taxon in sci_name:
        taxon = taxon.strip()  
        phyla.append([taxa_tree.phyla(taxon, default='unknown')])
    annotated = pd.DataFrame.from_dict({'taxa': sci_name, 'phyla': phyla}, orient='columns')
    annotated = annotated.set_index('taxa')
    annotated.to_csv(out)
	

if __name__ == '__main__':
    main()
