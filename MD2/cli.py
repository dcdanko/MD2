import click
import pandas as pd
from .taxa_tree import NCBITaxaTree


@click.group()
def main():
    pass
	
@main.command('data-table')
def data_table():
    """Merge existing data tables"""
    pass
	
@main.command('filter-microbes')
@click.argument('out', type=click.File('w'))
def filter_taxa(out):
    """Parse NCBI File to return the phyla classification of a species"""
    taxa_tree, sci_name = NCBITaxaTree.parse_files()
    taxonomy = {}
    for taxon in sci_name:
        taxon = taxon.strip()
        if taxa_tree.taxonomic_rank(taxon, default=None) != None:
            taxonomy[taxon] = taxa_tree.taxonomic_rank(taxon, default=None).values()
    col_names = ['taxonomic_id', 'species', 'genus', 'family', 'order', 'class', 'phylum']
    annotated = pd.DataFrame.from_dict(taxonomy, columns=col_names, orient='index')
    annotated.to_csv(out)
	
@main.command('annotate-taxa')
@click.option('--microbes', default=True, help='Whether limit to microbes i.e., Bacteria and Virus, only')
@click.argument('out', type=click.File('w'))
def annotate(microbes, out):
    """Parse NCBI File to return the rank for given scientific name"""
    taxa_tree, sci_name = NCBITaxaTree.parse_files()
    rank_file = {}
    for taxon in sci_name:
        taxon = taxon.strip()
        if microbes == False:        
            rank_file[taxon] = taxa_tree.rank_of_species(taxon).values()
        else:
            if taxa_tree.rank_microbes(taxon, default=None) != None:
                rank_file[taxon] = taxa_tree.rank_microbes(taxon, default=None).values()
    col_names = ['scientific name', 'taxonomic_id', 'rank']
    annotated = pd.DataFrame.from_dict(rank_file, columns=col_names, orient='index')
    annotated.to_csv(out)
	

if __name__ == '__main__':
    main()
