import click
import pandas as pd
import csv
from .taxa_tree import NCBITaxaTree
from .clean_table import (
    reduce_col,
    reduce_row,
    rename_col,
    rename_MD1_tables,
    metasub_process,
	)

@click.group()
def main():
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
            taxonomy[taxon] = taxa_tree.taxonomic_tree_species(taxon, default=None).values()
    col_names = ['taxonomic_id', 'species', 'genus', 'family', 'order', 'class', 'phylum']
    annotated = pd.DataFrame.from_dict(taxonomy, columns=col_names, orient='index')
    annotated.to_csv(out)
	
@main.command('annotate-taxa')
@click.option('--microbes', default=True, help='Whether limit to microbes i.e., Bacteria and Virus, only')
@click.argument('file_path', type=click.Path(exists=True))
def annotate(microbes, file_path):
    """Parse NCBI File to return the rank for given scientific name along with merged tables for all datasets"""
    taxa_tree, sci_name = NCBITaxaTree.parse_files()
    rank_file, tax_rank, rank = {}, list(), ''
    bacteria, viruses, fungi = {}, {}, {}
    for taxon in sci_name:
        taxon = taxon.strip()
        if microbes == False:        
            rank_file[taxon] = taxa_tree.rank_of_species(taxon).values()
        else:
            tax_rank, rank = taxa_tree.rank_microbes(taxon, default=None)
            if rank != None: 		
                if rank == 'Viruses':
                    viruses[taxon] = tax_rank
                elif rank == 'Bacteria':
                    bacteria[taxon] = tax_rank					
                elif rank == 'Fungi':
                    fungi[taxon] = tax_rank							
    col_names = ['scientific name', 'taxonomic_id', 'rank']
    annotated = pd.DataFrame.from_dict(viruses, columns=col_names, orient='index')
    annotated = taxa_tree.data_table(file_path, annotated)
    annotated.to_csv("NCBI_Virus_rank.csv")
    annotated = pd.DataFrame.from_dict(bacteria, columns=col_names, orient='index')
    annotate = taxa_tree.data_table(file_path, annotated)
    annotate.to_csv("NCBI_Bacteria_rank.csv")
    annotated = pd.DataFrame.from_dict(fungi, columns=col_names, orient='index')
    annotated = taxa_tree.data_table(file_path, annotated)
    annotated.to_csv("NCBI_Fungi_rank.csv")
	
@main.command('clean-file')
@click.option('--isvirus', default=False, help='Whether the file contains virus')
@click.argument('file', type=click.File('r'))
@click.argument('out', type=click.File('w'))
def clean_file(isvirus, file, out):
    """Clean up data-table to be used for Microbe Directory 2.0 and above"""
    file_name = pd.read_csv(file, index_col=False)
    header = list(file_name.columns.values)
    remove_col_file = reduce_col(isvirus, file_name)  
    remove_row_file = reduce_row(isvirus, remove_col_file)
    remove_row_file = remove_row_file.replace('nan', '')   
    remove_row_file.to_csv(out)
    
	
@main.command('metasub-preprocessing')
@click.option('--feature-name', default='city', help='The feature to consider')
@click.argument('file', type=click.File('r'))
@click.argument('metadata-file', type=click.File('r'))
@click.argument('out', type=click.File('w'))
def metasub_preprocess(feature_name, file, metadata_file, out):
    """Construct a table to integrate MetaSUB data based on chosen feature"""
    compiled_metasub = metasub_process(file, metadata_file, feature_name)
    compiled_metasub.to_csv(out)
    

if __name__ == '__main__':
    main()
