import click
import pandas as pd
import csv
from .taxa_tree import NCBITaxaTree


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
    merged_file = taxa_tree.data_table(file_path)
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
    annotated = pd.DataFrame.from_dict(bacteria, columns=col_names, orient='index')
    annotated = annotated.set_index('scientific name')
    annotated = annotated.join(merged_file, how='left')
    annotated.to_csv("NCBI_Bacteria_rank.csv")
    annotated = pd.DataFrame.from_dict(viruses, columns=col_names, orient='index')
    annotated = annotated.join(merged_file, how='left')
    annotated.to_csv("NCBI_Virus_rank.csv")
    annotated = pd.DataFrame.from_dict(fungi, columns=col_names, orient='index')
    annotated = annotated.join(merged_file, how='left')
    annotated.to_csv("NCBI_Fungi_rank.csv")
	

if __name__ == '__main__':
    main()
