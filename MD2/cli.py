import click
import pandas as pd
import numpy as np
import csv
from .taxa_tree import NCBITaxaTree
from .clean_table import (
    reduce_col,
    reduce_row,
    rename_col,
    rename_MD1_tables,
    modify_dataset_value,
	clean_count_datasets,
	)
from .dataset_modification import (
    metasub_process,
    convert_taxa_tree,
	)
from .dataset_stats import (
    verify_column_names,
    dataset_stats,
    column_compare,
    )
	
@click.group()
def main():
    pass
	
@main.command('filter-microbes')
@click.argument('out', type=click.File('w'))
def filter_microbes(out):
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
    col_names = ['scientific_name', 'taxonomic_id', 'rank']
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

@main.command('update-bacteria')
@click.argument('file', type=click.File('r'))
@click.argument('out', type=click.File('w'))	
def update_bacteria(file, out):
    """Parse NCBI File to fill a column which exhibits certain hierarchical traits"""
    file_name = pd.read_csv(file, index_col=False)
    taxa_tree, sci_name = NCBITaxaTree.parse_files()
    rank_list = ['subspecies', 'species', 'species group', 'species subgroup', 'subgenus', 'genus', 'subfamily', 'family', 'suborder', 'order', 
		        'subclass', 'class', 'subphylum', 'phylum', 'kingdom', 'superkingdom', 'no rank', 'varietas', 'forma', 'tribe']
    spore_forming_genus = ['Bacillus', 'Clostridium', 'Sporolactobacillus', 'Sporosarcina', 'Cerasibacillus', 'Alkalibacillus', 'Amphibacillus',
                            'Anoxybacillus', 'Filobacillus', 'Geobacillus', 'Gracilibacillus', 'Halobacillus', 'Lentibacillus', 'Oceanobacillus',
                            'Paraliobacillus', 'Pontibacillus', 'Tenuibacillus', 'Thalassobacillus', 'Virgibacillus', 'Jeotgalibacillus', 
                            'Marinibacillus', 'Planomicrobium', 'Ureibacillus', 'Sporobacterium', 'Desulfitobacterium', 'Desulfonispora', 
                            'Desulfosporosinus', 'Desulfotomaculum', 'Sporotomaculum', 'Syntrophobotulus', 'Pelotomaculum', 'Sporotomaculum', 
                            'Syntrophobotulus', 'Thermincola', 'Filifactor', 'Tepidibacter', 'Anaerotruncus', 'Oscillospira', 'Sporobacter']		
    not_spore_forming_genus = ['Halolactibacillus', 'Marinococcus', 'Saccharococcus', 'Planococcus', 'Caryophanon', 'Filibacter', 'Kurthia',
                                'Peptococcus', 'Cryptanaerobacter', 'Dehalobacter', 'Anaerofilum', 'Acetivibrio', 'Acetanaerobacterium', 
                                'Fastidiosipila', 'Papillibacter', 'Subdoligranulum']
    for index, rows in file_name.iterrows():
        #Update Spore_Forming in a Top-Down approach based on Genus
        if rank_list.index(rows['rank']) <= rank_list.index('genus'):
            if taxa_tree.genus(rows['scientific_name'], default=None) in (spore_forming_genus):
                file_name.loc[file_name['scientific_name']==rows['scientific_name'], 'spore_forming'] = 'Always'
            elif taxa_tree.genus(rows['scientific_name'], default=None) in (not_spore_forming_genus):
                file_name.loc[file_name['scientific_name']==rows['scientific_name'], 'spore_forming'] = 'Never'
        #Update Gram_Stain in a Mixed approach based on Genus
        if taxa_tree.ancestors_list('genus', rows['scientific_name'], default=None) != None:
            ancestors = taxa_tree.ancestors_list('genus', rows['scientific_name'], default=None).copy()
            df_val = file_name[file_name['scientific_name'].isin(ancestors)]
            if 'Positive' in df_val['gram_stain'].values:
                file_name.loc[file_name['scientific_name'].isin(ancestors), 'gram_stain'] = 'Positive'
            elif 'Negative'  in df_val['gram_stain'].values:
                file_name.loc[file_name['scientific_name'].isin(ancestors), 'gram_stain'] = 'Negative'           
    file_name.to_csv(out)
 
@main.command('stats-file')
@click.option('--microbe', default='bacteria', help='Whether the file contains bacteria, virus or fungi')
@click.argument('file', type=click.File('r'))
@click.argument('out', type=click.File('w'))
def clean_file(microbe, file, out):
    """Outputs the statistics of column that needs to be filled""" 
    file_name = pd.read_csv(file)
    final_file, stats_file = verify_column_names(file_name, microbe, default=None)
    final_file.to_csv(out)
    stats_file.to_csv(str(microbe + '_stats.csv'))
	
@main.command('metasub-preprocessing')
@click.option('--feature-name', default='city', help='The feature to consider')
@click.option('--subtext', default='metasub', help='Name of Study')
@click.argument('file', type=click.File('r'))
@click.argument('metadata-file', type=click.File('r'))
@click.argument('out', type=click.File('w'))
def metasub_preprocess(feature_name, subtext, file, metadata_file, out):
    """Construct a table to integrate MetaSUB data based on chosen feature"""
    file_name = tbl = pd.read_csv(filename, index_col=0) 
    compiled_metasub = metasub_process(file_name, metadata_file, feature_name, subtext)
    compiled_metasub.to_csv(out)
	
@main.command('dataset-preprocessing')
@click.option('--feature-name', default='city', help='The feature to consider')
@click.option('--subtext', default='MetaSUB', help='Name of Study')
@click.argument('file', type=click.File('r'))
@click.argument('biom-file', type=click.File('r'))
@click.argument('metadata-file', type=click.File('r'))
@click.argument('out', type=click.File('w'))
def dataset_preprocess(feature_name, subtext, file, biom_file, metadata_file, out):
    """Construct a table to integrate other datasets based on chosen features"""
    compiled_metasub = convert_taxa_tree(file, biom_file, metadata_file, feature_name, subtext)
    compiled_metasub.to_csv(out)
    
@main.command('column-compare')
@click.argument('file1', type=click.File('r'))
@click.argument('file2', type=click.File('r'))
def dataset_column_compare(file1, file2):
    """Compare the before and after statistics of a file"""
    file_data_1 = pd.read_csv(file1)
    file_data_2 = pd.read_csv(file2)
    file_data_1 = file_data_1.replace(r'^\s*$', np.nan, regex=True)
    file_data_2 = file_data_2.replace(r'^\s*$', np.nan, regex=True)
    stats1, stats2 = column_compare(file_data_1, file_data_2)
    stats1.to_csv('column_stats_before.csv')
    stats2.to_csv('column_stats_after.csv')

if __name__ == '__main__':
    main()
