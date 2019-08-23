import click
import pandas as pd
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
    dataset_stats
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
        if taxa_tree.taxonomic_tree_species(taxon, default=None) != None:
            taxonomy[taxon] = taxa_tree.taxonomic_tree_species(taxon, default=None).values()
    col_names = ['taxonomic_id', 'species', 'genus', 'family', 'order', 'class', 'phylum']
    annotated = pd.DataFrame.from_dict(taxonomy, columns=col_names, orient='index')
    annotated.to_csv(out)
	
@main.command('annotate-taxa')
@click.argument('file_path', type=click.Path(exists=True))
def annotate(file_path):
    """Parse NCBI File to return the rank for given scientific name along with merged tables for all datasets"""
    taxa_tree, sci_name = NCBITaxaTree.parse_files()
    rank_file, tax_rank, rank = {}, list(), ''
    bacteria, viruses, fungi = {}, {}, {}
    for taxon in sci_name:
        taxon = taxon.strip()
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
    cleaned_file = clean_file(True, annotated)	
    final_file, statistics = verify_column_names(cleaned_file, 'virus')
    final_file.to_csv("NCBI_Virus_rank.csv")	
    statistics.to_csv("Virus_statistics.csv")	       
    annotated = pd.DataFrame.from_dict(bacteria, columns=col_names, orient='index')
    annotated = taxa_tree.data_table(file_path, annotated)		
    add_new_datas = filter_taxa(annotated)			
    cleaned_file = clean_file(False, add_new_datas)	
    final_file, statistics = verify_column_names(cleaned_file, 'bacteria')
    final_file.to_csv("NCBI_Bacteria_rank.csv")	
    statistics.to_csv("Bacteria_statistics.csv")					       
    annotated = pd.DataFrame.from_dict(fungi, columns=col_names, orient='index')
    annotated = taxa_tree.data_table(file_path, annotated)
    cleaned_file = clean_file(False, annotated)	
    final_file, statistics = verify_column_names(cleaned_file, 'fungi')
    final_file.to_csv("NCBI_Fungi_rank.csv")	
    statistics.to_csv("Fungi_statistics.csv")
                    
	
@main.command('metasub-preprocessing')
@click.option('--feature-name', default='city', help='The feature to consider')
@click.option('--subtext', default='metasub', help='Name of Study')
@click.argument('file', type=click.File('r'))
@click.argument('metadata-file', type=click.File('r'))
@click.argument('out', type=click.File('w'))
def metasub_preprocess(feature_name, subtext, file, metadata_file, out):
    """Construct a table to integrate MetaSUB data based on chosen feature"""
    file_name = pd.read_csv(filename, index_col=0) 
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
	
def clean_file(isvirus, file_name):
    """Clean up data-table to be used for Microbe Directory 2.0 and above"""
    header = list(file_name.columns.values)
    remove_col_file = reduce_col(isvirus, file_name)
    remove_row_file = reduce_row(isvirus, remove_col_file) 
    remove_row_file = remove_row_file.replace('nan', '')
    return remove_col_file
	
def filter_taxa(file_name):
    """Parse NCBI File to fill a column which exhibits certain hierarchical traits"""
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
    return file_name
	

if __name__ == '__main__':
    main()
