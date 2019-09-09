import pandas as pd
import numpy as np


BACTERIA = ['gram_stain', 'antimicrobial_susceptibility', 'biofilm_forming', 'extremophile', 'extreme_environment', 'genome_availability', 'pathogenicity', 
            'animal_pathogen', 'plant_pathogen', 'fungi_pathogen', 'human_pathogen', 'human_microbiome_body_site', 'spore_forming', 'optimal_temperature', 
            'optimal_ph', 'virus_lineage', 'virus_name', 'commensal_microbe', 'commensal_host', 'metabolism_energy', 'metabolism_carbon', 
            'metabolism_reducing_equivalent', 'soil_location', 'extreme_environment_location', 'water_location']
			
FUNGI = ['gram_stain', 'antimicrobial_susceptibility', 'biofilm_forming', 'extremophile', 'extreme_environment', 'genome_availability', 'pathogenicity', 
        'animal_pathogen', 'plant_pathogen', 'fungi_pathogen', 'human_pathogen', 'human_microbiome_body_site', 'spore_forming', 'optimal_temperature', 
        'optimal_ph', 'virus_lineage', 'virus_name', 'commensal_microbe', 'commensal_host', 'metabolism_energy', 'metabolism_carbon', 
        'metabolism_reducing_equivalent', 'soil_location', 'extreme_environment_location', 'water_location']
		
VIRUS = ['antimicrobial_susceptibility', 'disease', 'extremophile', 'extreme_environment', 'genome_availability', 'genetic_material', 'strand, capsid', 
        'capsid_symmetry', 'host_lineage', 'host_name', 'pathogenicity', 'animal_pathogen', 'plant_pathogen', 'fungi_pathogen', 'bacteria_pathogen', 
        'protozoa_pathogen', 'human_pathogen', 'human_microbiome_body_site', 'optimal_temperature', 'optimal_ph', 'extreme_environment_location', 
        'water_location']
	
def verify_column_names(file, microbe_type, default=None):
    """Create columns if they do not exits"""
    if microbe_type=='bacteria':
        for bact in BACTERIA: 
            if bact not in file.columns:
                file[bact] = np.nan
        statistics = dataset_stats(file, microbe_type)
        return file, statistics
    elif microbe_type=='virus':
        for vir in VIRUS: 
            if vir not in file.columns:
                file[vir] = np.nan
        statistics = dataset_stats(file, microbe_type)
        return file, statistics
    elif microbe_type=='fungi':
        for fung in FUNGI: 
            if fung not in file.columns:
                file[fung] = np.nan
        statistics = dataset_stats(file, microbe_type)
        return file, statistics
    return default, default
	
def dataset_stats(file, microbe_type):
    """Calculate the number of non empty columns for each microbe"""
    if microbe_type=='bacteria':
        new_file = file.filter(BACTERIA, axis=1)
    elif microbe_type=='virus':
        new_file = file.filter(VIRUS, axis=1)
    elif microbe_type=='fungi':
        new_file = file.filter(FUNGI, axis=1)
    non_null = new_file.count(axis = 1) 
    stats_file = pd.DataFrame(list(zip(file['scientific_name'], non_null)), columns=['scientific_name', 'counts']) 
    return stats_file
    
def column_compare(file1, file2):
    """Before and after comparisons of non-na values for column"""
    stats1 = file1.notna().sum()
    stats2 = file2.notna().sum()
    return stats1, stats2
