import pandas as pd
import numpy as np
from scipy.stats import chisquare
from scipy import stats
from collections import Counter, defaultdict
from random import choices 
import csv

MICROBE_DIRECTORY = pd.read_csv(r'microbe-directory.csv', index_col=[7])
NUMERICAL_LIST = ["optimal_temperature", "optimal_ph", "pathogenicity"]
CATEGORICAL_LIST = ["gram_stain", "microbiome_location", "antimicrobial_susceptibility", "extreme_environment", "biofilm_forming", "animal_pathogen", "spore_forming", "plant_pathogen"]

def count_values_abundances(values, value_being_compared):
    x = defaultdict(float)
    for var in [True, False]:
        x[var] = 1 / (1000 * 1000)
    for var in values.keys():
        if var == value_being_compared:
            x[True] += values[var]
        else:
            x[False] += values[var]
    return x
        
def compare_categorical_abundances(value_being_compared, values_in_taxa_list_1, values_in_taxa_list_2):
    stats1 = count_values_abundances(values_in_taxa_list_1, value_being_compared)
    stats1 = pd.Series(stats1)
    stats2 = count_values_abundances(values_in_taxa_list_2, value_being_compared)
    stats2 = pd.Series(stats2)
    keyslist1 = list(values_in_taxa_list_1.keys())
    keyslist2 = list(values_in_taxa_list_2.keys())
    valueslist1 = list(values_in_taxa_list_1.values())
    valueslist2 = list(values_in_taxa_list_2.values())
    tenthousand_samples1 = choices(keyslist1, valueslist1, k = 10**4)
    tenthousand_samples2 = choices(keyslist2, valueslist2, k = 10**4)
    a = stats.ks_2samp(tenthousand_samples1, tenthousand_samples2)
    return pd.Series({
        'abundance_in': stats1,
        'abundance_out': stats2,
        'p-value': a.pvalue,
    })

def mean_ignore_nans(dictin):
    num = 0
    denom = 0
    for key, val in dictin.items():
        if not np.isnan(key):
            num += key*val
            denom += val
    return num/denom if denom != 0 else 0

def compare_numeric_abundances(values_in_taxa_list_1, values_in_taxa_list_2):
    """Retun a Pandas Series with [abundance-in, abundance-out, p-value]."""
    mean1 = mean_ignore_nans(values_in_taxa_list_1)
    mean2 = mean_ignore_nans(values_in_taxa_list_2)
    keyslist1 = list(values_in_taxa_list_1.keys())
    keyslist2 = list(values_in_taxa_list_2.keys())
    valueslist1 = list(values_in_taxa_list_1.values())
    valueslist2 = list(values_in_taxa_list_2.values())
    tenthousand_samples1 = choices(keyslist1, valueslist1, k = 10**4)
    tenthousand_samples2 = choices(keyslist2, valueslist2, k = 10**4)
    a = stats.ks_2samp(tenthousand_samples1, tenthousand_samples2)
    return pd.Series({
        'abundance_in': mean1, 
        'abundance_out': mean2,
        'p-value': a.pvalue,
    })

def compare_microbe_directory_dataframes(values_in_taxa_list_1, values_in_taxa_list_2):
    df_final = pd.DataFrame(columns = ['variable', 'type', 'dataset', 'value', 'abundance_in', 'abundance_out', 'p-value']) 
    for column_name in CATEGORICAL_LIST:
        taxa_list1 = values_in_taxa_list_1[column_name]
        taxa_list2 = values_in_taxa_list_2[column_name]
        values_list = (values_in_taxa_list_1[column_name] + values_in_taxa_list_1[column_name]).unique()
        dict1 = {} 
        dict2 = {}
        [dict1.update({values_in_taxa_list_1.at[key,column_name] : values_in_taxa_list_1.at[key,'WEIGHT']}) for key in values_in_taxa_list_1.index.tolist()]
        [dict2.update({values_in_taxa_list_2.at[key,column_name] : values_in_taxa_list_2.at[key,'WEIGHT']}) for key in values_in_taxa_list_2.index.tolist()]
        for var in values_list:
            categorical = compare_categorical_abundances(var, dict1, dict2)
            df_final = df_final.append({'variable': column_name, 'type': 'categorical', 'dataset': 'df', 'value': var, 'abundance_in': categorical[0], 'abundance_out': categorical[1], 'p-value': categorical[2]}, ignore_index = True)
    for column_name in NUMERICAL_LIST:
        dict1 = {}
        dict2 = {}
        [dict1.update({values_in_taxa_list_1.at[key,column_name] : values_in_taxa_list_1.at[key,'WEIGHT']}) for key in values_in_taxa_list_1.index.tolist()]
        [dict2.update({values_in_taxa_list_2.at[key,column_name] : values_in_taxa_list_2.at[key,'WEIGHT']}) for key in values_in_taxa_list_2.index.tolist()]
        numeric = compare_numeric_abundances(dict1, dict2)
        df_final = df_final.append({'variable': column_name, 'type': 'numerical', 'dataset': 'df', 'value': 'mean', 'abundance_in': numeric[0], 'abundance_out': numeric[1], 'p-value': numeric[2]}, ignore_index = True)
    return df_final

def compare_taxa_lists_abundances(values_in_taxa_list_1, values_in_taxa_list_2):
    """Return a Pandas DataFrame listing differences between two taxa lists.
    The DataFrame should have the following columns:
     1. variable, the thing being compared (e.g. spore forming)
     2. dataset, which dataset is being compared (e.g. taxa_list_1)
     3. value, the value being listed. For categorical variables this
        will be a value of the categorical variable (e.g. red from
        [red, green, blue]). For numeric variables this will be the
        name of a summary statistic (e.g. mean, median, mode)
     4. abundance-in, the value of (3) in (2) (e.g. taxa list 1)
     5. abundance-out, the value of (3) not in (2) (e.g. taxa list 2)
     6. p-value, the p-value of the difference
    This is a long format dataframe which means some of the data may
    be repeated.
    """

    df1 = MICROBE_DIRECTORY.loc[values_in_taxa_list_1.keys()]
    df1['WEIGHT'] = values_in_taxa_list_1.values 
    df2 = MICROBE_DIRECTORY.loc[values_in_taxa_list_2.keys()]
    df2['WEIGHT'] = values_in_taxa_list_2.values 
    return compare_microbe_directory_dataframes(df1, df2)

if __name__ == '__main__':
    # Run some simple tests
    cat_test_abundances = compare_categorical_abundances( 
            'A',
            {'A':0.2, 'B':0.3, 'D':0.5},
            {'B':0.25, 'C':0.4, 'D':0.25, 'E':0.1}
    )    
    print(cat_test_abundances)
    
    numeric_abundances_test = compare_numeric_abundances(
            {5:0.2, 6:0.25, 7:0.25, 8:0.1, 9:0.2},
            {4:0.2, 6:0.125, 7:0.3, 8:0.375},
    )
    print(numeric_abundances_test)
    
    taxa_lists_abundances_test = compare_taxa_lists_abundances(
            pd.Series({MICROBE_DIRECTORY.iloc[[0]].index[0]:0.2, MICROBE_DIRECTORY.iloc[[1]].index[0]:0.1, MICROBE_DIRECTORY.iloc[[2]].index[0]:0.3,   MICROBE_DIRECTORY.iloc[[3]].index[0]:0.25, MICROBE_DIRECTORY.iloc[[4]].index[0]:0.15}),
            pd.Series({MICROBE_DIRECTORY.iloc[[9]].index[0]:0.25, MICROBE_DIRECTORY.iloc[[10]].index[0]:0.3, MICROBE_DIRECTORY.iloc[[11]].index[0]:0.15, MICROBE_DIRECTORY.iloc[[12]].index[0]:0.2, MICROBE_DIRECTORY.iloc[[13]].index[0]:0.1})
    )
    print(taxa_lists_abundances_test)
