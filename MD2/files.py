import pandas as pd
import numpy as np
from scipy.stats import chisquare
from scipy import stats
from collections import Counter, defaultdict
import csv


def compare_datasets(taxa_list_1, taxa_list_2):
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
     

def count_values(values, value_being_compared):
    x=defaultdict(float)
    for var in [True, False]:
        x[var] = 1 / (1000 * 1000)
    for var in values:
        if var == value_being_compared:
            x[True] += 1
        else:
            x[False] +=1
    return x
        

def compare_categorical(value_being_compared, values_in_taxa_list_1, values_in_taxa_list_2):
    all_variables= set(values_in_taxa_list_1) | set(values_in_taxa_list_2)
    stats1 = count_values(values_in_taxa_list_1, value_being_compared)
    stats1 = pd.Series(stats1)
    stats2 = count_values(values_in_taxa_list_2, value_being_compared)
    stats2 = pd.Series(stats2)
    a = chisquare(stats1, stats2)
    return pd.Series({
        'abundance_in': stats1,
        'abundance_out': stats2,
        'p-value': a.pvalue,
    })


def compare_numeric(values_in_taxa_list_1, values_in_taxa_list_2):
    """Retun a Pandas Series with [abundance-in, abundance-out, p-value]."""
    mean1 = values_in_taxa_list_1.mean()
    mean2 = values_in_taxa_list_2.mean()
    a = stats.ttest_ind(values_in_taxa_list_1, values_in_taxa_list_2, equal_var=False)
    return pd.Series({
        'abundance_in': mean1, 
        'abundance_out': mean2,
        'p-value': a.pvalue,
    })


if __name__ == '__main__':
    # Run some simple tests
    categorical_test = compare_categorical(
        'yes',
         pd.Series(['yes', 'yes', 'no', 'yes', 'yes']),  # note that these are NOT the same length
         pd.Series(['no', 'no', 'no', 'yes', 'yes', 'no', 'no']),
    )
    print(categorical_test)

    cat_test_empty = compare_categorical(
        'A',
        pd.Series(['A', 'B', 'D']),
        pd.Series(['B', 'C', 'D', 'E'])
    )    
    print(cat_test_empty)
    
    numeric_test = compare_numeric(
        pd.Series([0, 1, 3, 0, 1, 1, 2, 2]),
        pd.Series([2, 2, 1, 3, 1, 3, 4]),
    )
    print(numeric_test)


df3 = pd.read_csv(r'MetaPhlan2_file2.csv')

df5 = pd.read_csv(r'MetaPhlan2_file1.csv')

MICROBE_DIRECTORY = pd.read_csv (r'microbe-directory.csv', index_col=7)


df = pd.DataFrame (MICROBE_DIRECTORY.iloc[0:5, 7:30])

df4 = pd.DataFrame (MICROBE_DIRECTORY.iloc[9:14, 7:30])

Taxa_List_1 = MICROBE_DIRECTORY.iloc[0:5].index.tolist()

Taxa_List_2 = MICROBE_DIRECTORY.iloc[9:14].index.tolist()

CATEGORICAL_LIST = ["gram_stain", "microbiome_location", "antimicrobial_susceptibility", "extreme_environment", "biofilm_forming", "animal_pathogen", "spore_forming", "plant_pathogen"]
NUMERICAL_LIST = ["optimal_temperature", "optimal_ph", "pathogenicity"]

def compare_microbe_directory_dataframes(values_in_taxa_list_1, values_in_taxa_list_2):
    df_final = pd.DataFrame (columns = ['variable', 'type', 'dataset', 'value', 'abundance_in', 'abundance_out', 'p-value'])
    for column_name in CATEGORICAL_LIST:
        a = values_in_taxa_list_1[column_name]
        b = values_in_taxa_list_2[column_name]
        c = (values_in_taxa_list_1[column_name] + values_in_taxa_list_1[column_name]).unique()
        for var in c:
            r = compare_categorical(var, values_in_taxa_list_1[column_name], values_in_taxa_list_2[column_name])
            df_final = df_final.append({'variable': column_name, 'type': 'categorical', 'dataset': 'df', 'value': var, 'abundance_in': r[0], 'abundance_out': r[1], 'p-value': r[2]}, ignore_index = True)
    for column_name in NUMERICAL_LIST:
        s = compare_numeric(values_in_taxa_list_1[column_name], values_in_taxa_list_2[column_name])
        df_final = df_final.append({'variable': column_name, 'type': 'numerical', 'dataset': 'df', 'value': 'mean', 'abundance_in': s[0], 'abundance_out': s[1], 'p-value': s[2]}, ignore_index = True)
    return df_final 

def compare_taxa_lists(values_in_taxa_list_1, values_in_taxa_list_2):
    df1= MICROBE_DIRECTORY.loc[values_in_taxa_list_1]
    df2 = MICROBE_DIRECTORY.loc[values_in_taxa_list_2]
    return compare_microbe_directory_dataframes(df1, df2)


if __name__ == '__main__':
    DataFrame_test = compare_microbe_directory_dataframes(
        df,
        df4,
    )
    print(DataFrame_test)

    Taxa_List_test = compare_taxa_lists(
        Taxa_List_1,
        Taxa_List_2,
    )
    print(Taxa_List_test)

