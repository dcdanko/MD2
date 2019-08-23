import pandas as pd
import numpy as np
from scipy.stats import chisquare
from scipy import stats
from collections import Counter, defaultdict
from random import choices 
import csv

def count_values_abundances(values, value_being_compared):
    x = defaultdict(float)
    for var in [True, False]:
        x[var] = 1 / (1000 * 1000)
  #  norm = /len(values)
    for var in values.keys():
        if var == value_being_compared:
            x[True] += values[var]
        else:
            x[False] += values[var]
    return x
        

def compare_categorical_abundances(value_being_compared, values_in_taxa_list_1, values_in_taxa_list_2):
    all_variables = set(values_in_taxa_list_1) | set(values_in_taxa_list_2)
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

def compare_numeric_abundances(values_in_taxa_list_1, values_in_taxa_list_2):
    """Retun a Pandas Series with [abundance-in, abundance-out, p-value]."""
    mean1 = sum([key * val for key, val in values_in_taxa_list_1.items()])
    mean2 = sum([key * val for key, val in values_in_taxa_list_2.items()])
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
        for var in values_list:
            categorical = compare_categorical(var, values_in_taxa_list_1[column_name], values_in_taxa_list_2[column_name])
            df_final = df_final.append({'variable': column_name, 'type': 'categorical', 'dataset': 'df', 'value': var, 'abundance_in': categorical[0], 'abundance_out': categorical[1], 'p-value': categorical[2]}, ignore_index = True)
    for column_name in NUMERICAL_LIST:
        numeric = compare_numeric(values_in_taxa_list_1[column_name], values_in_taxa_list_2[column_name])
        df_final = df_final.append({'variable': column_name, 'type': 'numerical', 'dataset': 'df', 'value': 'mean', 'abundance_in': numeric[0], 'abundance_out': numeric[1], 'p-value': numeric[2]}, ignore_index = True)
    return df_final


def compare_taxa_lists(values_in_taxa_list_1, values_in_taxa_list_2):
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
    df1 = MICROBE_DIRECTORY.loc[values_in_taxa_list_1]
    df2 = MICROBE_DIRECTORY.loc[values_in_taxa_list_2]
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

