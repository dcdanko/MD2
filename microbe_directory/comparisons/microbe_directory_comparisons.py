
import pandas as pd

from .statistics import (
    compare_categorical,
    compare_numeric,
    compare_categorical_abundances,
    compare_numeric_abundances,
)
from .constants import (
    MICROBE_DIRECTORY,
    NUMERICAL_LIST,
    CATEGORICAL_LIST,
)


def compare_microbe_directory_dataframes(values_in_taxa_list_1, values_in_taxa_list_2):
    df_final = pd.DataFrame(columns = ['variable', 'type', 'dataset', 'value', 'abundance_in', 'abundance_out', 'p-value'])
    for column_name in CATEGORICAL_LIST:
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


def compare_microbe_directory_dataframes_abundances(values_in_taxa_list_1, values_in_taxa_list_2):
    df_final = pd.DataFrame(columns = ['variable', 'type', 'dataset', 'value', 'abundance_in', 'abundance_out', 'p-value']) 
    for column_name in CATEGORICAL_LIST:
        values_list = (values_in_taxa_list_1[column_name] + values_in_taxa_list_1[column_name]).unique()
        dict1 = {} 
        dict2 = {}
        [dict1.update({values_in_taxa_list_1.at[key, column_name]: values_in_taxa_list_1.at[key, 'WEIGHT']}) for key in values_in_taxa_list_1.index.tolist()]
        [dict2.update({values_in_taxa_list_2.at[key, column_name]: values_in_taxa_list_2.at[key, 'WEIGHT']}) for key in values_in_taxa_list_2.index.tolist()]
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
    return compare_microbe_directory_dataframes_abundances(df1, df2)
