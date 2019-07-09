import pandas as pd
import numpy as np
from scipy.stats import chisquare


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
    pass


def compare_categorical(value_being_compared, values_in_taxa_list_1, values_in_taxa_list_2):
    """Return a Pandas Series with [abundance-in, abundance-out, p-value]."""
    pass


def compare_numeric(values_in_taxa_list_1, values_in_taxa_list_2):
    """Return a Pandas Series with [abundance-in, abundance-out, p-value]."""
    return pd.Series({
        'abundance_in': values_in_taxa_list_1.mean(),
        'abundance_out': values_in_taxa_list_2.mean(),
        'p-value': 1,  # TODO
    })


if __name__ == '__main__':
    # Run some simple tests
    categorical_test = compare_categorical(
        'yes',
         pd.Series(['yes', 'yes', 'no', 'yes', 'yes']),  # note that these are NOT the same length
         pd.Series(['no', 'no', 'no', 'yes', 'yes', 'no', 'no']),
    )
    print(categorical_test)
    numeric_test = compare_numeric(
        pd.Series([0, 1, 3, 0, 1, 1, 2, 2]),
        pd.Series([2, 2, , 3, 1, 3, 4]),
    )
    print(numeric_test)