
import pandas as pd
import numpy as np

from scipy.stats import chisquare
from scipy import stats
from collections import defaultdict
from random import choices


def compare_categorical(value_being_compared, values_in_taxa_list_1, values_in_taxa_list_2):
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


def count_values(values, value_being_compared):
    x = defaultdict(float)
    for var in [True, False]:
        x[var] = 1 / (1000 * 1000)
    for var in values:
        if var == value_being_compared:
            x[True] += 1
        else:
            x[False] += 1
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
    tenthousand_samples1 = choices(keyslist1, valueslist1, k=10**4)
    tenthousand_samples2 = choices(keyslist2, valueslist2, k=10**4)
    a = stats.ks_2samp(tenthousand_samples1, tenthousand_samples2)
    return pd.Series({
        'abundance_in': mean1,
        'abundance_out': mean2,
        'p-value': a.pvalue,
    })


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
