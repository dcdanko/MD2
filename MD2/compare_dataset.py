import pandas as pd
import numpy as np
from scipy.stats import chisquare
from scipy import stats
from collections import Counter, defaultdict

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
     

def count_values(values):
    variables = pd.Series(values).unique()
    x=defaultdict(float)
    for val in values:
        x[val] += 1/len(values)
    """
    for element in range(0,len(variables)):
        for i in values:
            if i == variables[element]:
                x[i] += 1./len(values)    
    """
    return x
        

    #value_being_compared = variables[1]
    #x = 0 
    #for element in values:
     #   if element == value_being_compared:
       #     x+=1
   # x = x/len(values)
    #return x,1-x

def compare_categorical(value_being_compared, values_in_taxa_list_1, values_in_taxa_list_2):
    """Return a Pandas Series with [abundance-in, abundance-out, p-value]."""
    stats1 = count_values(values_in_taxa_list_1)
    stats2 = count_values(values_in_taxa_list_2)

   # values_in_taxa_list_1.to_frame(name=values_in_taxa_list_1)
   # values_in_taxa_list_2.to_frame(name=values_in_taxa_list_2)
    return pd.Series({
        'abundance_in': stats1,
        'abundance_out': stats2,  # TODO
    })


def compare_numeric(values_in_taxa_list_1, values_in_taxa_list_2):
    """Retun a Pandas Series with [abundance-in, abundance-out, p-value]."""
#    len1 = len(values_in_taxa_list_1)
#    len2 = len(values_in_taxa_list_2)
    mean1 = values_in_taxa_list_1.mean()
    mean2 = values_in_taxa_list_2.mean()
#    var_1 = values_in_taxa_list_1.var(ddof=1)
#    var_2 = values_in_taxa_list_2.var(ddof=1)
#    meantot = (len1*mean1+len2*mean2)/(len1+len2)
#    s = (len1*(var_1+(mean1-meantot)**2)+len2*(var_2+(mean2-meantot)**2))/(len1+len2)
    """
    s = var_1/len1 + var_2/len2
    t = (mean1 - mean2)/np.sqrt(s)
    denom = (var_1/len1)**2/(len1-1) + (var_2/len2)**2/(len2-1)
    df = s*s/denom
    p = 1 - stats.t.cdf(t,df=df)
    """
    a = stats.ttest_ind(values_in_taxa_list_1, values_in_taxa_list_2, equal_var=False)
    return pd.Series({
        'abundance_in': mean1, 
        'abundance_out': mean2,
        'p-value': a.pvalue,  # TODO
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
        pd.Series([2, 2, 1, 3, 1, 3, 4]),
    )
    print(numeric_test)
