import numpy as np
from xarray import Dataset

def round_sig_figs(x: float, sig: int=2) -> int | float:

    '''
    Rounds float <x> to <sig> significant figures.
    '''

    return round(x, sig-int(np.floor(np.log10(abs(x))))-1)

def forceNamingConvention(ds: Dataset) -> Dataset:

    '''
    Handles cases with Datasets where variables are labelled using multiple possible
    common names, e.g. 'time', 't', 'year', etc.. Currently only does this for time
    and basin but could be expanded to more variables and their common name variants.

    inputs:
        - ds: Dataset with any old variable names
    output:
        - ds with variables (hopefully) renamed to match standard naming conventions
    '''

    LookupVariants = {
        'time':     {'t', 'Time', 'year', 'Year'},
        'basin':    {'basins', 'Basin', 'Basins'}
    }

    for standard_name, variants in LookupVariants.items():
        wrong_vars = variants.intersection(ds.variables)
        if len(wrong_vars) > 0 and standard_name not in ds.variables:
            wrong_name = wrong_vars.pop()
            ds = ds.rename({wrong_name: standard_name})
    
    return ds