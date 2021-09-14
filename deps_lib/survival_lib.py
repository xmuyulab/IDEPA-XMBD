import copy
import os
import datetime
import shutil
import random
import configparser

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from scipy import stats
from deps_lib import utils, raw_data, stable_pairs, methods_comp

import rpy2
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter



def col_adapt(samples_col, columns, characters = [' ', '-'], new_char = '.'):
    """
    调整clinial中samples列中，列名的格式
    """
    for ch in characters:
        for idx in samples_col.index:
            samples_col.loc[idx, columns] = samples_col.loc[idx, columns].replace(ch, new_char)
    return samples_col



def cox_analysis(_dep, clinial_select, THRES_SURV, cox_run_workdir):
    
    if os.path.exists(cox_run_workdir) == False:
        os.makedirs(cox_run_workdir)
    os.chdir(cox_run_workdir)
    
    _specific_prot = _dep[(_dep.sum(axis=1) >= _dep.shape[1] * THRES_SURV[0]) & (_dep.sum(axis=1) <= _dep.shape[1] * THRES_SURV[1])].index.tolist()

    for _dp in _specific_prot:
        clinial_select[_dp] = np.nan

    for _dp in _specific_prot:
        clinial_select.loc[:,_dp] = _dep.loc[_dp, :].astype('int').tolist()

    _specific_prot_pd = pd.DataFrame({'specific_prot': _specific_prot})

    with localconverter(ro.default_converter + pandas2ri.converter):
        r_clinial_select = ro.conversion.py2rpy(clinial_select)

    with localconverter(ro.default_converter + pandas2ri.converter):
        r_specific_dep_prot = ro.conversion.py2rpy(_specific_prot_pd)

    rf_cox_analysis = rpy2.robjects.globalenv['cox_analysis']

    r_result = rf_cox_analysis(r_clinial_select, 
                                   r_specific_dep_prot,
                                   cox_run_workdir)

    with localconverter(ro.default_converter + pandas2ri.converter):
        cox_result = ro.conversion.rpy2py(rpy2.robjects.globalenv['result_select'])
    
    cox_rate = cox_result.shape[0] / len(_specific_prot)

    return cox_result, cox_rate