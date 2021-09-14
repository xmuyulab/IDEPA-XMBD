import numpy as np
import pandas as pd

from scipy import stats
from deps_lib import utils
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector

r_stats = importr('stats')


def get_group_deps(individual_deps, THRES):
    
    n = individual_deps.shape[1]
    p = (individual_deps != 0).sum().sum() / (individual_deps.shape[0] * individual_deps.shape[1])

    p_list = []

    for idx in individual_deps.index:
        x = (individual_deps.loc[idx, :] != 0).sum()
        p_list.append(np.sum(stats.binom.pmf(np.arange(x,n+1), n, p)))
    q_list = list(r_stats.p_adjust(FloatVector(p_list), method = 'BH'))
    _label = np.array(q_list) < THRES
    group_deps = individual_deps[_label].index.tolist()
    
    return group_deps