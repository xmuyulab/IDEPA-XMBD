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


import rpy2
from rpy2.robjects.vectors import FloatVector
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

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

def get_kegg_result(up_path, down_path, GROUP_THRES, workdir, methods, R_dir, r_kegg_path):
    up = pd.read_csv(up_path, index_col=0).astype('int')
    down = pd.read_csv(down_path, index_col=0).astype('int')

    up[up == 1] = 2


    result = up + down

    group = get_group_deps(individual_deps = result, THRES = GROUP_THRES)

    _workdir = workdir + '/' + '/' + methods
    if os.path.exists(_workdir) == False:
        os.makedirs(_workdir)
    os.chdir(_workdir)

    group_pd = pd.DataFrame({'SYMBOL': group})
    group_pd.to_csv('./group_protein.txt', sep='\t', index=False, header=False)

    os.system('%s/Rscript %s %s %s'%(R_dir, r_kegg_path, _workdir+'/group_protein.txt', _workdir+'/kegg_result.csv'))

    kegg_result = pd.read_csv(_workdir+'/kegg_result.csv', index_col=0)
    kegg_result['methods'] = methods
    return kegg_result


def get_kegg_result_dataset(workdir, r_kegg_path, R_dir, GROUP_THRES, result_dir, kegg_methods_list, Q_VALUE_THRES):

    if os.path.exists(workdir) == False:
        os.makedirs(workdir)
    os.chdir(workdir)
    merge_data = []
    if 'Penda' in kegg_methods_list:
        penda_up_path = result_dir + '/penda/' + '/penda_up.csv'
        penda_down_path = result_dir + '/penda/' + '/penda_down.csv'
        
        kegg_result_penda = get_kegg_result(up_path = penda_up_path, 
                                            down_path = penda_down_path,
                                            GROUP_THRES = GROUP_THRES, 
                                            workdir = workdir,
                                            methods = 'penda',
                                            R_dir = R_dir, 
                                            r_kegg_path = r_kegg_path)
        kegg_result_penda = kegg_result_penda.loc[kegg_result_penda.sort_values(by='qvalue')[:Q_VALUE_THRES].index,:].reset_index(drop=True)
        if kegg_result_penda.shape[0] > 0 :
            merge_data.append(kegg_result_penda)

    if 'Peng method' in kegg_methods_list:
        peng_up_path = result_dir + '/peng/'  + '/peng_up.csv'
        peng_down_path = result_dir + '/peng/' + '/peng_down.csv'
        kegg_result_peng = get_kegg_result(up_path = peng_up_path, 
                                            down_path = peng_down_path,
                                            GROUP_THRES = GROUP_THRES, 
                                            workdir = workdir,
                                            methods = 'peng',
                                            R_dir = R_dir, 
                                            r_kegg_path = r_kegg_path)
        kegg_result_peng = kegg_result_peng.loc[kegg_result_peng.sort_values(by='qvalue')[:Q_VALUE_THRES].index,:].reset_index(drop=True)
        if kegg_result_peng.shape[0] > 0 :
            merge_data.append(kegg_result_peng)
        

    if 'Quantile' in kegg_methods_list:
        quantile_up_path = result_dir + '/quantile/' + '/quantile_up.csv'
        quantile_down_path = result_dir + '/quantile/' + '/quantile_down.csv'
        kegg_result_quantile = get_kegg_result(up_path = quantile_up_path, 
                                            down_path = quantile_down_path,
                                            GROUP_THRES = GROUP_THRES, 
                                            workdir = workdir,
                                            methods = 'quantile',
                                            R_dir = R_dir, 
                                            r_kegg_path = r_kegg_path)
        kegg_result_quantile = kegg_result_quantile.loc[kegg_result_quantile.sort_values(by='qvalue')[:Q_VALUE_THRES].index,:].reset_index(drop=True)
        if kegg_result_quantile.shape[0] > 0 :
            merge_data.append(kegg_result_quantile)

    if 'RankComp' in kegg_methods_list:
        rankc_v1_up_path = result_dir + '/rankcomp/'  + '/rankc_v1_up.csv'
        rankc_v1_down_path = result_dir + '/rankcomp/' + '/rankc_v1_down.csv'
        kegg_result_rankc_v1 = get_kegg_result(up_path = rankc_v1_up_path, 
                                            down_path = rankc_v1_down_path,
                                            GROUP_THRES = GROUP_THRES, 
                                            workdir = workdir,
                                            methods = 'rankc_v1',
                                            R_dir = R_dir, 
                                            r_kegg_path = r_kegg_path)
        kegg_result_rankc_v1 = kegg_result_rankc_v1.loc[kegg_result_rankc_v1.sort_values(by='qvalue')[:Q_VALUE_THRES].index,:].reset_index(drop=True)
        if kegg_result_rankc_v1.shape[0] > 0 :
            merge_data.append(kegg_result_rankc_v1)

        rankc_v2_up_path = result_dir + '/rankcomp/'  + '/rankc_v2_up.csv'
        rankc_v2_down_path = result_dir + '/rankcomp/' + '/rankc_v2_down.csv'
        kegg_result_rankc_v2 = get_kegg_result(up_path = rankc_v2_up_path, 
                                                down_path = rankc_v2_down_path,
                                                GROUP_THRES = GROUP_THRES, 
                                                workdir = workdir,
                                                methods = 'rankc_v2',
                                                R_dir = R_dir, 
                                                r_kegg_path = r_kegg_path)
        kegg_result_rankc_v2 = kegg_result_rankc_v2.loc[kegg_result_rankc_v2.sort_values(by='qvalue')[:Q_VALUE_THRES].index,:].reset_index(drop=True)
        if kegg_result_rankc_v2.shape[0] > 0 :
            merge_data.append(kegg_result_rankc_v2)

    kegg_result = pd.concat(merge_data, axis = 0)

    kegg_result.to_csv(workdir + '/kegg_result.csv')
    kegg_result_f = kegg_result[kegg_result['qvalue'] < GROUP_THRES].reset_index(drop=True)
    
    return kegg_result_f