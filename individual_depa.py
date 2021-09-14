#!/usr/bin/
# -*- coding: utf-8 -*-

__author__ = 'Liu Yachen'


import click
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
from matplotlib_venn import venn2,venn2_circles
from deps_lib import utils, raw_data, stable_pairs, methods_comp, penda_pro, methods_lib, similarity_lib, type1_error_lib, robustness_lib, kegg_lib, survival_lib

import rpy2
from rpy2.robjects.vectors import FloatVector
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
r_stats = importr('stats')

plt.rcParams['font.sans-serif'] = ['SimHei']  
sns.set(font_scale=1.5)  
plt.rcParams['axes.unicode_minus']=False 


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    pass

@click.version_option(version='1.0.0')
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--parameters_path', '-p', required=True, type=click.Path(), help='parameters file')
def stable(parameters_path):
    """ Get reversal pairs and consensus pairs """
    
    ########################## Read parameters ##############################
    params = configparser.ConfigParser()
    params.read(parameters_path)

    # workdir
    workdir = params.get('workdir', 'workdir')
    sp_workdir = params.get('workdir', 'sp_workdir')

    # script
    r_bpca_path = params.get('script path', 'r_bpca_path')
    r_penda_fdr_path = params.get('script path', 'r_penda_fdr_path')
    r_penda_path = params.get('script path', 'r_penda_path')
    reoa_path = params.get('script path', 'reoa_path')

    # preprocess parameters
    LOG_LABEL = eval(params.get('preprocess', 'log_label'))
    IMPUT_LABEL = eval(params.get('preprocess', 'imput_label'))
    NA_LABEL = eval(params.get('preprocess', 'na_label'))
    NA_RATE = eval(params.get('preprocess', 'na_rate'))
    NORMALIZATION = eval(params.get('preprocess', 'normalization'))

    INDEX_OLD_CHAR = eval(params.get('preprocess', 'index_old_char'))
    INDEX_NEW_CHAR = eval(params.get('preprocess', 'index_new_char'))
    COLUMNS_OLD_CHAR = eval(params.get('preprocess', 'columns_OLD_char'))
    COLUMNS_NEW_CHAR = eval(params.get('preprocess', 'columns_new_char'))

    # method parameters
    CYCLE_RANKC = eval(params.get('method parameters', 'cycle_rankc'))
    FDR = eval(params.get('method parameters', 'fdr'))
    MAX_EXCEPTION_RANKC = eval(params.get('method parameters', 'max_exception_rankc'))

    # stable pairs parameter
    data_path = params.get('stable pairs', 'data_path')
    normal_cohort_path = params.get('stable pairs', 'normal_cohort_path')
    tumor_cohort_path = params.get('stable pairs', 'tumor_cohort_path')

    HAS_SPECIFIC_PROTEIN = eval(params.get('stable pairs', 'label_specific_protein'))
    if HAS_SPECIFIC_PROTEIN:
        specific_protein_path = params.get('stable pairs', 'specific_protein_path')
    else:
        specific_protein_path = None

    SP_THRES = eval(params.get('stable pairs', 'sp_thres'))
    RANDOM_VISUAL = eval(params.get('stable pairs', 'random_visual'))
    N_VISUAL = eval(params.get('stable pairs', 'n_visual'))
    
    # read data
    data = utils.read_data(data_path)
    normal_cohort = utils.read_normal_cohort(normal_cohort_path)
    tumor_cohort = utils.read_tumor_cohort(tumor_cohort_path)

    if HAS_SPECIFIC_PROTEIN:
        specific_protein = utils.read_specific_protein(specific_protein_path)
    else:
        specific_protein = None

    rd = raw_data.InputData(data=data, 
                            normal_cohort=normal_cohort, 
                            tumor_cohort=tumor_cohort, 
                            specific_protein=specific_protein, 
                            HAS_SPECIFIC_PROTEIN=HAS_SPECIFIC_PROTEIN,
                            paired_data=None, 
                            paired_samples=None, 
                            HAS_PAIRED_DATA=None,  
                            INDEX_OLD_CHAR=INDEX_OLD_CHAR, 
                            INDEX_NEW_CHAR=INDEX_NEW_CHAR, 
                            COLUMNS_OLD_CHAR=COLUMNS_OLD_CHAR, 
                            COLUMNS_NEW_CHAR=COLUMNS_NEW_CHAR,
                            NORMALIZATION=NORMALIZATION, 
                            LOG_LABEL=LOG_LABEL, 
                            IMPUT_LABEL=IMPUT_LABEL,
                            NA_LABEL=NA_LABEL, 
                            NA_RATE=NA_RATE)

    preprocess_workdir = sp_workdir + '/preprocess'
    
    # data preprocess
    rd.data_preprocess(preprocess_workdir=preprocess_workdir,
                       r_bpca_path=r_bpca_path)

    sp = stable_pairs.StablePairs(rd = rd, 
                                 sp_workdir = sp_workdir, 
                                 SP_THRES = SP_THRES, 
                                 RANDOM_VISUAL = RANDOM_VISUAL, 
                                 NUM_VISUAL = N_VISUAL, 
                                 CYCLE_RANKC = CYCLE_RANKC, 
                                 FDR = FDR, 
                                 MAX_EXCEPTION_RANKC = MAX_EXCEPTION_RANKC, 
                                 fig_save_path=sp_workdir,
                                 reoa_path = reoa_path)

    sp.run_stablePairs()

    

@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--parameters_path', '-p', required=True, type=click.Path(), help='parameters file')
def comparison(parameters_path):
    ''' Compare the precision of DEAs with the number of DEPs '''

    ########################## parameters input ##############################
    params = configparser.ConfigParser()
    params.read(parameters_path)

    # workdir
    workdir = params.get('workdir', 'workdir')
    mc_workdir = params.get('workdir', 'mc_workdir')

    # script
    r_bpca_path = params.get('script path', 'r_bpca_path')
    r_penda_fdr_path = params.get('script path', 'r_penda_fdr_path')
    r_penda_path = params.get('script path', 'r_penda_path')
    reoa_path = params.get('script path', 'reoa_path')

    # preprocess parameters
    LOG_LABEL = eval(params.get('preprocess', 'log_label'))
    IMPUT_LABEL = eval(params.get('preprocess', 'imput_label'))
    NA_LABEL = eval(params.get('preprocess', 'na_label'))
    NA_RATE = eval(params.get('preprocess', 'na_rate'))
    NORMALIZATION = eval(params.get('preprocess', 'normalization'))

    INDEX_OLD_CHAR = eval(params.get('preprocess', 'index_old_char'))
    INDEX_NEW_CHAR = eval(params.get('preprocess', 'index_new_char'))
    COLUMNS_OLD_CHAR = eval(params.get('preprocess', 'columns_OLD_char'))
    COLUMNS_NEW_CHAR = eval(params.get('preprocess', 'columns_new_char'))

    # method parameters
    CYCLE_RANKC = eval(params.get('method parameters', 'cycle_rankc'))
    FDR = eval(params.get('method parameters', 'fdr'))
    MAX_EXCEPTION_RANKC = eval(params.get('method parameters', 'max_exception_rankc'))

    # methods comparison parameter
    data_path = params.get('methods comparison', 'data_path')
    normal_cohort_path = params.get('methods comparison', 'normal_cohort_path')
    tumor_cohort_path = params.get('methods comparison', 'tumor_cohort_path')
    paired_data_path = params.get('methods comparison', 'paired_data_path')
    paired_samples_path = params.get('methods comparison', 'paired_samples_path')
    METHODS_LIST = eval(params.get('methods comparison', 'methods_list'))

    # read data
    data = utils.read_data(data_path)
    normal_cohort = utils.read_normal_cohort(normal_cohort_path)
    tumor_cohort = utils.read_tumor_cohort(tumor_cohort_path)
    paired_data = utils.read_paired_data(paired_data_path)
    paired_samples = utils.read_paired_samples(paired_samples_path)

    rd = raw_data.InputData(data=data, 
                            normal_cohort=normal_cohort, 
                            tumor_cohort=tumor_cohort, 
                            specific_protein=None, 
                            HAS_SPECIFIC_PROTEIN=False,
                            paired_data=paired_data, 
                            paired_samples=paired_samples, 
                            HAS_PAIRED_DATA=True,  
                            INDEX_OLD_CHAR=INDEX_OLD_CHAR, 
                            INDEX_NEW_CHAR=INDEX_NEW_CHAR, 
                            COLUMNS_OLD_CHAR=COLUMNS_OLD_CHAR, 
                            COLUMNS_NEW_CHAR=COLUMNS_NEW_CHAR,
                            NORMALIZATION=NORMALIZATION, 
                            LOG_LABEL=LOG_LABEL, 
                            IMPUT_LABEL=IMPUT_LABEL,
                            NA_LABEL=NA_LABEL, 
                            NA_RATE=NA_RATE)

    preprocess_workdir = mc_workdir + '/preprocess'

    rd.data_preprocess(preprocess_workdir=preprocess_workdir,
                       r_bpca_path=r_bpca_path)

    mc = methods_comp.methodsComp(rd = rd, 
                                  r_penda_path = r_penda_path, 
                                  r_penda_fdr_path = r_penda_fdr_path,
                                  reoa_path = reoa_path,
                                  CYCLE_RANKC = CYCLE_RANKC, 
                                  FDR = FDR, 
                                  MAX_EXCEPTION_RANKC = MAX_EXCEPTION_RANKC, 
                                  PLOT_METHODS_COMPARE_RESULT = False,
                                  METHODS_LIST = METHODS_LIST)

    precision_mess, positive_mess = mc.run_methodsComp(mc_workdir = mc_workdir)
    
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--parameters_path', '-p', required=True, type=click.Path(), help='parameters file')
def evaluation(parameters_path):
    ''' Parameter evaluation '''
    
    ########################## parameters input ##############################
    params = configparser.ConfigParser()
    params.read(parameters_path)

    # workdir
    wordir = params.get('workdir', 'workdir')
    pe_workdir = params.get('workdir', 'pe_workdir')

    # script
    r_bpca_path = params.get('script path', 'r_bpca_path')
    r_penda_fdr_path = params.get('script path', 'r_penda_fdr_path')
    r_penda_path = params.get('script path', 'r_penda_path')
    reoa_path = params.get('script path', 'reoa_path')

    # preprocess parameters
    LOG_LABEL = eval(params.get('preprocess', 'log_label'))
    IMPUT_LABEL = eval(params.get('preprocess', 'imput_label'))
    NA_LABEL = eval(params.get('preprocess', 'na_label'))
    NA_RATE = eval(params.get('preprocess', 'na_rate'))
    NORMALIZATION = eval(params.get('preprocess', 'normalization'))

    INDEX_OLD_CHAR = eval(params.get('preprocess', 'index_old_char'))
    INDEX_NEW_CHAR = eval(params.get('preprocess', 'index_new_char'))
    COLUMNS_OLD_CHAR = eval(params.get('preprocess', 'columns_OLD_char'))
    COLUMNS_NEW_CHAR = eval(params.get('preprocess', 'columns_new_char'))

    # method parameters
    CYCLE_RANKC = eval(params.get('method parameters', 'cycle_rankc'))
    FDR = eval(params.get('method parameters', 'fdr'))
    MAX_EXCEPTION_RANKC = eval(params.get('method parameters', 'max_exception_rankc'))

    # parameters evalution
    data_path = params.get('parameters evalution', 'data_path')
    normal_cohort_path = params.get('parameters evalution', 'normal_cohort_path')
    tumor_cohort_path = params.get('parameters evalution', 'tumor_cohort_path')
    samples_path = params.get('parameters evalution', 'samples_path')

    SAMPLE_SIZE_LABEL = eval(params.get('parameters evalution', 'sample_size_label'))
    N_PROTEIN_LABEL = eval(params.get('parameters evalution', 'n_protein_label'))
    METHODS_LIST = eval(params.get('parameters evalution', 'methods_list'))
    SAMPLE_SIZE_LIST = eval(params.get('parameters evalution', 'sample_size_list'))
    N_PROTEIN_LIST = eval(params.get('parameters evalution', 'n_protein_list'))

    # read data
    data = utils.read_data(data_path)
    normal_cohort = utils.read_normal_cohort(normal_cohort_path)
    tumor_cohort = utils.read_tumor_cohort(tumor_cohort_path)
    samples = pd.read_csv(samples_path, sep='\t', index_col=0)
    samples = utils.col_adapt(samples_col = samples, characters = COLUMNS_OLD_CHAR, new_char = COLUMNS_NEW_CHAR)

    rd = raw_data.InputData(data=data, 
                            normal_cohort=normal_cohort, 
                            tumor_cohort=tumor_cohort, 
                            specific_protein=None, 
                            HAS_SPECIFIC_PROTEIN=False,
                            paired_data=None, 
                            paired_samples=None, 
                            HAS_PAIRED_DATA=False,  
                            INDEX_OLD_CHAR=INDEX_OLD_CHAR, 
                            INDEX_NEW_CHAR=INDEX_NEW_CHAR, 
                            COLUMNS_OLD_CHAR=COLUMNS_OLD_CHAR, 
                            COLUMNS_NEW_CHAR=COLUMNS_NEW_CHAR,
                            NORMALIZATION=NORMALIZATION, 
                            LOG_LABEL=LOG_LABEL, 
                            IMPUT_LABEL=IMPUT_LABEL,
                            NA_LABEL=NA_LABEL, 
                            NA_RATE=NA_RATE)

    preprocess_workdir = pe_workdir + '/preprocess'
    rd.data_preprocess(preprocess_workdir=preprocess_workdir,
                       r_bpca_path=r_bpca_path)

    samples.to_csv(preprocess_workdir+'/samples.txt', sep='\t')
    data = rd.data_imput
    normal_cohort = rd.normal_cohort_adapt
    tumor_cohort = rd.tumor_cohort_adapt

    if SAMPLE_SIZE_LABEL:
        workdir_ss = os.path.abspath(pe_workdir) + '/sample_size'

        if os.path.exists(workdir_ss) == False:
            os.makedirs(workdir_ss)

        ss_path_list = []

        ss_prec_mess = {}
        ss_posi_mess = {}

        i = 0
        for ss in SAMPLE_SIZE_LIST:
            i += 1
            t = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print('\nSAMPLE SIZE LIST: %s/%s    %s\n'%(i, len(SAMPLE_SIZE_LIST),t))

            _path = workdir_ss + '/%s'%ss
            ss_path_list.append(_path)

            if os.path.exists(_path) == False:
                os.makedirs(_path)
            os.chdir(_path)

            _samples_part = samples.loc[random.sample(samples.index.tolist(), ss * 2), :].reset_index(drop=True)

            _paired_index = random.sample(_samples_part.index.tolist(), ss)
            _paired_index.sort()

            _unpaired_index = []
            for _s in _samples_part.index:
                if _s not in _paired_index:
                    _unpaired_index.append(_s)

            _paired_samples = _samples_part.loc[_paired_index, :].reset_index(drop=True)
            _unpaired_samples = _samples_part.loc[_unpaired_index, :].reset_index(drop=True)

            _tumor_cohort = pd.DataFrame({'tumor': _paired_samples['tumor']})
            _normal_cohort = pd.DataFrame({'normal': _unpaired_samples['normal']})
            _paired_data = data.loc[:, _paired_samples['normal'].tolist() + _paired_samples['tumor'].tolist()]
            _data = data.loc[:, _paired_samples['tumor'].tolist() + _unpaired_samples['normal'].tolist()]

            _paired_samples.to_csv('./paired_samples.txt', sep='\t', index=False)
            _paired_data.to_csv('./paired_data.csv')
            _tumor_cohort.to_csv('./tumor.txt', sep='\t', index=False)
            _normal_cohort.to_csv('./normal.txt', sep='\t', index=False)
            _data.to_csv('./data.csv')

            rd = raw_data.InputData(data=_data, 
                                    normal_cohort=_normal_cohort, 
                                    tumor_cohort=_tumor_cohort, 
                                    specific_protein=None, 
                                    HAS_SPECIFIC_PROTEIN=False,
                                    paired_data=_paired_data, 
                                    paired_samples=_paired_samples, 
                                    HAS_PAIRED_DATA=True,  
                                    INDEX_OLD_CHAR=['-', ' '], 
                                    INDEX_NEW_CHAR='.', 
                                    COLUMNS_OLD_CHAR=['-', ' '], 
                                    COLUMNS_NEW_CHAR='.',
                                    NORMALIZATION=False, 
                                    LOG_LABEL=False, 
                                    NA_LABEL='False', 
                                    NA_RATE=0)

            _preprocess_workdir = _path + '/preprocess'

            rd.data_preprocess(preprocess_workdir=_preprocess_workdir,
                               r_bpca_path=r_bpca_path)

            mc = methods_comp.methodsComp(rd = rd, 
                                          r_penda_path = r_penda_path, 
                                          r_penda_fdr_path = r_penda_fdr_path,
                                          reoa_path = reoa_path,
                                          CYCLE_RANKC = CYCLE_RANKC, 
                                          FDR = FDR, 
                                          MAX_EXCEPTION_RANKC = MAX_EXCEPTION_RANKC, 
                                          PLOT_METHODS_COMPARE_RESULT = False,
                                          METHODS_LIST = METHODS_LIST)

            _mc_workdir = _path + '/methods_compare'
            precision_mess, positive_mess = mc.run_methodsComp(mc_workdir = _mc_workdir)

            ss_prec_mess[ss] = precision_mess
            ss_posi_mess[ss] = positive_mess

        prec_set = ss_prec_mess
        num_set = ss_posi_mess

        x_x = range(len(prec_set.keys()))
        x = list(prec_set.keys())
        y_p, y_n = {}, {}

        for _m in list(prec_set.values())[0].keys():
            _m_p = []
            _m_n = []

            for _x in x:
                _m_p.append(np.mean(prec_set[_x][_m]))
                _m_n.append(np.mean(num_set[_x][_m]))
            y_p[_m] = _m_p
            y_n[_m] = _m_n


        bar_width = 0.1
        tick_label = list(prec_set.keys())
        fig = plt.figure(figsize=(12,8))

        ax1 = fig.add_subplot(111)

        if 'RankComp' in METHODS_LIST:
            ax1.bar(np.array(x_x)+(0 * bar_width), y_p['rankc_v1'], bar_width, align="center", label='prec_rankc_v1', alpha=1, color='b')
            ax1.bar(np.array(x_x)+(1 * bar_width), y_p['rankc_v2'], bar_width, align="center", label='prec_rankc_v2', alpha=1, color='c')
        if 'Penda' in METHODS_LIST:
            ax1.bar(np.array(x_x)+(2 * bar_width), y_p['penda'], bar_width, align="center", label='prec_penda', alpha=1, color='g')
        if 'T-test' in METHODS_LIST:
            ax1.bar(np.array(x_x)+(3 * bar_width), y_p['ttest'], bar_width, align="center", label='Prec_ttest', alpha=1, color='m')
        if 'Wilcoxon' in METHODS_LIST:
            ax1.bar(np.array(x_x)+(4 * bar_width), y_p['wilcox'], bar_width, align="center", label='prec_wilcox', alpha=1, color='r')
        if 'Quantile' in METHODS_LIST:
            ax1.bar(np.array(x_x)+(5 * bar_width), y_p['quantile'], bar_width, align="center", label='prec_quantile', alpha=1, color='y')
        if 'Peng_methods' in METHODS_LIST:
            ax1.bar(np.array(x_x)+(6 * bar_width), y_p['peng_method'], bar_width, align="center", label='prec_peng', alpha=1, color='k')

        ax1.set_xlabel("Sample Size", fontsize=20)
        ax1.set_ylabel('Precision', fontsize=20)
        ax1.legend(loc=[1.2, 0.5])

        ax2 = ax1.twinx()
        if 'RankComp' in METHODS_LIST:
            ax2.plot(np.array(x_x)+(3 * bar_width), y_n['rankc_v1'],  label='num_rankc_v1', alpha=1, ms=10, lw=3, marker='o', c='b')
            ax2.plot(np.array(x_x)+(3 * bar_width), y_n['rankc_v2'], label='num_rankc_v2',  alpha=1, ms=10, lw=3, marker='o', c='c')
        if 'Penda' in METHODS_LIST:
            ax2.plot(np.array(x_x)+(3 * bar_width), y_n['penda'], label='num_penda', alpha=1, ms=10, lw=3, marker='o', c='g')
        if 'T-test' in METHODS_LIST:
            ax2.plot(np.array(x_x)+(3 * bar_width), y_n['ttest'], label='num_ttest', alpha=1, ms=10, lw=3, marker='o', c='m')
        if 'Wilcoxon' in METHODS_LIST:
            ax2.plot(np.array(x_x)+(3 * bar_width), y_n['wilcox'], label='num_wilcox', alpha=1, ms=10, lw=3, marker='o', c='r')
        if 'Quantile' in METHODS_LIST:
            ax2.plot(np.array(x_x)+(3 * bar_width), y_n['quantile'], label='num_quantile', alpha=1, ms=10, lw=3, marker='o', c='y')
        if 'Peng_methods' in METHODS_LIST:
            ax2.plot(np.array(x_x)+(3 * bar_width), y_n['peng_method'], label='num_peng', alpha=1, ms=10, lw=3, marker='o', c='k')

        ax2.set_ylabel('Number of DE protein', fontsize='20')
        sns.despine(left=True, bottom=True)
        ax2.tick_params(labelsize=15)
        ax2.legend(loc=[1.2, 0])

        plt.xticks(np.array(x_x)+ 3*bar_width/2, tick_label)

        plt.savefig(pe_workdir + '/pe_sample_size.pdf', dpi=800, bbox_inches='tight')
    
    if N_PROTEIN_LABEL:
        workdir_np = os.path.abspath(pe_workdir) + '/n_protein'

        if os.path.exists(workdir_np) == False:
            os.makedirs(workdir_np)

        np_path_list = []
        np_prec_mess = {}
        np_posi_mess = {}

        i = 0

        for pn in N_PROTEIN_LIST:
            i += 1
            t = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print('\nNUM PROTEIN LIST: %s/%s    %s\n'%(i, len(N_PROTEIN_LIST),t))

            _path = workdir_np + '/%s'%pn
            np_path_list.append(_path)

            if os.path.exists(_path) == False:
                os.makedirs(_path)
            os.chdir(_path)

            _data_part = data.loc[random.sample(data.index.tolist(), pn), :]

            _paired_index = random.sample(samples.index.tolist(), int(samples.shape[0] / 2))
            _paired_index.sort()

            _unpaired_index = []
            for _s in samples.index:
                if _s not in _paired_index:
                    _unpaired_index.append(_s)

            _paired_samples = samples.loc[_paired_index, :].reset_index(drop=True)
            _unpaired_samples = samples.loc[_unpaired_index, :].reset_index(drop=True)

            _tumor_cohort = pd.DataFrame({'tumor': _paired_samples['tumor']})
            _normal_cohort = pd.DataFrame({'normal': _unpaired_samples['normal']})
            _paired_data = _data_part.loc[:, _paired_samples['normal'].tolist() + _paired_samples['tumor'].tolist()]
            _data = _data_part.loc[:, _paired_samples['tumor'].tolist() + _unpaired_samples['normal'].tolist()]

            _paired_samples.to_csv('./paired_samples.txt', sep='\t', index=False)
            _paired_data.to_csv('./paired_data.csv')
            _tumor_cohort.to_csv('./tumor.txt', sep='\t', index=False)
            _normal_cohort.to_csv('./normal.txt', sep='\t', index=False)
            _data.to_csv('./data.csv')

            rd = raw_data.InputData(data=_data, 
                                    normal_cohort=_normal_cohort, 
                                    tumor_cohort=_tumor_cohort, 
                                    specific_protein=None, 
                                    HAS_SPECIFIC_PROTEIN=False,
                                    paired_data=_paired_data, 
                                    paired_samples=_paired_samples, 
                                    HAS_PAIRED_DATA=True,  
                                    INDEX_OLD_CHAR=['-', ' '], 
                                    INDEX_NEW_CHAR='.', 
                                    COLUMNS_OLD_CHAR=['-', ' '], 
                                    COLUMNS_NEW_CHAR='.',
                                    NORMALIZATION=False, 
                                    LOG_LABEL=False, 
                                    NA_LABEL='False', 
                                    NA_RATE=0)

            _preprocess_workdir = _path + '/preprocess'

            rd.data_preprocess(preprocess_workdir=_preprocess_workdir,
                               r_bpca_path=r_bpca_path)

            mc = methods_comp.methodsComp(rd = rd, 
                                          r_penda_path = r_penda_path, 
                                          r_penda_fdr_path = r_penda_fdr_path,
                                          reoa_path = reoa_path,
                                          CYCLE_RANKC = CYCLE_RANKC, 
                                          FDR = FDR, 
                                          MAX_EXCEPTION_RANKC = MAX_EXCEPTION_RANKC, 
                                          PLOT_METHODS_COMPARE_RESULT = False,
                                          METHODS_LIST = METHODS_LIST)

            _mc_workdir = _path + '/methods_compare'
            precision_mess, positive_mess = mc.run_methodsComp(mc_workdir = _mc_workdir)

            np_prec_mess[pn] = precision_mess
            np_posi_mess[pn] = positive_mess

        prec_set = np_prec_mess
        num_set = np_posi_mess

        x_x = range(len(prec_set.keys()))
        x = list(prec_set.keys())
        y_p, y_n = {}, {}

        for _m in list(prec_set.values())[0].keys():
            _m_p = []
            _m_n = []

            for _x in x:
                _m_p.append(np.mean(prec_set[_x][_m]))
                _m_n.append(np.mean(num_set[_x][_m]))
            y_p[_m] = _m_p
            y_n[_m] = _m_n


        bar_width = 0.1
        tick_label = list(prec_set.keys())
        fig = plt.figure(figsize=(12,8))


        ax1 = fig.add_subplot(111)

        if 'RankComp' in METHODS_LIST:
            ax1.bar(np.array(x_x)+(0 * bar_width), y_p['rankc_v1'], bar_width, align="center", label='prec_rankc_v1', alpha=1, color='b')
            ax1.bar(np.array(x_x)+(1 * bar_width), y_p['rankc_v2'], bar_width, align="center", label='prec_rankc_v2', alpha=1, color='c')
        if 'Penda' in METHODS_LIST:
            ax1.bar(np.array(x_x)+(2 * bar_width), y_p['penda'], bar_width, align="center", label='prec_penda', alpha=1, color='g')
        if 'T-test' in METHODS_LIST:
            ax1.bar(np.array(x_x)+(3 * bar_width), y_p['ttest'], bar_width, align="center", label='Prec_ttest', alpha=1, color='m')
        if 'Wilcoxon' in METHODS_LIST:
            ax1.bar(np.array(x_x)+(4 * bar_width), y_p['wilcox'], bar_width, align="center", label='prec_wilcox', alpha=1, color='r')
        if 'Quantile' in METHODS_LIST:
            ax1.bar(np.array(x_x)+(5 * bar_width), y_p['quantile'], bar_width, align="center", label='prec_quantile', alpha=1, color='y')
        if 'Peng_methods' in METHODS_LIST:
            ax1.bar(np.array(x_x)+(6 * bar_width), y_p['peng_method'], bar_width, align="center", label='prec_peng', alpha=1, color='k')

        ax1.set_xlabel("The number of protein", fontsize=20)
        ax1.set_ylabel('Precision', fontsize=20)
        ax1.legend(loc=[1.2, 0.5])


        ax2 = ax1.twinx()   
        if 'RankComp' in METHODS_LIST:
            ax2.plot(np.array(x_x)+(3 * bar_width), y_n['rankc_v1'],  label='num_rankc_v1', alpha=1, ms=10, lw=3, marker='o', c='b')
            ax2.plot(np.array(x_x)+(3 * bar_width), y_n['rankc_v2'], label='num_rankc_v2',  alpha=1, ms=10, lw=3, marker='o', c='c')
        if 'Penda' in METHODS_LIST:
            ax2.plot(np.array(x_x)+(3 * bar_width), y_n['penda'], label='num_penda', alpha=1, ms=10, lw=3, marker='o', c='g')
        if 'T-test' in METHODS_LIST:
            ax2.plot(np.array(x_x)+(3 * bar_width), y_n['ttest'], label='num_ttest', alpha=1, ms=10, lw=3, marker='o', c='m')
        if 'Wilcoxon' in METHODS_LIST:
            ax2.plot(np.array(x_x)+(3 * bar_width), y_n['wilcox'], label='num_wilcox', alpha=1, ms=10, lw=3, marker='o', c='r')
        if 'Quantile' in METHODS_LIST:
            ax2.plot(np.array(x_x)+(3 * bar_width), y_n['quantile'], label='num_quantile', alpha=1, ms=10, lw=3, marker='o', c='y')
        if 'Peng_methods' in METHODS_LIST:
            ax2.plot(np.array(x_x)+(3 * bar_width), y_n['peng_method'], label='num_peng', alpha=1, ms=10, lw=3, marker='o', c='k')

        ax2.set_ylabel('Number of DE protein', fontsize='20')
        sns.despine(left=True, bottom=True)   
        ax2.tick_params(labelsize=15)
        ax2.legend(loc=[1.2, 0])

        plt.xticks(np.array(x_x)+ 3*bar_width/2, tick_label)

        plt.savefig(pe_workdir + '/pe_num_protein.pdf', dpi=800, bbox_inches='tight')

        
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--parameters_path', '-p', required=True, type=click.Path(), help='parameters file')
def similarity(parameters_path):
    ''' Compare the similarity of DEAs '''
    
    ########################## parameters input ##############################
    params = configparser.ConfigParser()
    params.read(parameters_path)

    # workdir
    similarity_workdir = params.get('workdir', 'similarity_workdir')

    # script
    r_bpca_path = params.get('script path', 'r_bpca_path')
    r_penda_fdr_path = params.get('script path', 'r_penda_fdr_path')
    r_penda_path = params.get('script path', 'r_penda_path')
    reoa_path = params.get('script path', 'reoa_path')

    # preprocess parameters
    LOG_LABEL = eval(params.get('preprocess', 'log_label'))
    IMPUT_LABEL = eval(params.get('preprocess', 'imput_label'))
    NA_LABEL = eval(params.get('preprocess', 'na_label'))
    NA_RATE = eval(params.get('preprocess', 'na_rate'))
    NORMALIZATION = eval(params.get('preprocess', 'normalization'))

    INDEX_OLD_CHAR = eval(params.get('preprocess', 'index_old_char'))
    INDEX_NEW_CHAR = eval(params.get('preprocess', 'index_new_char'))
    COLUMNS_OLD_CHAR = eval(params.get('preprocess', 'columns_OLD_char'))
    COLUMNS_NEW_CHAR = eval(params.get('preprocess', 'columns_new_char'))

    # method parameters
    CYCLE_RANKC = eval(params.get('method parameters', 'cycle_rankc'))
    FDR = eval(params.get('method parameters', 'fdr'))
    MAX_EXCEPTION_RANKC = eval(params.get('method parameters', 'max_exception_rankc'))

    # methods comparison parameter
    data_path = params.get('similarity', 'data_path')
    normal_cohort_path = params.get('similarity', 'normal_cohort_path')
    tumor_cohort_path = params.get('similarity', 'tumor_cohort_path')
    METHODS_LIST = eval(params.get('similarity', 'methods_list'))
    N_CC = eval(params.get('similarity', 'n_cc'))

    # read data
    data = utils.read_data(data_path)
    normal_cohort = utils.read_normal_cohort(normal_cohort_path)
    tumor_cohort = utils.read_tumor_cohort(tumor_cohort_path)

    rd = raw_data.InputData(data=data, 
                            normal_cohort=normal_cohort, 
                            tumor_cohort=tumor_cohort, 
                            specific_protein=None, 
                            HAS_SPECIFIC_PROTEIN=False,
                            paired_data=None, 
                            paired_samples=None, 
                            HAS_PAIRED_DATA=False,  
                            INDEX_OLD_CHAR=INDEX_OLD_CHAR, 
                            INDEX_NEW_CHAR=INDEX_NEW_CHAR, 
                            COLUMNS_OLD_CHAR=COLUMNS_OLD_CHAR, 
                            COLUMNS_NEW_CHAR=COLUMNS_NEW_CHAR,
                            NORMALIZATION=NORMALIZATION, 
                            LOG_LABEL=LOG_LABEL, 
                            IMPUT_LABEL=IMPUT_LABEL,
                            NA_LABEL=NA_LABEL, 
                            NA_RATE=NA_RATE)

    preprocess_workdir = similarity_workdir + '/preprocess'

    rd.data_preprocess(preprocess_workdir=preprocess_workdir,
                       r_bpca_path=r_bpca_path)

    mc = methods_comp.methodsComp(rd = rd, 
                                  r_penda_path = r_penda_path, 
                                  r_penda_fdr_path = r_penda_fdr_path,
                                  reoa_path = reoa_path,
                                  CYCLE_RANKC = CYCLE_RANKC, 
                                  FDR = FDR, 
                                  MAX_EXCEPTION_RANKC = MAX_EXCEPTION_RANKC, 
                                  PLOT_METHODS_COMPARE_RESULT = False,
                                  METHODS_LIST = METHODS_LIST)

    mc.run_methodsComp(mc_workdir = similarity_workdir+'/methods', method_comp_label = False)

    result_dir = similarity_workdir+'/methods'

    similarity_lib.get_algorithm_similarity(result_dir = result_dir, 
                                            workdir = similarity_workdir, 
                                            N = N_CC)
    
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--parameters_path', '-p', required=True, type=click.Path(), help='parameters file')
def type1error(parameters_path):
    ''' Type one error '''
    ########################## parameters input ##############################
    params = configparser.ConfigParser()
    params.read(parameters_path)

    # workdir
    # workdir = params.get('workdir', 'workdir')
    type1_error_workdir = params.get('workdir', 'type1_error_workdir')

    # script
    r_bpca_path = params.get('script path', 'r_bpca_path')
    r_penda_fdr_path = params.get('script path', 'r_penda_fdr_path')
    r_penda_path = params.get('script path', 'r_penda_path')
    reoa_path = params.get('script path', 'reoa_path')

    # preprocess parameters
    LOG_LABEL = eval(params.get('preprocess', 'log_label'))
    IMPUT_LABEL = eval(params.get('preprocess', 'imput_label'))
    NA_LABEL = eval(params.get('preprocess', 'na_label'))
    NA_RATE = eval(params.get('preprocess', 'na_rate'))
    NORMALIZATION = eval(params.get('preprocess', 'normalization'))

    INDEX_OLD_CHAR = eval(params.get('preprocess', 'index_old_char'))
    INDEX_NEW_CHAR = eval(params.get('preprocess', 'index_new_char'))
    COLUMNS_OLD_CHAR = eval(params.get('preprocess', 'columns_OLD_char'))
    COLUMNS_NEW_CHAR = eval(params.get('preprocess', 'columns_new_char'))

    # method parameters
    CYCLE_RANKC = eval(params.get('method parameters', 'cycle_rankc'))
    FDR = eval(params.get('method parameters', 'fdr'))
    MAX_EXCEPTION_RANKC = eval(params.get('method parameters', 'max_exception_rankc'))

    # methods comparison parameter
    data_path_list = eval(params.get('type 1 error', 'data_path_list'))
    normal_col_path_list = eval(params.get('type 1 error', 'normal_col_path_list'))
    tumor_col_path_list = eval(params.get('type 1 error', 'tumor_col_path_list'))
    log_label_list = eval(params.get('type 1 error', 'log_label_list'))
    na_label_list = eval(params.get('type 1 error', 'na_label_list'))
    data_label_list = eval(params.get('type 1 error', 'data_label_list'))
    nd_sample_size = eval(params.get('type 1 error', 'nd_sample_size'))
    METHODS_LIST = eval(params.get('type 1 error', 'methods_list'))

    workdir = os.path.abspath(type1_error_workdir)
    imputation_workspace = type1_error_workdir + '/imputation'

    data_imput_set = {}
    data_normal_imput_set = {}
    data_tumor_imput_set = {}
    normal_col_adapt_set = {}
    tumor_col_adapt_set = {}

    print('Preprocess')
    for i in range(len(data_label_list)):
        t = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print('# %s: %s'%(data_label_list[i], t))
        data_imput, data_normal_imput, data_tumor_imput, normal_col_adapt, tumor_col_adapt = type1_error_lib.data_preprocess(data_path = data_path_list[i],
                                                                                                                             normal_col_path = normal_col_path_list[i],
                                                                                                                             tumor_col_path = tumor_col_path_list[i],
                                                                                                                             r_bpca_path = r_bpca_path,
                                                                                                                             workdir = imputation_workspace,
                                                                                                                             data_label = data_label_list[i],
                                                                                                                             LOG_LABEL = log_label_list[i], 
                                                                                                                             NA_LABEL = na_label_list[i], 
                                                                                                                             NA_RATE = NA_RATE,
                                                                                                                             INDEX_OLD_CHAR = INDEX_OLD_CHAR,
                                                                                                                             INDEX_NEW_CHAR = INDEX_NEW_CHAR,
                                                                                                                             COLUMNS_OLD_CHAR = COLUMNS_OLD_CHAR,
                                                                                                                             COLUMNS_NEW_CHAR = COLUMNS_NEW_CHAR)
    


        data_imput_set[data_label_list[i]] = data_imput
        data_normal_imput_set[data_label_list[i]] = data_normal_imput
        data_tumor_imput_set[data_label_list[i]] = data_tumor_imput
        normal_col_adapt_set[data_label_list[i]] = normal_col_adapt
        tumor_col_adapt_set[data_label_list[i]] = tumor_col_adapt

    nd_data_list = []
    nd_cohort_1_list = []
    nd_cohort_2_list = []
    for k in data_normal_imput_set.keys():
        _nd_data_list, _nd_cohort_1_list, _nd_cohort_2_list = type1_error_lib.create_nulldata(data_normal_imput = data_normal_imput_set[k], 
                                                                                             normal_col_adapt = normal_col_adapt_set[k], 
                                                                                             ND_SAMPLE_SIZE = nd_sample_size)
        nd_data_list = nd_data_list + _nd_data_list
        nd_cohort_1_list = nd_cohort_1_list + _nd_cohort_1_list
        nd_cohort_2_list = nd_cohort_2_list + _nd_cohort_2_list

    null_data_workspace = workdir + '/null_data'
    if os.path.exists(null_data_workspace) == False:
        os.makedirs(null_data_workspace)
    os.chdir(null_data_workspace)

    # null data
    tmp_workspace_list = []
    for i in range(len(nd_data_list)):
        tmp_workspace = null_data_workspace + '/' + '%s'%i
        if os.path.exists(tmp_workspace) == False:
            os.makedirs(tmp_workspace)
        os.chdir(tmp_workspace)

        nd_cohort_1_list[i].to_csv('./cohort_1.dat', sep='\t', index=0, header=0)
        nd_cohort_2_list[i].to_csv('./cohort_2.dat', sep='\t', index=0, header=0)
        nd_cohort_1_list[i].to_csv('./cohort_1.csv')
        nd_cohort_2_list[i].to_csv('./cohort_2.csv')
        tmp_workspace_list.append(tmp_workspace)

    nd_rankc = []
    nd_rankc_v1 = []
    nd_rankc_v2 = []
    nd_penda_fdr_result = []
    nd_penda_result = []
    nd_ttest_result = []
    nd_wilcox_result = []
    nd_penda_pro_result = []
    nd_quantile_result = []
    nd_peng_result = []

    print('Analysis')
    for i in range(len(nd_data_list)):

        t = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print('\nData index: %s (%s)'%(i, t))
        _path = tmp_workspace_list[i]

        os.chdir(_path)

        if 'RankComp' in METHODS_LIST:
            # rankc
            t = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print('# Rankc %s'%(t))
            type1_error_lib.rankc_j2(reoa_path = reoa_path, normal_data = nd_cohort_1_list[i], tumor_data = nd_cohort_2_list[i], 
                                     workdir=_path, cycle_rankc=CYCLE_RANKC, fdr=FDR, max_exception_rankc=MAX_EXCEPTION_RANKC)
            _dep = type1_error_lib.get_rankc_j2_results(_path)
            _dep = _dep.iloc[:, :-1]
            _dep.index = nd_cohort_1_list[i].index
            _dep.columns = nd_cohort_1_list[i].columns
            nd_rankc.append(_dep)
            nd_rankc_v1.append(_dep // 4)
            nd_rankc_v2.append(_dep % 4)

        if 'Penda' in METHODS_LIST:
            # Penda
            t = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print('# Penda %s'%(t))
        #     _penda_fdr_up, _penda_fdr_down, _penda_fdr_result = run_penda_fdr(rscript_path = r_penda_fdr_path, workdir = _path, fdr = fdr)
            _penda_up, _penda_down, _penda_result = type1_error_lib.run_penda(rscript_path = r_penda_path, workdir = _path)
        #     nd_penda_fdr_result.append(_penda_fdr_result)
            nd_penda_result.append(_penda_result)

        if 'Peng method' in METHODS_LIST:
            #Peng
            t = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print('# Peng %s'%(t))
            _peng_up, _peng_down = methods_lib.run_peng_method(normal_data = nd_cohort_1_list[i], 
                                                               tumor_data = nd_cohort_2_list[i])
            _peng_up.to_csv('./peng_up.csv')
            _peng_down.to_csv('./peng_down.csv')
            _peng_result = _peng_up | _peng_down
            nd_peng_result.append(_peng_result)

        if 'T-test' in METHODS_LIST:
            # T-test
            t = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print('# T-test %s'%(t))
            _ttest_up, _ttest_down = type1_error_lib.run_ttest(normal_run_path='./cohort_1.csv', tumor_test_path='./cohort_2.csv', workdir=_path, 
                                                               Q_THRES = FDR, SAVE_OUT = True)
            _ttest_result = _ttest_up + _ttest_down
            nd_ttest_result.append(_ttest_result)

        if 'Wilcoxon' in METHODS_LIST:
            # wilcoxon
            t = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print('# wilcox %s'%(t))
            _wilcox_up, _wilcox_down = type1_error_lib.run_wilcox(normal_run_path='./cohort_1.csv', tumor_test_path='./cohort_2.csv', workdir=_path, 
                                                                  Q_THRES = FDR, SAVE_OUT = True)
            _wilcox_result = _wilcox_up + _wilcox_down
            nd_wilcox_result.append(_wilcox_result)

        if 'Penda pro' in METHODS_LIST:
            # Penda pro
            t = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print('# Penda pro %s'%(t))
            _penda_pro_up_qv, _penda_pro_down_qv = penda_pro.run_penda_pro(normal_path = './cohort_1.csv', tumor_path = './cohort_2.csv', 
                                                                       FDR = fdr, CONVERGENCE_THRESHOLD = 0.99, 
                                                                       MAX_CYCLE = cycle_rankc, K = 30, THRESHOLD_LH = 0.99)
            _penda_pro_up = _penda_pro_up_qv < fdr
            _penda_pro_down = _penda_pro_down_qv < fdr
            nd_penda_pro_result = _penda_pro_up | _penda_pro_down

        if 'Quantile' in METHODS_LIST:
            # quantile
            t = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print('# Quantile %s'%(t))
            _quantile_up, _quantile_down = methods_lib.quantile_dep(normal_data_path = './cohort_1.csv', tumor_data_path = './cohort_2.csv', 
                                                        workdir = _path, quantile = 0.05, factor = 1.2, SAVE_OUT = True)
            _quantile_result = _quantile_up + _quantile_down
            nd_quantile_result.append(_quantile_result)

    nd_samples_fdp = []
    nd_method = []

    if 'RankComp' in METHODS_LIST:
        for i in range(len(nd_rankc_v1)):
            nd_samples_fdp.extend(((nd_rankc_v1[i] != 0).sum(axis=0) / nd_rankc_v1[i].shape[0]).values)
            for col in nd_rankc_v1[i].columns:
                nd_method.append('rankc_v1')

        for i in range(len(nd_rankc_v2)):
            nd_samples_fdp.extend(((nd_rankc_v2[i] != 0).sum(axis=0) / nd_rankc_v2[i].shape[0]).values)
            for col in nd_rankc_v2[i].columns:
                nd_method.append('rankc_v2')

    if 'Penda' in METHODS_LIST:
        for i in range(len(nd_penda_result)):
            nd_samples_fdp.extend(((nd_penda_result[i] != 0).sum(axis=0) / nd_penda_result[i].shape[0]).values)
            for col in nd_penda_result[i].columns:
                nd_method.append('penda')

    if 'Peng method' in METHODS_LIST:
        for i in range(len(nd_peng_result)):
            nd_samples_fdp.extend(((nd_peng_result[i] != 0).sum(axis=0) / nd_peng_result[i].shape[0]).values)
            for col in nd_peng_result[i].columns:
                nd_method.append('peng_method')

    if 'Penda pro' in METHODS_LIST:
        for i in range(len(nd_penda_pro_result)):
            nd_samples_fdp.extend(((nd_penda_pro_result[i] != 0).sum(axis=0) / nd_penda_pro_result[i].shape[0]).values)
            for col in nd_penda_pro_result[i].columns:
                nd_method.append('penda_pro')

    if 'T-test' in METHODS_LIST:
        for i in range(len(nd_ttest_result)):
            nd_samples_fdp.extend(((nd_ttest_result[i] != 0).sum(axis=0) / nd_ttest_result[i].shape[0]).values)
            for col in nd_ttest_result[i].columns:
                nd_method.append('T-test')

    if 'Wilcoxon' in METHODS_LIST:
        for i in range(len(nd_wilcox_result)):
            nd_samples_fdp.extend(((nd_wilcox_result[i] != 0).sum(axis=0) / nd_wilcox_result[i].shape[0]).values)
            for col in nd_wilcox_result[i].columns:
                nd_method.append('wilcoxon')

    if 'Quantile' in METHODS_LIST:
        for i in range(len(nd_quantile_result)):
            nd_samples_fdp.extend(((nd_quantile_result[i] != 0).sum(axis=0) / nd_quantile_result[i].shape[0]).values)
            for col in nd_quantile_result[i].columns:
                nd_method.append('quantile')

    nd_fpr_result = pd.DataFrame({'FDP': nd_samples_fdp, 'method': nd_method})

    plt.rcParams['font.sans-serif'] = ['SimHei']  
    sns.set(font_scale=1.5)  
    # plt.rcParams['savefig.dpi'] = 600 

    # plt.rcParams['figure.dpi'] = 300 
    plt.rcParams['axes.unicode_minus']=False 

    plt.figure(figsize=(12, 8))

    sns.boxplot(x = 'method', y='FDP', data=nd_fpr_result)
    sns.swarmplot(x = 'method', y="FDP", data=nd_fpr_result, color ='k',size = 3,alpha = 0.8)
    plt.hlines(0.05, -9, 9,color="red")

    plt.xlabel('methods', fontsize=20)
    plt.ylabel('FPR', fontsize=20)
    plt.title("Null data with FPR", fontsize=25)
    # plt.title("The precisions of DE protein by RankComp(v1-v2)-Penda-Wilcox: %s"%title, fontsize=25)
    # plt.legend(fontsize=20)
    plt.xticks(fontsize = 15)
    plt.yticks(fontsize = 15)

    plt.savefig(workdir +'/type_one_error.pdf', dpi=800, format='pdf') 
    plt.show()

    merge_list = []
    if 'RankComp' in METHODS_LIST:
        nd_param_sn_rankc_v1 = type1_error_lib.get_nd_param_sn(nd_result = nd_rankc_v1, 
                                                               nd_set = nd_data_list, 
                                                               method = 'rankc_v1')

        nd_param_sn_rankc_v2 = type1_error_lib.get_nd_param_sn(nd_result = nd_rankc_v2, 
                                                               nd_set = nd_data_list, 
                                                               method = 'rankc_v2')
        merge_list.append(nd_param_sn_rankc_v1)
        merge_list.append(nd_param_sn_rankc_v2)

    if 'Penda' in METHODS_LIST:
        nd_param_sn_penda = type1_error_lib.get_nd_param_sn(nd_result = nd_penda_result, 
                                                           nd_set = nd_data_list, 
                                                           method = 'penda')
        merge_list.append(nd_param_sn_penda)

    if 'Peng method' in METHODS_LIST:
        nd_param_sn_peng = type1_error_lib.get_nd_param_sn(nd_result = nd_peng_result, 
                                                           nd_set = nd_data_list, 
                                                           method = 'peng_method')
        merge_list.append(nd_param_sn_peng)

    if 'Penda pro' in METHODS_LIST:
        nd_param_sn_penda_pro = type1_error_lib.get_nd_param_sn(nd_result = nd_penda_pro_result, 
                                                               nd_set = nd_data_list, 
                                                               method = 'penda_pro')
        merge_list.append(nd_param_sn_penda_pro)

    if 'T-test' in METHODS_LIST:
        nd_param_sn_ttest = type1_error_lib.get_nd_param_sn(nd_result = nd_ttest_result, 
                                                           nd_set = nd_data_list, 
                                                           method = 'T-test')
        merge_list.append(nd_param_sn_ttest)

    if 'wilcoxon' in METHODS_LIST:
        nd_param_sn_wilcoxon = type1_error_lib.get_nd_param_sn(nd_result = nd_wilcox_result, 
                                                               nd_set = nd_data_list, 
                                                               method = 'wilcoxon')
        merge_list.append(nd_param_sn_wilcoxon)

    if 'Quantile' in METHODS_LIST:
        nd_param_sn_quantile = type1_error_lib.get_nd_param_sn(nd_result = nd_quantile_result, 
                                                               nd_set = nd_data_list, 
                                                               method = 'quantile')
        merge_list.append(nd_param_sn_quantile)


    nd_param_sn_result = pd.concat(merge_list, axis=0)

    plt.rcParams['font.sans-serif'] = ['SimHei']  
    sns.set(font_scale=1.5)  
    # plt.rcParams['savefig.dpi'] = 600

    # plt.rcParams['figure.dpi'] = 300 
    plt.rcParams['axes.unicode_minus']=False 

    plt.figure(figsize=(12, 8))

    sns.boxplot(x = 'method', y='sn', data=nd_param_sn_result[nd_param_sn_result['param_type'] == 'mean'])
    sns.swarmplot(x = 'method', y="sn", data=nd_param_sn_result[nd_param_sn_result['param_type'] == 'mean'], color ='k',size = 3,alpha = 0.8)
    plt.hlines(0, -9, 9,color="red")

    plt.xlabel('methods', fontsize=20)
    plt.ylabel('sn', fontsize=20)
    plt.title("Characteristics of false-positive genes: mean", fontsize=25)
    # plt.title("The precisions of DE protein by RankComp(v1-v2)-Penda-Wilcox: %s"%title, fontsize=25)
    # plt.legend(fontsize=20)
    plt.xticks(fontsize = 15)
    plt.yticks(fontsize = 15)
    plt.ylim(-10,20)
    plt.savefig(workdir +'/mean_influence.pdf', dpi=800, format='pdf') 
    plt.show()

    plt.rcParams['font.sans-serif'] = ['SimHei']
    sns.set(font_scale=1.5)  
    plt.rcParams['axes.unicode_minus']=False 

    plt.figure(figsize=(12, 8))

    sns.boxplot(x = 'method', y='sn', data=nd_param_sn_result[nd_param_sn_result['param_type'] == 'var'])
    sns.swarmplot(x = 'method', y="sn", data=nd_param_sn_result[nd_param_sn_result['param_type'] == 'var'], color ='k',size = 3,alpha = 0.8)
    plt.hlines(0, -9, 9,color="red")

    plt.xlabel('methods', fontsize=20)
    plt.ylabel('sn', fontsize=20)
    plt.title("Characteristics of false-positive genes: var", fontsize=25)
    # plt.title("The precisions of DE protein by RankComp(v1-v2)-Penda-Wilcox: %s"%title, fontsize=25)
    # plt.legend(fontsize=20)
    plt.xticks(fontsize = 15)
    plt.yticks(fontsize = 15)
    plt.ylim(-10,10)

    plt.savefig(workdir +'/var_influence.pdf', dpi=800, format='pdf') 
    plt.show()

    plt.rcParams['font.sans-serif'] = ['SimHei']  
    sns.set(font_scale=1.5)  
    plt.rcParams['axes.unicode_minus']=False 

    plt.figure(figsize=(12, 8))

    sns.boxplot(x = 'method', y='sn', data=nd_param_sn_result[nd_param_sn_result['param_type'] == 'cv'])
    sns.swarmplot(x = 'method', y="sn", data=nd_param_sn_result[nd_param_sn_result['param_type'] == 'cv'], color ='k',size = 3,alpha = 0.8)
    plt.hlines(0, -9, 9,color="red")

    plt.xlabel('methods', fontsize=20)
    plt.ylabel('sn', fontsize=20)
    plt.title("Characteristics of false-positive genes: cv", fontsize=25)
    # plt.title("The precisions of DE protein by RankComp(v1-v2)-Penda-Wilcox: %s"%title, fontsize=25)
    # plt.legend(fontsize=20)
    plt.xticks(fontsize = 15)
    plt.yticks(fontsize = 15)
    plt.ylim(-0.1,0.2)

    plt.savefig(workdir +'/cv_influence.pdf', dpi=800, format='pdf') 
    plt.show()
    
    
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--parameters_path', '-p', required=True, type=click.Path(), help='parameters file')
def robustness_individual(parameters_path):
    ''' Evaluation of algorithm robustness at the individual level '''
    ######################### parameters input #############################
    params = configparser.ConfigParser()
    params.read(parameters_path)

    # workdir
    robustness_workdir = params.get('workdir', 'ar_workdir')

    # script
    r_bpca_path = params.get('script path', 'r_bpca_path')
    r_penda_fdr_path = params.get('script path', 'r_penda_fdr_path')
    r_penda_path = params.get('script path', 'r_penda_path')
    reoa_path = params.get('script path', 'reoa_path')

    # preprocess parameters
    LOG_LABEL = eval(params.get('preprocess', 'log_label'))
    IMPUT_LABEL = eval(params.get('preprocess', 'imput_label'))
    NA_LABEL = eval(params.get('preprocess', 'na_label'))
    NA_RATE = eval(params.get('preprocess', 'na_rate'))
    NORMALIZATION = eval(params.get('preprocess', 'normalization'))

    INDEX_OLD_CHAR = eval(params.get('preprocess', 'index_old_char'))
    INDEX_NEW_CHAR = eval(params.get('preprocess', 'index_new_char'))
    COLUMNS_OLD_CHAR = eval(params.get('preprocess', 'columns_OLD_char'))
    COLUMNS_NEW_CHAR = eval(params.get('preprocess', 'columns_new_char'))

    # method parameters
    CYCLE_RANKC = eval(params.get('method parameters', 'cycle_rankc'))
    FDR = eval(params.get('method parameters', 'fdr'))
    MAX_EXCEPTION_RANKC = eval(params.get('method parameters', 'max_exception_rankc'))

    # methods comparison parameter
    data_1_path = params.get('robustness', 'data_1_path')
    tumor_cohort_1_path = params.get('robustness', 'tumor_cohort_1_path')
    normal_cohort_1_path = params.get('robustness', 'normal_cohort_1_path')

    data_2_path = params.get('robustness', 'data_2_path')
    tumor_cohort_2_path = params.get('robustness', 'tumor_cohort_2_path')
    normal_cohort_2_path = params.get('robustness', 'normal_cohort_2_path')

    log_label_list = eval(params.get('robustness', 'log_label_list'))
    na_label_list = eval(params.get('robustness', 'na_label_list'))
    imput_label_list = eval(params.get('robustness', 'imput_label_list'))

    robustness_individual_list = eval(params.get('robustness', 'robustness_individual_list'))
    robustness_group_list = eval(params.get('robustness', 'robustness_group_list'))

    N_CC = eval(params.get('robustness', 'n_cc'))

    data_1 = utils.read_data(data_1_path)
    normal_cohort_1 = utils.read_normal_cohort(normal_cohort_1_path)
    tumor_cohort_1 = utils.read_tumor_cohort(tumor_cohort_1_path)
    data_2 = utils.read_data(data_2_path)
    normal_cohort_2 = utils.read_normal_cohort(normal_cohort_2_path)
    tumor_cohort_2 = utils.read_tumor_cohort(tumor_cohort_2_path)

    rd1 = raw_data.InputData(data=data_1, 
                            normal_cohort=normal_cohort_1, 
                            tumor_cohort=tumor_cohort_1, 
                            specific_protein=None, 
                            HAS_SPECIFIC_PROTEIN=False,
                            paired_data=None, 
                            paired_samples=None, 
                            HAS_PAIRED_DATA=False,  
                            INDEX_OLD_CHAR=INDEX_OLD_CHAR, 
                            INDEX_NEW_CHAR=INDEX_NEW_CHAR, 
                            COLUMNS_OLD_CHAR=COLUMNS_OLD_CHAR, 
                            COLUMNS_NEW_CHAR=COLUMNS_NEW_CHAR,
                            NORMALIZATION=NORMALIZATION, 
                            LOG_LABEL=log_label_list[0], 
                            IMPUT_LABEL=imput_label_list[0],
                            NA_LABEL=na_label_list[0], 
                            NA_RATE=NA_RATE)

    rd2 = raw_data.InputData(data=data_2, 
                            normal_cohort=normal_cohort_2, 
                            tumor_cohort=tumor_cohort_2, 
                            specific_protein=None, 
                            HAS_SPECIFIC_PROTEIN=False,
                            paired_data=None, 
                            paired_samples=None, 
                            HAS_PAIRED_DATA=False,  
                            INDEX_OLD_CHAR=INDEX_OLD_CHAR, 
                            INDEX_NEW_CHAR=INDEX_NEW_CHAR, 
                            COLUMNS_OLD_CHAR=COLUMNS_OLD_CHAR, 
                            COLUMNS_NEW_CHAR=COLUMNS_NEW_CHAR,
                            NORMALIZATION=NORMALIZATION, 
                            LOG_LABEL=log_label_list[1], 
                            IMPUT_LABEL=imput_label_list[1],
                            NA_LABEL=na_label_list[1], 
                            NA_RATE=NA_RATE)

    preprocess_workdir_1 = robustness_workdir + '/preprocess_1'

    rd1.data_preprocess(preprocess_workdir=preprocess_workdir_1,
                       r_bpca_path=r_bpca_path)

    preprocess_workdir_2 = robustness_workdir + '/preprocess_2'

    rd2.data_preprocess(preprocess_workdir=preprocess_workdir_2,
                       r_bpca_path=r_bpca_path)
    
    new_index = list(set(rd1.data_imput.index).intersection(set(rd2.data_imput.index)))
    
    data_1_f = rd1.data_imput.loc[new_index, :]
    data_2_f = rd2.data_imput.loc[new_index, :]

    normal_data_1 = data_1_f.loc[:, rd1.normal_cohort_adapt['normal']]
    tumor_data_1 = data_1_f.loc[:, rd1.tumor_cohort_adapt['tumor']]
    normal_data_2 = data_2_f.loc[:, rd2.normal_cohort_adapt['normal']]
    tumor_data_2 = data_2_f.loc[:, rd2.tumor_cohort_adapt['tumor']]

    tumor_data = tumor_data_1

    # workdir
    individual_workdir = robustness_workdir + '/individual'

    if os.path.exists(individual_workdir) == False:
        os.makedirs(individual_workdir)
    os.chdir(individual_workdir)

    if 'RankComp' in robustness_individual_list:
        t = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print('# Rankcom %s'%(t))
        rankc_workspace1 = individual_workdir + '/rankcomp/1'
        if os.path.exists(rankc_workspace1) == True:
            shutil.rmtree(rankc_workspace1)
        os.makedirs(rankc_workspace1)
        os.chdir(rankc_workspace1)


        utils.write_normal_dat(normal_data_1, rankc_workspace1)
        utils.write_tumor_dat(tumor_data, rankc_workspace1)

        # run rankcomp v1 and v2
        utils.rankc_j2_qvalues(reoa_path, normal_data_1, tumor_data, rankc_workspace1, 
                               CYCLE_RANKC, FDR, MAX_EXCEPTION_RANKC)
        # get rankcomp result
        rankc_v1_up, rankc_v1_down, rankc_v2_up, rankc_v2_down = utils.get_rankc_j2_qvalues_result(rankc_workspace1, tumor_data)

        rankc_workspace2 = individual_workdir + '/rankcomp/2'
        if os.path.exists(rankc_workspace2) == True:
            shutil.rmtree(rankc_workspace2)
        os.makedirs(rankc_workspace2)
        os.chdir(rankc_workspace2)

        utils.write_normal_dat(normal_data_2, rankc_workspace2)
        utils.write_tumor_dat(tumor_data, rankc_workspace2)

        # run rankcomp v1 and v2
        utils.rankc_j2_qvalues(reoa_path, normal_data_2, tumor_data, rankc_workspace2, 
                               CYCLE_RANKC, FDR, MAX_EXCEPTION_RANKC)
        # get rankcomp result
        rankc_v1_up, rankc_v1_down, rankc_v2_up, rankc_v2_down = utils.get_rankc_j2_qvalues_result(rankc_workspace2, tumor_data)

    if 'Penda' in robustness_individual_list:
        t = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print('# Penda %s'%(t))
        penda_workspace1 = individual_workdir + '/penda/1'
        if os.path.exists(penda_workspace1) == True:
            shutil.rmtree(penda_workspace1)
        os.makedirs(penda_workspace1)
        os.chdir(penda_workspace1)

        # run penda
        utils.write_penda_data(normal_data_1, tumor_data, penda_workspace1)
        utils.run_penda(r_penda_path)

        penda_workspace2 = individual_workdir + '/penda/2'
        if os.path.exists(penda_workspace2) == True:
            shutil.rmtree(penda_workspace2)
        os.makedirs(penda_workspace2)
        os.chdir(penda_workspace2)

        # run penda
        utils.write_penda_data(normal_data_2, tumor_data, penda_workspace2)
        utils.run_penda(r_penda_path)

    if 'Penda pro' in robustness_individual_list:
        t = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print('# Penda pro%s'%(t))
        penda_pro_workspace1 = individual_workdir + '/penda_pro/1'
        if os.path.exists(penda_pro_workspace1) == True:
            shutil.rmtree(penda_pro_workspace1)
        os.makedirs(penda_pro_workspace1)
        os.chdir(penda_pro_workspace1)

        # run penda pro
        utils.write_penda_data(normal_data_1, tumor_data, penda_pro_workspace1)
        penda_up_pro_qv, penda_down_pro_qv = penda_pro.run_penda_pro(normal_path = './normal_run.csv', 
                                                                     tumor_path = './tumor_run.csv', 
                                                                     FDR = 0.05, 
                                                                     CONVERGENCE_THRESHOLD = 0.95, 
                                                                     MAX_CYCLE = 48, 
                                                                     K = 20, 
                                                                     THRESHOLD_LH = 0.99)
        penda_up_pro_qv.to_csv('./penda_up_pro_qv.csv')
        penda_down_pro_qv.to_csv('./penda_down_pro_qv.csv')

        penda_pro_workspace2 = individual_workdir + '/penda_pro/2'
        if os.path.exists(penda_pro_workspace2) == True:
            shutil.rmtree(penda_pro_workspace2)
        os.makedirs(penda_pro_workspace2)
        os.chdir(penda_pro_workspace2)

        # run penda pro
        utils.write_penda_data(normal_data_2, tumor_data, penda_pro_workspace2)
        penda_up_pro_qv, penda_down_pro_qv = penda_pro.run_penda_pro(normal_path = './normal_run.csv', 
                                                                     tumor_path = './tumor_run.csv', 
                                                                     FDR = 0.05, 
                                                                     CONVERGENCE_THRESHOLD = 0.95, 
                                                                     MAX_CYCLE = 48, 
                                                                     K = 20, 
                                                                     THRESHOLD_LH = 0.99)
        penda_up_pro_qv.to_csv('./penda_up_pro_qv.csv')
        penda_down_pro_qv.to_csv('./penda_down_pro_qv.csv')

    if 'T-test' in robustness_individual_list:
        t = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print('# T-test %s'%(t))
        ttest_workspace1 = individual_workdir + '/ttest/1'
        if os.path.exists(ttest_workspace1) == True:
            shutil.rmtree(ttest_workspace1)
        os.makedirs(ttest_workspace1)
        os.chdir(ttest_workspace1)

        utils.write_penda_data(normal_data_1, tumor_data, ttest_workspace1)

        ttest_up, ttest_down = methods_lib.run_ttest(normal_run_path = './normal_run.csv',
                                                     tumor_test_path = './tumor_run.csv',
                                                     workdir = ttest_workspace1,
                                                     Q_THRES = FDR,
                                                     D_THRES = 0,
                                                     SAVE_OUT = True)

        ttest_workspace2 = individual_workdir + '/ttest/2'
        if os.path.exists(ttest_workspace2) == True:
            shutil.rmtree(ttest_workspace2)
        os.makedirs(ttest_workspace2)
        os.chdir(ttest_workspace2)

        utils.write_penda_data(normal_data_2, tumor_data, ttest_workspace2)

        ttest_up, ttest_down = methods_lib.run_ttest(normal_run_path = './normal_run.csv',
                                                     tumor_test_path = './tumor_run.csv',
                                                     workdir = ttest_workspace2,
                                                     Q_THRES = FDR,
                                                     D_THRES = 0,
                                                     SAVE_OUT = True)

    if 'Wilcoxon' in robustness_individual_list:
        t = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print('# Wilcoxon %s'%(t))
        wilcox_workspace1 = individual_workdir + '/wilcoxon/1'
        if os.path.exists(wilcox_workspace1) == True:
            shutil.rmtree(wilcox_workspace1)
        os.makedirs(wilcox_workspace1)
        os.chdir(wilcox_workspace1)

        utils.write_penda_data(normal_data_1, tumor_data, wilcox_workspace1)

        wilcox_up, wilcox_down = methods_lib.run_wilcox(normal_run_path = './normal_run.csv', 
                                                        tumor_test_path = './tumor_run.csv',
                                                        workdir = wilcox_workspace1, 
                                                        Q_THRES = FDR, D_THRES = 0,
                                                        SAVE_OUT = True)


        wilcox_workspace2 = individual_workdir + '/wilcoxon/2'
        if os.path.exists(wilcox_workspace2) == True:
            shutil.rmtree(wilcox_workspace2)
        os.makedirs(wilcox_workspace2)
        os.chdir(wilcox_workspace2)

        utils.write_penda_data(normal_data_2, tumor_data, wilcox_workspace2)

        wilcox_up, wilcox_down = methods_lib.run_wilcox(normal_run_path = './normal_run.csv', 
                                                        tumor_test_path = './tumor_run.csv',
                                                        workdir = wilcox_workspace2, 
                                                        Q_THRES = FDR, D_THRES = 0,
                                                        SAVE_OUT = True)

    if 'Quantile' in robustness_individual_list:
        t = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print('# quantile %s'%(t))
        quantile_workspace1 = individual_workdir + '/quantile/1'
        if os.path.exists(quantile_workspace1) == True:
            shutil.rmtree(quantile_workspace1)
        os.makedirs(quantile_workspace1)
        os.chdir(quantile_workspace1)

        utils.write_penda_data(normal_data_1, tumor_data, quantile_workspace1)

        quantile_up, quantile_down = methods_lib.quantile_dep(normal_data_path = './normal_run.csv',
                                                              tumor_data_path = './tumor_run.csv',
                                                              workdir = quantile_workspace1, 
                                                              quantile = 0.05, 
                                                              factor = 1.2,
                                                              SAVE_OUT = True)


        quantile_workspace2 = individual_workdir + '/quantile/2'
        if os.path.exists(quantile_workspace2) == True:
            shutil.rmtree(quantile_workspace2)
        os.makedirs(quantile_workspace2)
        os.chdir(quantile_workspace2)

        utils.write_penda_data(normal_data_2, tumor_data, quantile_workspace2)

        quantile_up, quantile_down = methods_lib.quantile_dep(normal_data_path = './normal_run.csv',
                                                              tumor_data_path = './tumor_run.csv',
                                                              workdir = quantile_workspace2, 
                                                              quantile = 0.05, 
                                                              factor = 1.2,
                                                              SAVE_OUT = True)

    if 'RankComp' in robustness_individual_list:

        rankc_v1_1_up_qv_path = individual_workdir + '/rankcomp/1/rankc_v1_up_qvalues.csv'
        rankc_v1_2_up_qv_path = individual_workdir + '/rankcomp/2/rankc_v1_up_qvalues.csv'
        rankc_v1_1_down_qv_path = individual_workdir + '/rankcomp/1/rankc_v1_down_qvalues.csv'
        rankc_v1_2_down_qv_path = individual_workdir + '/rankcomp/2/rankc_v1_down_qvalues.csv'

        rankc_v2_1_up_qv_path = individual_workdir + '/rankcomp/1/rankc_v2_up_qvalues.csv'
        rankc_v2_2_up_qv_path = individual_workdir + '/rankcomp/2/rankc_v2_up_qvalues.csv'
        rankc_v2_1_down_qv_path = individual_workdir + '/rankcomp/1/rankc_v2_down_qvalues.csv'
        rankc_v2_2_down_qv_path = individual_workdir + '/rankcomp/2/rankc_v2_down_qvalues.csv'

        rankc_v1_1_up_qv = pd.read_csv(rankc_v1_1_up_qv_path, index_col=0)
        rankc_v1_2_up_qv = pd.read_csv(rankc_v1_2_up_qv_path, index_col=0)
        rankc_v1_1_down_qv = pd.read_csv(rankc_v1_1_down_qv_path, index_col=0)
        rankc_v1_2_down_qv = pd.read_csv(rankc_v1_2_down_qv_path, index_col=0)

        rankc_v2_1_up_qv = pd.read_csv(rankc_v2_1_up_qv_path, index_col=0)
        rankc_v2_2_up_qv = pd.read_csv(rankc_v2_2_up_qv_path, index_col=0)
        rankc_v2_1_down_qv = pd.read_csv(rankc_v2_1_down_qv_path, index_col=0)
        rankc_v2_2_down_qv = pd.read_csv(rankc_v2_2_down_qv_path, index_col=0)

        up_rankc_v1_rb = similarity_lib.get_concordance_curve_ar(rankc_v1_1_up_qv, rankc_v1_2_up_qv, N_CC)
        up_rankc_v2_rb = similarity_lib.get_concordance_curve_ar(rankc_v2_1_up_qv, rankc_v2_2_up_qv, N_CC)
        down_rankc_v1_rb = similarity_lib.get_concordance_curve_ar(rankc_v1_1_down_qv, rankc_v1_2_down_qv, N_CC)
        down_rankc_v2_rb = similarity_lib.get_concordance_curve_ar(rankc_v2_1_down_qv, rankc_v2_2_down_qv, N_CC)


    if 'Penda pro' in robustness_individual_list:

        penda_pro_1_up_qv_path = individual_workdir + '/penda_pro/1/penda_up_pro_qv.csv'
        penda_pro_2_up_qv_path = individual_workdir + '/penda_pro/2/penda_up_pro_qv.csv'
        penda_pro_1_down_qv_path = individual_workdir + '/penda_pro/1/penda_down_pro_qv.csv'
        penda_pro_2_down_qv_path = individual_workdir + '/penda_pro/2/penda_down_pro_qv.csv'

        penda_pro_1_up_qv = pd.read_csv(penda_pro_1_up_qv_path, index_col=0)
        penda_pro_2_up_qv = pd.read_csv(penda_pro_2_up_qv_path, index_col=0)
        penda_pro_1_down_qv = pd.read_csv(penda_pro_1_down_qv_path, index_col=0)
        penda_pro_2_down_qv = pd.read_csv(penda_pro_2_down_qv_path, index_col=0)

        up_pp_rb = similarity_lib.get_concordance_curve_ar(penda_pro_1_up_qv, penda_pro_2_up_qv, N_CC)
        down_pp_rb = similarity_lib.get_concordance_curve_ar(penda_pro_1_down_qv, penda_pro_2_down_qv, N_CC)

    if 'T-test' in robustness_individual_list:

        ttest_up_1_path = individual_workdir + '/ttest/1/ttest_up.csv'
        ttest_down_1_path = individual_workdir + '/ttest/1/ttest_down.csv'
        ttest_qvalues_1_path = individual_workdir + '/ttest/1/ttest_qvalues.csv'
        ttest_up_2_path = individual_workdir + '/ttest/2/ttest_up.csv'
        ttest_down_2_path = individual_workdir + '/ttest/2/ttest_down.csv'
        ttest_qvalues_2_path = individual_workdir + '/ttest/2/ttest_qvalues.csv'

        ttest_up_1 = pd.read_csv(ttest_up_1_path, index_col=0)
        ttest_down_1 = pd.read_csv(ttest_down_1_path, index_col=0)
        ttest_qvalues_1 = pd.read_csv(ttest_qvalues_1_path, index_col=0)
        ttest_up_2 = pd.read_csv(ttest_up_2_path, index_col=0)
        ttest_down_2 = pd.read_csv(ttest_down_2_path, index_col=0)
        ttest_qvalues_2 = pd.read_csv(ttest_qvalues_2_path, index_col=0)

        # get ttest_up and ttest_down qvalues
        ttest_up_1_qv = copy.deepcopy(ttest_up_1)
        ttest_down_1_qv = copy.deepcopy(ttest_down_1)

        ttest_up_1_qv.iloc[:,:] = 1
        ttest_down_1_qv.iloc[:,:] = 1

        for idx in ttest_up_1.index:
            if ttest_up_1.loc[idx, :].sum() == ttest_up_1.shape[1]:
                ttest_up_1_qv.loc[idx, :] = ttest_qvalues_1.loc[idx, :]
            elif ttest_down_1.loc[idx, :].sum() == ttest_down_1.shape[1]:
                ttest_down_1_qv.loc[idx, :] = ttest_qvalues_1.loc[idx, :]

        ttest_up_2_qv = copy.deepcopy(ttest_up_2)
        ttest_down_2_qv = copy.deepcopy(ttest_down_2)

        ttest_up_2_qv.iloc[:,:] = 1
        ttest_down_2_qv.iloc[:,:] = 1

        for idx in ttest_up_2.index:
            if ttest_up_2.loc[idx, :].sum() == ttest_up_2.shape[1]:
                ttest_up_2_qv.loc[idx, :] = ttest_qvalues_2.loc[idx, :]
            elif ttest_down_2.loc[idx, :].sum() == ttest_down_2.shape[1]:
                ttest_down_2_qv.loc[idx, :] = ttest_qvalues_2.loc[idx, :]

        up_tt_rb = similarity_lib.get_concordance_curve_ar(ttest_up_1_qv, ttest_up_2_qv, N_CC)
        down_tt_rb = similarity_lib.get_concordance_curve_ar(ttest_down_1_qv, ttest_down_2_qv, N_CC)

    if 'Wilcoxon' in robustness_individual_list:

        wilcox_up_1_path = individual_workdir + '/wilcoxon/1/wilcox_up.csv'
        wilcox_down_1_path = individual_workdir + '/wilcoxon/1/wilcox_down.csv'
        wilcox_qvalues_1_path = individual_workdir + '/wilcoxon/1/wilcox_qvalues.csv'
        wilcox_up_2_path = individual_workdir + '/wilcoxon/2/wilcox_up.csv'
        wilcox_down_2_path = individual_workdir + '/wilcoxon/2/wilcox_down.csv'
        wilcox_qvalues_2_path = individual_workdir + '/wilcoxon/2/wilcox_qvalues.csv'

        wilcox_up_1 = pd.read_csv(wilcox_up_1_path, index_col=0)
        wilcox_down_1 = pd.read_csv(wilcox_down_1_path, index_col=0)
        wilcox_qvalues_1 = pd.read_csv(wilcox_qvalues_1_path, index_col=0)
        wilcox_up_2 = pd.read_csv(wilcox_up_2_path, index_col=0)
        wilcox_down_2 = pd.read_csv(wilcox_down_2_path, index_col=0)
        wilcox_qvalues_2 = pd.read_csv(wilcox_qvalues_2_path, index_col=0)

        wilcox_up_1_qv = copy.deepcopy(wilcox_up_1)
        wilcox_down_1_qv = copy.deepcopy(wilcox_down_1)

        wilcox_up_1_qv.iloc[:,:] = 1
        wilcox_down_1_qv.iloc[:,:] = 1

        for idx in wilcox_up_1.index:
            if wilcox_up_1.loc[idx, :].sum() == wilcox_up_1.shape[1]:
                wilcox_up_1_qv.loc[idx, :] = wilcox_qvalues_1.loc[idx, :]
            elif wilcox_down_1.loc[idx, :].sum() == wilcox_down_1.shape[1]:
                wilcox_down_1_qv.loc[idx, :] = wilcox_qvalues_1.loc[idx, :]

        wilcox_up_2_qv = copy.deepcopy(wilcox_up_2)
        wilcox_down_2_qv = copy.deepcopy(wilcox_down_2)

        wilcox_up_2_qv.iloc[:,:] = 1
        wilcox_down_2_qv.iloc[:,:] = 1

        for idx in wilcox_up_2.index:
            if wilcox_up_2.loc[idx, :].sum() == wilcox_up_2.shape[1]:
                wilcox_up_2_qv.loc[idx, :] = wilcox_qvalues_2.loc[idx, :]
            elif wilcox_down_2.loc[idx, :].sum() == wilcox_down_2.shape[1]:
                wilcox_down_2_qv.loc[idx, :] = wilcox_qvalues_2.loc[idx, :]

        up_wil_rb = similarity_lib.get_concordance_curve_ar(wilcox_up_1_qv, wilcox_up_2_qv, N_CC)
        down_wil_rb = similarity_lib.get_concordance_curve_ar(wilcox_down_1_qv, wilcox_down_2_qv, N_CC)

    import matplotlib.pyplot as plt

    plt.figure(figsize=(12,8))
    if 'RankComp' in robustness_individual_list:
        plt.plot(up_rankc_v1_rb , linewidth=4, label='up_rankc_v1_rb: %.4f'%(similarity_lib.get_concordance_score(up_rankc_v1_rb)))
        plt.plot(up_rankc_v2_rb , linewidth=4, label='up_rankc_v2_rb: %.4f'%(similarity_lib.get_concordance_score(up_rankc_v2_rb)))
    if 'Penda pro' in robustness_individual_list:
        plt.plot(up_pp_rb , linewidth=4, label='up_pp_rb: %.4f'%(similarity_lib.get_concordance_score(up_pp_rb)))
    if 'T-test' in robustness_individual_list:
        plt.plot(up_tt_rb , linewidth=4, label='up_tt_rb: %.4f'%(similarity_lib.get_concordance_score(up_tt_rb)))
    if 'Wilcoxon' in robustness_individual_list:
        plt.plot(up_wil_rb , linewidth=4, label='up_wil_rb: %.4f'%(similarity_lib.get_concordance_score(up_wil_rb)))

    plt.legend(fontsize=15)
    plt.xlabel("Top N protein", fontsize=15)
    plt.ylabel("Concordance", fontsize=15)
    plt.title("Robustness of the algorithm", fontsize=20)
    plt.savefig(individual_workdir +'/up_algo_roubustness.pdf', dpi=800, format='pdf', bbox_inches='tight') #指定分辨
    # plt.savefig(individual_workdir + '/up_algo_roubustness.png', dpi=800)
    plt.show()

    import matplotlib.pyplot as plt

    plt.figure(figsize=(12,8))
    if 'Penda pro' in robustness_individual_list:
        plt.plot(down_pp_rb , linewidth=4, label='down_pp_rb: %.4f'%(similarity_lib.get_concordance_score(down_pp_rb)))
    if 'RankComp' in robustness_individual_list:
        plt.plot(down_rankc_v1_rb , linewidth=4, label='down_rankc_v1_rb: %.4f'%(similarity_lib.get_concordance_score(down_rankc_v1_rb)))
        plt.plot(down_rankc_v2_rb , linewidth=4, label='down_rankc_v2_rb: %.4f'%(similarity_lib.get_concordance_score(down_rankc_v2_rb)))
    if 'T-test' in robustness_individual_list:
        plt.plot(down_tt_rb , linewidth=4, label='down_tt_rb: %.4f'%(similarity_lib.get_concordance_score(down_tt_rb)))
    if 'Wilcoxon' in robustness_individual_list:
        plt.plot(down_wil_rb , linewidth=4, label='down_wil_rb: %.4f'%(similarity_lib.get_concordance_score(down_wil_rb)))

    plt.legend(fontsize=15)
    plt.xlabel("Top N protein", fontsize=15)
    plt.ylabel("Concordance", fontsize=15)
    plt.title("Robustness of the algorithm", fontsize=20)
    plt.savefig(individual_workdir +'/down_algo_roubustness.pdf', dpi=800, format='pdf', bbox_inches='tight') #指定分辨
    # plt.savefig(individual_workdir + '/down_algo_roubustness.png', dpi=800)
    plt.show()
    
    
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--parameters_path', '-p', required=True, type=click.Path(), help='parameters file')
def robustness_group(parameters_path):
    ''' Evaluation of algorithm robustness at the individual level '''
    ########################## parameters input ##############################
    params = configparser.ConfigParser()
    params.read(parameters_path)

    # workdir
    # workdir = params.get('workdir', 'workdir')
    robustness_workdir = params.get('workdir', 'ar_workdir')

    # script
    r_bpca_path = params.get('script path', 'r_bpca_path')
    r_penda_fdr_path = params.get('script path', 'r_penda_fdr_path')
    r_penda_path = params.get('script path', 'r_penda_path')
    reoa_path = params.get('script path', 'reoa_path')

    # preprocess parameters
    LOG_LABEL = eval(params.get('preprocess', 'log_label'))
    IMPUT_LABEL = eval(params.get('preprocess', 'imput_label'))
    NA_LABEL = eval(params.get('preprocess', 'na_label'))
    NA_RATE = eval(params.get('preprocess', 'na_rate'))
    NORMALIZATION = eval(params.get('preprocess', 'normalization'))

    INDEX_OLD_CHAR = eval(params.get('preprocess', 'index_old_char'))
    INDEX_NEW_CHAR = eval(params.get('preprocess', 'index_new_char'))
    COLUMNS_OLD_CHAR = eval(params.get('preprocess', 'columns_OLD_char'))
    COLUMNS_NEW_CHAR = eval(params.get('preprocess', 'columns_new_char'))

    # method parameters
    CYCLE_RANKC = eval(params.get('method parameters', 'cycle_rankc'))
    FDR = eval(params.get('method parameters', 'fdr'))
    MAX_EXCEPTION_RANKC = eval(params.get('method parameters', 'max_exception_rankc'))

    # methods comparison parameter
    data_1_path = params.get('robustness', 'data_1_path')
    tumor_cohort_1_path = params.get('robustness', 'tumor_cohort_1_path')
    normal_cohort_1_path = params.get('robustness', 'normal_cohort_1_path')

    data_2_path = params.get('robustness', 'data_2_path')
    tumor_cohort_2_path = params.get('robustness', 'tumor_cohort_2_path')
    normal_cohort_2_path = params.get('robustness', 'normal_cohort_2_path')

    log_label_list = eval(params.get('robustness', 'log_label_list'))
    na_label_list = eval(params.get('robustness', 'na_label_list'))
    imput_label_list = eval(params.get('robustness', 'imput_label_list'))

    robustness_individual_list = eval(params.get('robustness', 'robustness_individual_list'))
    robustness_group_list = eval(params.get('robustness', 'robustness_group_list'))

    N_CC = eval(params.get('robustness', 'n_cc'))
    GROUP_THRES = eval(params.get('robustness', 'group_thres'))

    # read data
    data_1 = utils.read_data(data_1_path)
    normal_cohort_1 = utils.read_normal_cohort(normal_cohort_1_path)
    tumor_cohort_1 = utils.read_tumor_cohort(tumor_cohort_1_path)
    data_2 = utils.read_data(data_2_path)
    normal_cohort_2 = utils.read_normal_cohort(normal_cohort_2_path)
    tumor_cohort_2 = utils.read_tumor_cohort(tumor_cohort_2_path)

    rd1 = raw_data.InputData(data=data_1, 
                            normal_cohort=normal_cohort_1, 
                            tumor_cohort=tumor_cohort_1, 
                            specific_protein=None, 
                            HAS_SPECIFIC_PROTEIN=False,
                            paired_data=None, 
                            paired_samples=None, 
                            HAS_PAIRED_DATA=False,  
                            INDEX_OLD_CHAR=INDEX_OLD_CHAR, 
                            INDEX_NEW_CHAR=INDEX_NEW_CHAR, 
                            COLUMNS_OLD_CHAR=COLUMNS_OLD_CHAR, 
                            COLUMNS_NEW_CHAR=COLUMNS_NEW_CHAR,
                            NORMALIZATION=NORMALIZATION, 
                            LOG_LABEL=log_label_list[0], 
                            IMPUT_LABEL=imput_label_list[0],
                            NA_LABEL=na_label_list[0], 
                            NA_RATE=NA_RATE)

    rd2 = raw_data.InputData(data=data_2, 
                            normal_cohort=normal_cohort_2, 
                            tumor_cohort=tumor_cohort_2, 
                            specific_protein=None, 
                            HAS_SPECIFIC_PROTEIN=False,
                            paired_data=None, 
                            paired_samples=None, 
                            HAS_PAIRED_DATA=False,  
                            INDEX_OLD_CHAR=INDEX_OLD_CHAR, 
                            INDEX_NEW_CHAR=INDEX_NEW_CHAR, 
                            COLUMNS_OLD_CHAR=COLUMNS_OLD_CHAR, 
                            COLUMNS_NEW_CHAR=COLUMNS_NEW_CHAR,
                            NORMALIZATION=NORMALIZATION, 
                            LOG_LABEL=log_label_list[1], 
                            IMPUT_LABEL=imput_label_list[1],
                            NA_LABEL=na_label_list[1], 
                            NA_RATE=NA_RATE)

    preprocess_workdir_1 = robustness_workdir + '/preprocess_1'

    rd1.data_preprocess(preprocess_workdir=preprocess_workdir_1,
                       r_bpca_path=r_bpca_path)

    preprocess_workdir_2 = robustness_workdir + '/preprocess_2'

    rd2.data_preprocess(preprocess_workdir=preprocess_workdir_2,
                       r_bpca_path=r_bpca_path)

    new_index = list(set(rd1.data_imput.index).intersection(set(rd2.data_imput.index)))
    
    data_1_f = rd1.data_imput.loc[new_index, :]
    data_2_f = rd2.data_imput.loc[new_index, :]

    normal_data_1 = data_1_f.loc[:, rd1.normal_cohort_adapt['normal']]
    tumor_data_1 = data_1_f.loc[:, rd1.tumor_cohort_adapt['tumor']]
    normal_data_2 = data_2_f.loc[:, rd2.normal_cohort_adapt['normal']]
    tumor_data_2 = data_2_f.loc[:, rd2.tumor_cohort_adapt['tumor']]

    # workdir
    group_workdir = robustness_workdir + '/group'

    if os.path.exists(group_workdir) == False:
        os.makedirs(group_workdir)
    os.chdir(group_workdir)

    common_rate = {}

    if 'RankComp' in robustness_group_list:
        t = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print('# Rankcom %s'%(t))
        rankc_workspace1 = group_workdir + '/rankcomp/1'
        if os.path.exists(rankc_workspace1) == True:
            shutil.rmtree(rankc_workspace1)
        os.makedirs(rankc_workspace1)
        os.chdir(rankc_workspace1)

        utils.write_normal_dat(normal_data_1, rankc_workspace1)
        utils.write_tumor_dat(tumor_data_1, rankc_workspace1)

        # run rankcomp v1 and v2
        utils.rankc_j2_qvalues(reoa_path, normal_data_1, tumor_data_1, rankc_workspace1, 
                               CYCLE_RANKC, FDR, MAX_EXCEPTION_RANKC)
        # get rankcomp result
        rankc_v1_up, rankc_v1_down, rankc_v2_up, rankc_v2_down = utils.get_rankc_j2_qvalues_result(rankc_workspace1, tumor_data_1)

        rankc_workspace2 = group_workdir + '/rankcomp/2'
        if os.path.exists(rankc_workspace2) == True:
            shutil.rmtree(rankc_workspace2)
        os.makedirs(rankc_workspace2)
        os.chdir(rankc_workspace2)

        utils.write_normal_dat(normal_data_2, rankc_workspace2)
        utils.write_tumor_dat(tumor_data_2, rankc_workspace2)

        # run rankcomp v1 and v2
        utils.rankc_j2_qvalues(reoa_path, normal_data_2, tumor_data_2, rankc_workspace2, 
                               CYCLE_RANKC, FDR, MAX_EXCEPTION_RANKC)
        # get rankcomp result
        rankc_v1_up, rankc_v1_down, rankc_v2_up, rankc_v2_down = utils.get_rankc_j2_qvalues_result(rankc_workspace2, tumor_data_2)

        rankc_v1_down_1 = pd.read_csv(rankc_workspace1 + '/rankc_v1_down.csv', index_col = 0).astype('int')
        rnakc_v2_down_1 = pd.read_csv(rankc_workspace1 + '/rankc_v2_down.csv', index_col = 0).astype('int')
        rankc_v1_up_1 = pd.read_csv(rankc_workspace1 + '/rankc_v1_up.csv', index_col = 0).astype('int')
        rnakc_v2_up_1 = pd.read_csv(rankc_workspace1 + '/rankc_v2_up.csv', index_col = 0).astype('int')
        rankc_v1_up_1[rankc_v1_up_1 == 1] = 2
        rnakc_v2_up_1[rnakc_v2_up_1 == 1] = 2
        rankc_v1_1 = rankc_v1_up_1 + rankc_v1_down_1
        rankc_v2_1 = rnakc_v2_up_1 + rnakc_v2_down_1

        rankc_v1_down_2 = pd.read_csv(rankc_workspace2 + '/rankc_v1_down.csv', index_col = 0).astype('int')
        rnakc_v2_down_2 = pd.read_csv(rankc_workspace2 + '/rankc_v2_down.csv', index_col = 0).astype('int')
        rankc_v1_up_2 = pd.read_csv(rankc_workspace2 + '/rankc_v1_up.csv', index_col = 0).astype('int')
        rnakc_v2_up_2 = pd.read_csv(rankc_workspace2 + '/rankc_v2_up.csv', index_col = 0).astype('int')
        rankc_v1_up_2[rankc_v1_up_2 == 1] = 2
        rnakc_v2_up_2[rnakc_v2_up_2 == 1] = 2
        rankc_v1_2 = rankc_v1_up_2 + rankc_v1_down_2
        rankc_v2_2 = rnakc_v2_up_2 + rnakc_v2_down_2

        rankc_v1_1_dep = robustness_lib.get_group_deps(individual_deps = rankc_v1_1, THRES = GROUP_THRES)
        rankc_v1_2_dep = robustness_lib.get_group_deps(individual_deps = rankc_v1_2, THRES = GROUP_THRES)
        rankc_v2_1_dep = robustness_lib.get_group_deps(individual_deps = rankc_v2_1, THRES = GROUP_THRES)
        rankc_v2_2_dep = robustness_lib.get_group_deps(individual_deps = rankc_v2_2, THRES = GROUP_THRES)

        common_rate['rankc_v1'] = 2 * len(set(rankc_v1_1_dep).intersection(set(rankc_v1_2_dep))) / (len(rankc_v1_1_dep) + len(rankc_v1_2_dep))
        common_rate['rankc_v2'] = 2 * len(set(rankc_v2_1_dep).intersection(set(rankc_v2_2_dep))) / (len(rankc_v2_1_dep) + len(rankc_v2_2_dep))


    if 'Penda' in robustness_group_list:
        t = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print('# Penda %s'%(t))
        penda_workspace1 = group_workdir + '/penda/1'
        if os.path.exists(penda_workspace1) == True:
            shutil.rmtree(penda_workspace1)
        os.makedirs(penda_workspace1)
        os.chdir(penda_workspace1)

        # run penda
        utils.write_penda_data(normal_data_1, tumor_data_1, penda_workspace1)
        utils.run_penda(r_penda_path)

        penda_workspace2 = group_workdir + '/penda/2'
        if os.path.exists(penda_workspace2) == True:
            shutil.rmtree(penda_workspace2)
        os.makedirs(penda_workspace2)
        os.chdir(penda_workspace2)

        # run penda
        utils.write_penda_data(normal_data_2, tumor_data_2, penda_workspace2)
        utils.run_penda(r_penda_path)

        penda_down_1 = pd.read_csv(penda_workspace1 + '/penda_down.csv', index_col = 0).astype('int')
        penda_up_1 = pd.read_csv(penda_workspace1 + '/penda_up.csv', index_col = 0).astype('int')
        penda_up_1[penda_up_1 == 1] = 2
        penda_1 = penda_up_1 + penda_down_1

        penda_down_2 = pd.read_csv(penda_workspace2 + '/penda_down.csv', index_col = 0).astype('int')
        penda_up_2 = pd.read_csv(penda_workspace2 + '/penda_up.csv', index_col = 0).astype('int')
        penda_up_2[penda_up_2 == 1] = 2
        penda_2 = penda_up_2 + penda_down_2

        penda_1_dep = robustness_lib.get_group_deps(individual_deps = penda_1, THRES = GROUP_THRES)
        penda_2_dep = robustness_lib.get_group_deps(individual_deps = penda_2, THRES = GROUP_THRES)

        common_rate['penda'] = 2 * len(set(penda_1_dep).intersection(set(penda_2_dep))) / (len(penda_1_dep) + len(penda_2_dep))


    if 'T-test' in robustness_group_list:
        t = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print('# T-test %s'%(t))
        ttest_workspace1 = group_workdir + '/ttest/1'
        if os.path.exists(ttest_workspace1) == True:
            shutil.rmtree(ttest_workspace1)
        os.makedirs(ttest_workspace1)
        os.chdir(ttest_workspace1)

        utils.write_penda_data(normal_data_1, tumor_data_1, ttest_workspace1)

        ttest_up, ttest_down = methods_lib.run_ttest(normal_run_path = './normal_run.csv',
                                                     tumor_test_path = './tumor_run.csv',
                                                     workdir = ttest_workspace1,
                                                     Q_THRES = FDR,
                                                     D_THRES = 0,
                                                     SAVE_OUT = True)

        ttest_workspace2 = group_workdir + '/ttest/2'
        if os.path.exists(ttest_workspace2) == True:
            shutil.rmtree(ttest_workspace2)
        os.makedirs(ttest_workspace2)
        os.chdir(ttest_workspace2)

        utils.write_penda_data(normal_data_2, tumor_data_2, ttest_workspace2)

        ttest_up, ttest_down = methods_lib.run_ttest(normal_run_path = './normal_run.csv',
                                                     tumor_test_path = './tumor_run.csv',
                                                     workdir = ttest_workspace2,
                                                     Q_THRES = FDR,
                                                     D_THRES = 0,
                                                     SAVE_OUT = True)

        ttest_down_1 = pd.read_csv(ttest_workspace1 + '/ttest_down.csv', index_col = 0).astype('int')
        ttest_up_1 = pd.read_csv(ttest_workspace1 + '/ttest_up.csv', index_col = 0).astype('int')
        ttest_up_1[ttest_up_1 == 1] = 2
        ttest_1 = ttest_up_1 + ttest_down_1

        ttest_down_2 = pd.read_csv(ttest_workspace2 + '/ttest_down.csv', index_col = 0).astype('int')
        ttest_up_2 = pd.read_csv(ttest_workspace2 + '/ttest_up.csv', index_col = 0).astype('int')
        ttest_up_2[ttest_up_2 == 1] = 2
        ttest_2 = ttest_up_2 + ttest_down_2

        ttest_1_dep = robustness_lib.get_group_deps(individual_deps = ttest_1, THRES = GROUP_THRES)
        ttest_2_dep = robustness_lib.get_group_deps(individual_deps = ttest_2, THRES = GROUP_THRES)

        common_rate['ttest'] = 2 * len(set(ttest_1_dep).intersection(set(ttest_2_dep))) / (len(ttest_1_dep) + len(ttest_2_dep))


    if 'Wilcoxon' in robustness_group_list:
        t = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print('# Wilcoxon %s'%(t))
        wilcox_workspace1 = group_workdir + '/wilcoxon/1'
        if os.path.exists(wilcox_workspace1) == True:
            shutil.rmtree(wilcox_workspace1)
        os.makedirs(wilcox_workspace1)
        os.chdir(wilcox_workspace1)

        utils.write_penda_data(normal_data_1, tumor_data_1, wilcox_workspace1)

        wilcox_up, wilcox_down = methods_lib.run_wilcox(normal_run_path = './normal_run.csv', 
                                                        tumor_test_path = './tumor_run.csv',
                                                        workdir = wilcox_workspace1, 
                                                        Q_THRES = FDR, D_THRES = 0,
                                                        SAVE_OUT = True)


        wilcox_workspace2 = group_workdir + '/wilcoxon/2'
        if os.path.exists(wilcox_workspace2) == True:
            shutil.rmtree(wilcox_workspace2)
        os.makedirs(wilcox_workspace2)
        os.chdir(wilcox_workspace2)

        utils.write_penda_data(normal_data_2, tumor_data_2, wilcox_workspace2)

        wilcox_up, wilcox_down = methods_lib.run_wilcox(normal_run_path = './normal_run.csv', 
                                                        tumor_test_path = './tumor_run.csv',
                                                        workdir = wilcox_workspace2, 
                                                        Q_THRES = FDR, D_THRES = 0,
                                                        SAVE_OUT = True)

        wilcox_down_1 = pd.read_csv(wilcox_workspace1 + '/wilcox_down.csv', index_col = 0).astype('int')
        wilcox_up_1 = pd.read_csv(wilcox_workspace1 + '/wilcox_up.csv', index_col = 0).astype('int')
        wilcox_up_1[wilcox_up_1 == 1] = 2
        wilcox_1 = wilcox_up_1 + wilcox_down_1

        wilcox_down_2 = pd.read_csv(wilcox_workspace2 + '/wilcox_down.csv', index_col = 0).astype('int')
        wilcox_up_2 = pd.read_csv(wilcox_workspace2 + '/wilcox_up.csv', index_col = 0).astype('int')
        wilcox_up_2[wilcox_up_2 == 1] = 2
        wilcox_2 = wilcox_up_2 + wilcox_down_2

        wilcox_1_dep = robustness_lib.get_group_deps(individual_deps = wilcox_1, THRES = GROUP_THRES)
        wilcox_2_dep = robustness_lib.get_group_deps(individual_deps = wilcox_2, THRES = GROUP_THRES)

        common_rate['wilcoxon'] = 2 * len(set(wilcox_1_dep).intersection(set(wilcox_2_dep))) / (len(wilcox_1_dep) + len(wilcox_2_dep))


    if 'Peng method' in robustness_group_list:
        t = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print('# Peng method %s'%(t))
        peng_workspace1 = group_workdir + '/peng/1'
        if os.path.exists(peng_workspace1) == True:
            shutil.rmtree(peng_workspace1)
        os.makedirs(peng_workspace1)
        os.chdir(peng_workspace1)

        # run penda
        peng_up, peng_down = methods_lib.run_peng_method(normal_data_1, tumor_data_1)
        peng_up.to_csv('./peng_up.csv')
        peng_down.to_csv('./peng_down.csv')

        peng_workspace2 = group_workdir + '/peng/2'
        if os.path.exists(peng_workspace2) == True:
            shutil.rmtree(peng_workspace2)
        os.makedirs(peng_workspace2)
        os.chdir(peng_workspace2)

        # run penda
        peng_up, peng_down = methods_lib.run_peng_method(normal_data_2, tumor_data_2)
        peng_up.to_csv('./peng_up.csv')
        peng_down.to_csv('./peng_down.csv')

        peng_down_1 = pd.read_csv(peng_workspace1 + '/peng_down.csv', index_col = 0).astype('int')
        peng_up_1 = pd.read_csv(peng_workspace1 + '/peng_up.csv', index_col = 0).astype('int')
        peng_up_1[peng_up_1 == 1] =2
        peng_1 = peng_up_1 + peng_down_1

        peng_down_2 = pd.read_csv(peng_workspace2 + '/peng_down.csv', index_col = 0).astype('int')
        peng_up_2 = pd.read_csv(peng_workspace2 + '/peng_up.csv', index_col = 0).astype('int')
        peng_up_2[peng_up_2 == 1] =2
        peng_2 = peng_up_2 + peng_down_2

        peng_1_dep = robustness_lib.get_group_deps(individual_deps = peng_1, THRES = GROUP_THRES)
        peng_2_dep = robustness_lib.get_group_deps(individual_deps = peng_2, THRES = GROUP_THRES)

        common_rate['peng'] = 2 * len(set(peng_1_dep).intersection(set(peng_2_dep))) / (len(peng_1_dep) + len(peng_2_dep))


    if 'Quantile' in robustness_group_list:
        t = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print('# quantile %s'%(t))
        quantile_workspace1 = group_workdir + '/quantile/1'
        if os.path.exists(quantile_workspace1) == True:
            shutil.rmtree(quantile_workspace1)
        os.makedirs(quantile_workspace1)
        os.chdir(quantile_workspace1)

        utils.write_penda_data(normal_data_1, tumor_data_1, quantile_workspace1)

        quantile_up, quantile_down = methods_lib.quantile_dep(normal_data_path = './normal_run.csv',
                                                              tumor_data_path = './tumor_run.csv',
                                                              workdir = quantile_workspace1, 
                                                              quantile = 0.05, 
                                                              factor = 1.2,
                                                              SAVE_OUT = True)


        quantile_workspace2 = group_workdir + '/quantile/2'
        if os.path.exists(quantile_workspace2) == True:
            shutil.rmtree(quantile_workspace2)
        os.makedirs(quantile_workspace2)
        os.chdir(quantile_workspace2)

        utils.write_penda_data(normal_data_2, tumor_data_2, quantile_workspace2)

        quantile_up, quantile_down = methods_lib.quantile_dep(normal_data_path = './normal_run.csv',
                                                              tumor_data_path = './tumor_run.csv',
                                                              workdir = quantile_workspace2, 
                                                              quantile = 0.05, 
                                                              factor = 1.2,
                                                              SAVE_OUT = True)

        quantile_down_1 = pd.read_csv(quantile_workspace1 + '/quantile_down.csv', index_col = 0).astype('int')
        quantile_up_1 = pd.read_csv(quantile_workspace1 + '/quantile_up.csv', index_col = 0).astype('int')
        quantile_up_1[quantile_up_1 == 1] = 2
        quantile_1 = quantile_up_1 + quantile_down_1

        quantile_down_2 = pd.read_csv(quantile_workspace2 + '/quantile_down.csv', index_col = 0).astype('int')
        quantile_up_2 = pd.read_csv(quantile_workspace2 + '/quantile_up.csv', index_col = 0).astype('int')
        quantile_up_2[quantile_up_2 == 1] = 2
        quantile_2 = quantile_up_2 + quantile_down_2

        quantile_1_dep = robustness_lib.get_group_deps(individual_deps = quantile_1, THRES = GROUP_THRES)
        quantile_2_dep = robustness_lib.get_group_deps(individual_deps = quantile_2, THRES = GROUP_THRES)

        common_rate['quantile'] = 2 * len(set(quantile_1_dep).intersection(set(quantile_2_dep))) / (len(quantile_1_dep) + len(quantile_2_dep))


    from matplotlib.backends.backend_pdf import PdfPages
    os.chdir(group_workdir)
    with PdfPages('ar_lung_common_deps.pdf') as pdf:
        if 'RankComp' in robustness_group_list:
            plt.figure(figsize=(8, 8))#控
            g=venn2(subsets = [set(rankc_v1_1_dep),set(rankc_v1_2_dep)], 
                    set_labels = ('1-lung', '5-lung'), 
                    set_colors=("#098154","#c72e29"),
                    alpha=0.6,
                    normalize_to=1.0,
                   )

            g=venn2_circles(subsets = [set(rankc_v1_1_dep),set(rankc_v1_2_dep)], 
                    linestyle='--', linewidth=2.0, color="black"
                   )
            plt.title('RankComp v1', fontsize=25)
            pdf.savefig()  # saves the current figure into a pdf page
            plt.close()


            plt.figure(figsize=(8, 8))#控
            g=venn2(subsets = [set(rankc_v2_1_dep),set(rankc_v2_2_dep)], 
                    set_labels = ('1-lung', '5-lung'), 
                    set_colors=("#098154","#c72e29"),
                    alpha=0.6,
                    normalize_to=1.0,
                   )

            g=venn2_circles(subsets = [set(rankc_v2_1_dep),set(rankc_v2_2_dep)], 
                    linestyle='--', linewidth=2.0, color="black"
                   )
            plt.title('RankComp v2', fontsize=25)
            pdf.savefig()  # saves the current figure into a pdf page
            plt.close()

        if 'Penda' in robustness_group_list:
            plt.figure(figsize=(8, 8))
            g=venn2(subsets = [set(penda_1_dep),set(penda_2_dep)], 
                    set_labels = ('1-lung', '5-lung'), 
                    set_colors=("#098154","#c72e29"),
                    alpha=0.6,
                    normalize_to=1.0,
                   )

            g=venn2_circles(subsets = [set(penda_1_dep),set(penda_2_dep)], 
                    linestyle='--', linewidth=2.0, color="black"
                   )
            plt.title('Penda', fontsize=25)
            pdf.savefig()  # saves the current figure into a pdf page
            plt.close()

        if 'Peng method' in robustness_group_list:
            plt.figure(figsize=(8, 8))
            g=venn2(subsets = [set(peng_1_dep),set(peng_2_dep)],
                    set_labels = ('1-lung', '5-lung'), 
                    set_colors=("#098154","#c72e29"),
                    alpha=0.6,
                    normalize_to=1.0,
                   )

            g=venn2_circles(subsets = [set(peng_1_dep),set(peng_2_dep)], 
                    linestyle='--', linewidth=2.0, color="black"
                   )
            plt.title('Peng method', fontsize=25)
            pdf.savefig()  # saves the current figure into a pdf page
            plt.close()

        if 'Quantile' in robustness_group_list:
            plt.figure(figsize=(8, 8))
            g=venn2(subsets = [set(quantile_1_dep),set(quantile_2_dep)], 
                    set_labels = ('1-lung', '5-lung'), 
                    set_colors=("#098154","#c72e29"),
                    alpha=0.6,
                    normalize_to=1.0,
                   )

            g=venn2_circles(subsets = [set(quantile_1_dep),set(quantile_2_dep)], 
                    linestyle='--', linewidth=2.0, color="black"
                   )
            plt.title('Quantile', fontsize=25)
            pdf.savefig()  # saves the current figure into a pdf page
            plt.close()

        if 'T-test' in robustness_group_list:
            plt.figure(figsize=(8, 8))
            g=venn2(subsets = [set(ttest_1_dep),set(ttest_2_dep)],
                    set_labels = ('1-lung', '5-lung'), 
                    set_colors=("#098154","#c72e29"),
                    alpha=0.6,
                    normalize_to=1.0,
                   )

            g=venn2_circles(subsets = [set(ttest_1_dep),set(ttest_2_dep)], 
                    linestyle='--', linewidth=2.0, color="black"
                   )
            plt.title('T-test', fontsize=25)
            pdf.savefig()  # saves the current figure into a pdf page
            plt.close()

        if 'Wilcoxon' in robustness_group_list:
            plt.figure(figsize=(8, 8))
            g=venn2(subsets = [set(wilcox_1_dep),set(wilcox_2_dep)],
                    set_labels = ('1-lung', '5-lung'), 
                    set_colors=("#098154","#c72e29"),
                    alpha=0.6,
                    normalize_to=1.0,
                   )

            g=venn2_circles(subsets = [set(wilcox_1_dep),set(wilcox_2_dep)], 
                    linestyle='--', linewidth=2.0, color="black"
                   )
            plt.title('Wilcoxon', fontsize=25)
            pdf.savefig()  # saves the current figure into a pdf page
            plt.close()



    common_rate_pd = pd.DataFrame({'methods':common_rate.keys(), 'common_rate': common_rate.values()})

    plt.figure(figsize=(12,8))


    plt.xlabel("Methods", fontsize=20)
    plt.ylabel('Common rate', fontsize=20)

    plt.plot(common_rate_pd['methods'], common_rate_pd['common_rate'], alpha=1, ms=10, lw=3, marker='o', c='b')


    sns.despine(left=True, bottom=True)   
    plt.tick_params(labelsize=15)
    plt.legend(loc=[1.02, 0.5])

    # plt.xticks(np.array(x_x)+ 3*bar_width/2, tick_label)

    plt.savefig(group_workdir + '/ar_common_rate.pdf', dpi=800, bbox_inches='tight')
    plt.show()
    
    
    
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--parameters_path', '-p', required=True, type=click.Path(), help='parameters file')
def kegg(parameters_path):
    ''' Pathway enrichment analysis '''
    ########################## parameters input ##############################
    params = configparser.ConfigParser()
    params.read(parameters_path)

    # workdir
    # workdir = params.get('workdir', 'workdir')
    kegg_workdir = params.get('workdir', 'kegg_workdir')

    # script
    r_bpca_path = params.get('script path', 'r_bpca_path')
    r_penda_fdr_path = params.get('script path', 'r_penda_fdr_path')
    r_penda_path = params.get('script path', 'r_penda_path')
    reoa_path = params.get('script path', 'reoa_path')

    # preprocess parameters
    LOG_LABEL = eval(params.get('preprocess', 'log_label'))
    IMPUT_LABEL = eval(params.get('preprocess', 'imput_label'))
    NA_LABEL = eval(params.get('preprocess', 'na_label'))
    NA_RATE = eval(params.get('preprocess', 'na_rate'))
    NORMALIZATION = eval(params.get('preprocess', 'normalization'))

    INDEX_OLD_CHAR = eval(params.get('preprocess', 'index_old_char'))
    INDEX_NEW_CHAR = eval(params.get('preprocess', 'index_new_char'))
    COLUMNS_OLD_CHAR = eval(params.get('preprocess', 'columns_OLD_char'))
    COLUMNS_NEW_CHAR = eval(params.get('preprocess', 'columns_new_char'))

    # method parameters
    CYCLE_RANKC = eval(params.get('method parameters', 'cycle_rankc'))
    FDR = eval(params.get('method parameters', 'fdr'))
    MAX_EXCEPTION_RANKC = eval(params.get('method parameters', 'max_exception_rankc'))

    # methods comparison parameter
    data_path = params.get('kegg', 'data_path')
    normal_cohort_path = params.get('kegg', 'normal_cohort_path')
    tumor_cohort_path = params.get('kegg', 'tumor_cohort_path')
    KEGG_METHODS_LIST = eval(params.get('kegg', 'kegg_methods_list'))
    GROUP_THRES = eval(params.get('kegg', 'group_thres'))
    r_kegg_path = eval(params.get('kegg', 'r_kegg_path'))
    R_DIR = eval(params.get('kegg', 'r_dir'))
    Q_VALUE_THRES = eval(params.get('kegg', 'qvalue_thres'))


    # read data
    data = utils.read_data(data_path)
    normal_cohort = utils.read_normal_cohort(normal_cohort_path)
    tumor_cohort = utils.read_tumor_cohort(tumor_cohort_path)

    rd = raw_data.InputData(data=data, 
                            normal_cohort=normal_cohort, 
                            tumor_cohort=tumor_cohort, 
                            specific_protein=None, 
                            HAS_SPECIFIC_PROTEIN=False,
                            paired_data=None, 
                            paired_samples=None, 
                            HAS_PAIRED_DATA=False,  
                            INDEX_OLD_CHAR=INDEX_OLD_CHAR, 
                            INDEX_NEW_CHAR=INDEX_NEW_CHAR, 
                            COLUMNS_OLD_CHAR=COLUMNS_OLD_CHAR, 
                            COLUMNS_NEW_CHAR=COLUMNS_NEW_CHAR,
                            NORMALIZATION=NORMALIZATION, 
                            LOG_LABEL=LOG_LABEL, 
                            IMPUT_LABEL=IMPUT_LABEL,
                            NA_LABEL=NA_LABEL, 
                            NA_RATE=NA_RATE)

    preprocess_workdir = kegg_workdir + '/preprocess'

    rd.data_preprocess(preprocess_workdir=preprocess_workdir,
                       r_bpca_path=r_bpca_path)

    mc = methods_comp.methodsComp(rd = rd, 
                                  r_penda_path = r_penda_path, 
                                  r_penda_fdr_path = r_penda_fdr_path,
                                  reoa_path = reoa_path,
                                  CYCLE_RANKC = CYCLE_RANKC, 
                                  FDR = FDR, 
                                  MAX_EXCEPTION_RANKC = MAX_EXCEPTION_RANKC, 
                                  PLOT_METHODS_COMPARE_RESULT = False,
                                  METHODS_LIST = KEGG_METHODS_LIST)
    methods_workdir = kegg_workdir + '/methods_result'
    mc.run_methodsComp(mc_workdir = methods_workdir, method_comp_label=False)

    kegg_result_workdir = kegg_workdir + '/kegg_result'

    kegg_result = kegg_lib.get_kegg_result_dataset(workdir = kegg_result_workdir,
                                                     r_kegg_path = r_kegg_path,
                                                     R_dir = R_DIR,
                                                     GROUP_THRES = GROUP_THRES,
                                                     result_dir = methods_workdir,
                                                     kegg_methods_list = KEGG_METHODS_LIST,
                                                     Q_VALUE_THRES = Q_VALUE_THRES)


    dt = pd.pivot_table(kegg_result, index=['Description'], columns=['methods'], values=['qvalue'], aggfunc=[max])
    dt = -dt.apply(np.log10)

    max_value, max_idx = dt.stack().max(), dt.stack().idxmax()

    new_col = dt.columns.levels[2]
    dt.columns = new_col

    dt.sort_index(axis=0, ascending=False, inplace=True)  
    dt.sort_index(axis=1, ascending=False, inplace=True)  

    
    plt.figure(figsize=(8, 13))
    g = sns.heatmap(dt, vmin=0.0, fmt='.2g', cmap='Blues', cbar=True, cbar_kws={'label': '-log10(p)'})
    plt.ylabel("Description",fontsize=15)
    plt.xlabel("Methods",fontsize=15)

    plt.savefig(kegg_workdir + '/heatmap_individual.pdf', format='pdf', bbox_inches='tight') 
    plt.show()


    
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--parameters_path', '-p', required=True, type=click.Path(), help='parameters file')
def survival(parameters_path):
    ''' survival analysis '''
    
    rpy2.robjects.r('''
        library('survival')
        library('survminer')
        library('SuperExactTest')

        cox_analysis = function(clinical_select, specific_prot, 
             wdir,
             p_thres = 0.05
             )
    {
      specific_prot = as.character(specific_prot$specific_prot)
      specific_prot_fix = c()

      clinical_select <<- clinical_select
      for(i in specific_prot){
        specific_prot_fix<-c(specific_prot_fix,gsub('-','.',i))
      }


      covariates <- specific_prot_fix
      univ_formulas <- sapply(covariates,
                              function(x) as.formula(paste('Surv(time, status)~', x)))

      univ_models <- lapply( univ_formulas, function(x){coxph(x, data = clinical_select)})
      # Extract data 
      univ_results <- lapply(univ_models,
                             function(x){ 
                               x <- summary(x)
                               p.value<-signif(x$wald["pvalue"], digits=2)
                               wald.test<-signif(x$wald["test"], digits=2)
                               beta<-signif(x$coef[1], digits=2);#coeficient beta
                               HR <-signif(x$coef[2], digits=2);#exp(beta)
                               HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                               HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                               HR <- paste0(HR, " (", 
                                            HR.confint.lower, "-", HR.confint.upper, ")")
                               res<-c(beta, HR, wald.test, p.value)
                               names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                             "p.value")
                               return(res)
                               #return(exp(cbind(coef(x),confint(x))))
                             })
      res <- t(as.data.frame(univ_results, check.names = FALSE))
      result = as.data.frame(res)
      result_f = result[complete.cases(result),]

      result_f$p.value = as.numeric(as.character(result_f$p.value))

      write.csv(result_f, 
                paste(wdir,'/cox_result.csv', sep=""),
                sep=',', col.names = T, row.names = T)


      result_select = result_f[result_f$p.value < p_thres,]

      write.csv(result_select, 
                paste(wdir,'/cox_result_filter.csv', sep=""),
                sep=',', col.names = T, row.names = T)

      rownames_rs = rownames(result_select)

      pdf(paste(wdir,'/prot_cox.pdf', sep=""))
      for(n in seq_along(rownames_rs)){

        prot <<- rownames_rs[n]

        fit = survfit(Surv(time, status)~get(prot), data = clinical_select)

        #  pdf(paste(get('prot'),".pdf",sep = ""),width=15,height=10)
        p = ggsurvplot(fit, data = clinical_select,
                        conf.int = TRUE,
                        pval = TRUE,
                        fun = "pct",
                        risk.table = T,
                        size = 2,
                        linetype = "strata",
                        palette = c("#E7B800",
                                    "#2E9FDF"),
                        legend = "bottom",
                        legend.title = prot,
                        legend.labs = c("UNDEP",
                                        "DEP"))
        print(p)


      }  
      dev.off()

      result_select <<- data.frame(result_select)

      return(result_select)

    }

        ''')

    
    ########################## parameters input ##############################
    params = configparser.ConfigParser()
    params.read(parameters_path)

    # workdir
    # workdir = params.get('workdir', 'workdir')
    surv_workdir = params.get('workdir', 'survival_workdir')

    # script
    r_bpca_path = params.get('script path', 'r_bpca_path')
    r_penda_fdr_path = params.get('script path', 'r_penda_fdr_path')
    r_penda_path = params.get('script path', 'r_penda_path')
    reoa_path = params.get('script path', 'reoa_path')

    # preprocess parameters
    LOG_LABEL = eval(params.get('preprocess', 'log_label'))
    IMPUT_LABEL = eval(params.get('preprocess', 'imput_label'))
    NA_LABEL = eval(params.get('preprocess', 'na_label'))
    NA_RATE = eval(params.get('preprocess', 'na_rate'))
    NORMALIZATION = eval(params.get('preprocess', 'normalization'))

    INDEX_OLD_CHAR = eval(params.get('preprocess', 'index_old_char'))
    INDEX_NEW_CHAR = eval(params.get('preprocess', 'index_new_char'))
    COLUMNS_OLD_CHAR = eval(params.get('preprocess', 'columns_OLD_char'))
    COLUMNS_NEW_CHAR = eval(params.get('preprocess', 'columns_new_char'))

    # method parameters
    CYCLE_RANKC = eval(params.get('method parameters', 'cycle_rankc'))
    FDR = eval(params.get('method parameters', 'fdr'))
    MAX_EXCEPTION_RANKC = eval(params.get('method parameters', 'max_exception_rankc'))

    # survival parameters
    data_path = params.get('survival analysis', 'data_path')
    normal_cohort_path = params.get('survival analysis', 'normal_cohort_path')
    tumor_cohort_path = params.get('survival analysis', 'tumor_cohort_path')
    clinial_path = params.get('survival analysis', 'clinial_path')

    SURVIVAL_METHODS_LIST = eval(params.get('survival analysis', 'survival_methods_list'))
    SURV_THRES = eval(params.get('survival analysis', 'surv_thres'))

    # read data
    data = utils.read_data(data_path)
    normal_cohort = utils.read_normal_cohort(normal_cohort_path)
    tumor_cohort = utils.read_tumor_cohort(tumor_cohort_path)
    clinial = pd.read_csv(clinial_path, index_col=0)

    clinial = survival_lib.col_adapt(clinial, columns = 'samples', characters = COLUMNS_OLD_CHAR, new_char = COLUMNS_NEW_CHAR)

    rd = raw_data.InputData(data=data, 
                            normal_cohort=normal_cohort, 
                            tumor_cohort=tumor_cohort, 
                            specific_protein=None, 
                            HAS_SPECIFIC_PROTEIN=False,
                            paired_data=None, 
                            paired_samples=None, 
                            HAS_PAIRED_DATA=False,  
                            INDEX_OLD_CHAR=INDEX_OLD_CHAR, 
                            INDEX_NEW_CHAR=INDEX_NEW_CHAR, 
                            COLUMNS_OLD_CHAR=COLUMNS_OLD_CHAR, 
                            COLUMNS_NEW_CHAR=COLUMNS_NEW_CHAR,
                            NORMALIZATION=NORMALIZATION, 
                            LOG_LABEL=LOG_LABEL, 
                            IMPUT_LABEL=IMPUT_LABEL,
                            NA_LABEL=NA_LABEL, 
                            NA_RATE=NA_RATE)

    preprocess_workdir = surv_workdir + '/preprocess'

    rd.data_preprocess(preprocess_workdir=preprocess_workdir,
                       r_bpca_path=r_bpca_path)

    mc = methods_comp.methodsComp(rd = rd, 
                                  r_penda_path = r_penda_path, 
                                  r_penda_fdr_path = r_penda_fdr_path,
                                  reoa_path = reoa_path,
                                  CYCLE_RANKC = CYCLE_RANKC, 
                                  FDR = FDR, 
                                  MAX_EXCEPTION_RANKC = MAX_EXCEPTION_RANKC, 
                                  PLOT_METHODS_COMPARE_RESULT = False,
                                  METHODS_LIST = SURVIVAL_METHODS_LIST)
    methods_workdir = surv_workdir + '/methods_result'
    mc.run_methodsComp(mc_workdir = methods_workdir, method_comp_label=False)

    survival_workdir = surv_workdir + '/survive_analysis'

    if os.path.exists(survival_workdir) == False:
        os.makedirs(survival_workdir)
    os.chdir(survival_workdir)

    merge_data = []

    if 'Penda' in SURVIVAL_METHODS_LIST:
        ### penda
        penda_up_path = surv_workdir + '/methods_result/penda/penda_up.csv'
        penda_down_path = surv_workdir + '/methods_result/penda/penda_down.csv'

        _penda_up = pd.read_csv(penda_up_path, index_col=0)
        _penda_down = pd.read_csv(penda_down_path, index_col=0)

        penda_dep = _penda_up | _penda_down
        penda_dep[penda_dep] = 1

        penda_cox_workdir = survival_workdir + '/penda'
        clinial_select = clinial.set_index('samples').loc[penda_dep.columns.tolist(), :].reset_index()
        penda_cox_result, penda_cox_rate = survival_lib.cox_analysis(_dep = penda_dep, 
                                                                      clinial_select = clinial_select, 
                                                                      THRES_SURV = SURV_THRES, 
                                                                      cox_run_workdir = penda_cox_workdir)
        penda_cox_result['methods'] = 'penda'
        if penda_cox_result.shape[0] > 0:
            merge_data.append(penda_cox_result)

    if 'Quantile' in SURVIVAL_METHODS_LIST:
        ### quantile
        quantile_up_path = surv_workdir + '/methods_result/quantile/quantile_up.csv'
        quantile_down_path = surv_workdir + '/methods_result/quantile/quantile_down.csv'

        _quantile_up = pd.read_csv(quantile_up_path, index_col=0)
        _quantile_down = pd.read_csv(quantile_down_path, index_col=0)

        _quantile_up = _quantile_up.astype('int')
        _quantile_down = _quantile_down.astype('int')
        quantile_dep = _quantile_up | _quantile_down

        quantile_cox_workdir = survival_workdir + '/quantile'
        clinial_select = clinial.set_index('samples').loc[quantile_dep.columns.tolist(), :].reset_index()
        quantile_cox_result, quantile_cox_rate = survival_lib.cox_analysis(_dep = quantile_dep, 
                                                                           clinial_select = clinial_select, 
                                                                           THRES_SURV = SURV_THRES, 
                                                                           cox_run_workdir = quantile_cox_workdir)
        quantile_cox_result['methods'] = 'quantile'
        if quantile_cox_result.shape[0] > 0:
            merge_data.append(quantile_cox_result)

    if 'RankComp' in SURVIVAL_METHODS_LIST:
        ### Rankcomp v1
        rankc_v1_up_path = surv_workdir + '/methods_result/rankcomp/rankc_v1_up.csv'
        rankc_v1_down_path = surv_workdir + '/methods_result/rankcomp/rankc_v1_down.csv'

        _rankc_v1_up = pd.read_csv(rankc_v1_up_path, index_col=0)
        _rankc_v1_down = pd.read_csv(rankc_v1_down_path, index_col=0)

        rankc_v1_dep = _rankc_v1_up | _rankc_v1_down
        rankc_v1_dep[rankc_v1_dep] = 1

        rankc_v1_cox_workdir = survival_workdir + '/rankc_v1'
        clinial_select = clinial.set_index('samples').loc[rankc_v1_dep.columns.tolist(), :].reset_index()
        rankc_v1_cox_result, rankc_v1_cox_rate = survival_lib.cox_analysis(_dep = rankc_v1_dep, 
                                                                          clinial_select = clinial_select, 
                                                                          THRES_SURV = SURV_THRES, 
                                                                          cox_run_workdir = rankc_v1_cox_workdir)
        rankc_v1_cox_result['methods'] = 'rankc_v1'
        if rankc_v1_cox_result.shape[0] > 0:
            merge_data.append(rankc_v1_cox_result)

        rankc_v2_up_path = surv_workdir + '/methods_result/rankcomp/rankc_v2_up.csv'
        rankc_v2_down_path = surv_workdir + '/methods_result/rankcomp/rankc_v2_down.csv'

        _rankc_v2_up = pd.read_csv(rankc_v2_up_path, index_col=0)
        _rankc_v2_down = pd.read_csv(rankc_v2_down_path, index_col=0)

        rankc_v2_dep = _rankc_v2_up | _rankc_v2_down
        rankc_v2_dep[rankc_v2_dep] = 1

        rankc_v2_cox_workdir = survival_workdir + '/rankc_v2'
        clinial_select = clinial.set_index('samples').loc[rankc_v2_dep.columns.tolist(), :].reset_index()
        rankc_v2_cox_result, rankc_v2_cox_rate = survival_lib.cox_analysis(_dep = rankc_v2_dep, 
                                                                          clinial_select = clinial_select, 
                                                                          THRES_SURV = SURV_THRES, 
                                                                          cox_run_workdir = rankc_v2_cox_workdir)
        rankc_v2_cox_result['methods'] = 'rankc_v2'
        if rankc_v2_cox_result.shape[0] > 0:
            merge_data.append(rankc_v2_cox_result)

    if 'Peng method' in SURVIVAL_METHODS_LIST:
        peng_up_path = surv_workdir + '/methods_result/peng_method/peng_up.csv'
        peng_down_path = surv_workdir + '/methods_result/peng_method/peng_down.csv'

        _peng_up = pd.read_csv(peng_up_path, index_col=0)
        _peng_down = pd.read_csv(peng_down_path, index_col=0)

        peng_dep = _peng_up | _peng_down
        peng_dep[peng_dep] = 1

        peng_cox_workdir = survival_workdir + '/peng'
        clinial_select = clinial.set_index('samples').loc[peng_dep.columns.tolist(), :].reset_index()
        peng_cox_result, peng_cox_rate = survival_lib.cox_analysis(_dep = peng_dep, 
                                                                  clinial_select = clinial_select, 
                                                                  THRES_SURV = SURV_THRES, 
                                                                  cox_run_workdir = peng_cox_workdir)
        peng_cox_result['methods'] = 'peng'
        if peng_cox_result.shape[0] > 0:
            merge_data.append(peng_cox_result)

    cox_result = pd.concat(merge_data, axis=0)
    cox_result.reset_index(inplace=True)
    cox_result.columns = ['protein', 'bata', 'HR', 'wald.test', 'p_value', 'methods']
    cox_result.to_csv(surv_workdir + '/cox_result.csv')

    methods_list = cox_result['methods'].unique()
    n_protein = []
    for ml in methods_list:
        n_protein.append(cox_result[cox_result['methods'] == ml].shape[0])

    hr_color = []
    for hr in cox_result['HR']:
        _hr = eval(hr.split(' ')[0])
        if _hr > 1:
            hr_color.append('red')
        else:
            hr_color.append('blue')

    import matplotlib.pyplot as plt
    import numpy as np
    plt.figure(figsize=(12,8))

    size = list(1 / cox_result['p_value'].values)
    x = cox_result['methods']
    y = cox_result['p_value']

    # colors = np.random.rand(len(tem)) 

    fig = plt.figure(figsize=(12,8))

    ax1 = fig.add_subplot(111)
    ax1.bar(methods_list, n_protein, alpha=0.5, color='grey')
    ax1.set_ylabel('Number of DE protein', fontsize='20')
    ax1.set_xlabel("Methods", fontsize=20)

    ax2 = ax1.twinx()   
    ax2.scatter(x, y, s=size, c=hr_color, alpha=0.5)  
    ax2.set_ylabel('P values', fontsize=20)
    sns.despine(left=True, bottom=True)   
    ax2.tick_params(labelsize=15)
    plt.savefig(surv_workdir + '/survival.pdf', format='pdf', bbox_inches='tight') 
    plt.show()



    
    
    
    

cli.add_command(stable)
cli.add_command(comparison)
cli.add_command(evaluation)
cli.add_command(similarity)
cli.add_command(type1error)
cli.add_command(robustness_individual)
cli.add_command(robustness_group)
cli.add_command(kegg)
cli.add_command(survival)

if __name__ == '__main__':

    cli()
