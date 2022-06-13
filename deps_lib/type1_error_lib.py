import copy
import os
import datetime
import random

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from scipy import stats

from deps_lib import methods_lib, penda_pro, utils

from deps_lib import utils
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector

r_stats = importr('stats')


def data_col_adapt(data, characters = [' ', '-'], new_char = '.'):
    """
    由于后续需要用到R语言工具，对data列名中不合法的字符进行替换。
    :param data: protein abundance
    :param character: characters that need to be replaced
    :param new_char: new character
    :return data: renamed data
    """
    for ch in characters:
        new_col = []
        for col in data.columns:
            new_col.append(col.replace(ch, new_char))
        data.columns = new_col
        
    return data

def data_index_adapt(data, characters = [' ', '-'], new_char = '.'):
    """
    由于后续需要用到R语言工具，对data行名中不合法的字符进行替换。
    :param data: protein abundance
    :param character: characters that need to be replaced
    :param new_char: new character
    :return data: renamed data
    """
    for ch in characters:
        new_index = []
        for idx in data.index:
            new_index.append(idx.replace(ch, new_char))
        data.index = new_index
        
    return data

def data_col_index_adapt(data, col_characters = [' ', '-'], index_characters = ['-'], new_char = '.'):
    """
    由于后续需要用到R语言工具，对data行列名中不合法的字符进行替换。
    :param data: protein abundance
    :param col_characters
    :param index_characters
    :param new_char: new character
    :return data: renamed data
    """
    data_tmp1 = data_col_adapt(data, characters = col_characters, new_char = new_char)
    data_tmp2 = data_index_adapt(data_tmp1, characters = index_characters, new_char = new_char)
    return data_tmp2

def data_index_recover(data, index_characters = '-', new_char = '.'):
    """
    由于 下游分析的需要，需要将转换后的行名恢复符合国际命名规则的蛋白质名称
    :param data
    :param index_characters
    :param new_char
    :return data
    """
    new_index = []
    for idx in data.index:
        new_index.append(idx.replace(new_char, index_characters))
    data.index = new_index
        
    return data

def col_adapt(samples_col, characters = [' ', '-'], new_char = '.'):
    """
    由于后续需要用到R语言工具，对samples样品名中不合法的字符进行替换。
    :param samples_col: sample columns
    :param character: characters that need to be replaced
    :param new_char: new character
    :return samples_col: renamed columns
    """
    for ch in characters:
        for idx in range(samples_col.shape[0]):
            for col in range(samples_col.shape[1]):
                samples_col.iloc[idx, col] = samples_col.iloc[idx, col].replace(ch, new_char)
                
    return samples_col

def data_preprocess(data_path, normal_col_path, tumor_col_path, r_bpca_path, workdir, data_label, 
                    LOG_LABEL=False, NA_LABEL = '', NA_RATE=0.3, INDEX_OLD_CHAR = ['-', ' '],
                    INDEX_NEW_CHAR = '.', COLUMNS_OLD_CHAR = ['-', ' '],
                    COLUMNS_NEW_CHAR = '.'):
    # 数据读入
    data = pd.read_csv(data_path, index_col=0)
    normal_col = pd.read_csv(normal_col_path, sep='\t')
    tumor_col = pd.read_csv(tumor_col_path, sep='\t')

    # 将代表NA的标记修改为np.nan
    if NA_LABEL == '':
        pass
    else:
        data = data.replace(NA_LABEL, np.nan)

    # 质量控制
    data_tumor = data.loc[:, tumor_col['tumor']]
    data_normal = data.loc[:, normal_col['normal']]
    data_tumor_fna = data_tumor[data_tumor.isna().sum(axis=1) <= NA_RATE * data_tumor.shape[1]]
    data_normal_fna = data_normal[data_normal.isna().sum(axis=1) <= NA_RATE * data_normal.shape[1]]
    data_fna_index = set(data_tumor_fna.index).intersection(set(data_normal_fna.index))
    data_fna = data.loc[data_fna_index,:]

    # LT
    if LOG_LABEL:
        data_fna_lt = copy.deepcopy(data_fna)
        for col in data_fna.columns:
            data_fna_lt[col] = np.log2(data_fna[col].values)
    else:
        data_fna_lt = copy.deepcopy(data_fna)

    # Z-score normalization
    data_fna_lt_zn = pd.DataFrame(stats.zscore(data_fna_lt, axis=0, nan_policy='omit'), index=data_fna_lt.index, columns=data_fna_lt.columns)

    # adapt columns
    data_fna_lt_zn_da = data_col_index_adapt(data = data_fna_lt_zn,
                                            col_characters = COLUMNS_OLD_CHAR,
                                            index_characters = INDEX_OLD_CHAR, 
                                            new_char = INDEX_NEW_CHAR)

    normal_col_adapt = col_adapt(normal_col, characters = COLUMNS_OLD_CHAR, new_char = COLUMNS_NEW_CHAR)
    tumor_col_adapt = col_adapt(tumor_col, characters = INDEX_OLD_CHAR, new_char = INDEX_NEW_CHAR)

    # imputation
    if os.path.exists(workdir) == False:
        os.makedirs(workdir)
    os.chdir(workdir)

    labeldir = workdir + '/' + data_label
    if os.path.exists(labeldir) == False:
        os.makedirs(labeldir)
    os.chdir(labeldir)

    tumor_col_adapt.to_csv('./tumor_col.txt', sep='\t', index=False)
    normal_col_adapt.to_csv('./normal_col.txt', sep='\t', index=False)
    data_fna_lt_zn_da.to_csv('./data.csv')

    # R script
    os.system('Rscript %s ./data.csv ./tumor_col.txt ./normal_col.txt ./data_imput.csv'%(r_bpca_path))

    data_imput = pd.read_csv('./data_imput.csv', index_col=0)

    data_tumor_imput = data_imput.loc[:, tumor_col_adapt['tumor']]
    data_normal_imput = data_imput.loc[:, normal_col_adapt['normal']]
    return data_imput, data_normal_imput, data_tumor_imput, normal_col_adapt, tumor_col_adapt


def create_nulldata(data_normal_imput, normal_col_adapt, ND_SAMPLE_SIZE):
    """
    利用normal data创建null data
    :param data_normal_imput: 预处理后的蛋白质丰度数据
    :param normal_col_adapt: normal队列列名
    :param ND_SAMPLE_SIZE: null data sample size
    :return nd_data_list: null data 
    :return nd_cohort_1_list: null data(cohort 1)
    :return nd_cohort_2_list: null data(cohort 2)
    """
    # 利用normal data 构建null data
    nd_data_list = []
    nd_cohort_1_list = []
    nd_cohort_2_list = []

    n_nd_normal = normal_col_adapt.shape[0] // (2 * ND_SAMPLE_SIZE)

    normal_index_list = random.sample(range(normal_col_adapt.shape[0]), n_nd_normal * 2 * ND_SAMPLE_SIZE)

    for idx in range(n_nd_normal):
        nd_cohort_1_list.append(data_normal_imput.iloc[:, normal_index_list[2 * idx * ND_SAMPLE_SIZE : (2 * idx + 1) * ND_SAMPLE_SIZE]])
        nd_cohort_2_list.append(data_normal_imput.iloc[:, normal_index_list[(2 * idx + 1) * ND_SAMPLE_SIZE : (2 * idx + 2) * ND_SAMPLE_SIZE]])
        nd_data_list.append(data_normal_imput.iloc[:, normal_index_list[2 * idx * ND_SAMPLE_SIZE : (2 * idx + 2) * ND_SAMPLE_SIZE]])
    
    return nd_data_list, nd_cohort_1_list, nd_cohort_2_list


def run_ttest(normal_run_path, tumor_test_path, workdir, Q_THRES = 0.05, D_THRES = 0, SAVE_OUT = False):
    """进行ttest 差异蛋白分析， p value 版本

    :param normal_run_path: normal file path.
    :param tumor_test_path: tumor file path.
    :param Q_THRES: 假说检验的阈值
    :param D_THRES: 差异蛋白的阈值
    :param SAVE_OUT: 是否保存输出结果
    :param save_dir: 输出目录
    :return: ttest_up, ttest_down(ttest 得到的调控蛋白)

    """
    os.chdir(workdir)
    normal_run = pd.read_csv(normal_run_path, index_col=0)
    tumor_test = pd.read_csv(tumor_test_path, index_col=0)

    P_THRES = 0.05

    ttest_list = []
    for ids in normal_run.index:
        _,p = stats.levene(normal_run.loc[ids,:], tumor_test.loc[ids, :])
        ttest_list.append(p)

    ttest_rate = (np.array(ttest_list) > P_THRES).sum() / len(ttest_list)

    _equal_val = []
    for tt in ttest_list:
        if tt > 0.05:
            _equal_val.append(False)
        else:
            _equal_val.append(True)

    direct_list = (tumor_test.mean(axis=1) - normal_run.mean(axis=1)).values
    # 计算差异
    p_value_list = []
    for idx in range(tumor_test.shape[0]):
        p_v = stats.ttest_ind(tumor_test.iloc[idx,:].values,normal_run.iloc[idx,:].values, equal_var = _equal_val[idx]).pvalue
        p_value_list.append(p_v)

    q_value_list = list(r_stats.p_adjust(FloatVector(p_value_list), method = 'BH'))

    tumor_ttest_up = copy.deepcopy(tumor_test)
    tumor_ttest_down = copy.deepcopy(tumor_test)
    tumor_ttest_qvalues = copy.deepcopy(tumor_test)
    for idx in range(tumor_test.shape[0]):
        tumor_ttest_qvalues.iloc[idx, :] = q_value_list[idx]

    fdr_down, fdr_up = [], []
    num_down, num_up = [], []

    tumor_ttest_up.iloc[:,:] = 0
    tumor_ttest_down.iloc[:,:] = 0

    tumor_ttest_up[(np.array(p_value_list) < Q_THRES) & (direct_list > D_THRES)] = 1
    tumor_ttest_down[(np.array(p_value_list) < Q_THRES) & (direct_list < -D_THRES)] = 1

    if SAVE_OUT:

        tumor_ttest_up.to_csv('./ttest_up.csv')
        tumor_ttest_down.to_csv('./ttest_down.csv')
        tumor_ttest_qvalues.to_csv('./ttest_qvalues.csv')

    return tumor_ttest_up, tumor_ttest_down


def run_wilcox(normal_run_path, tumor_test_path, workdir, Q_THRES = 0.05, D_THRES = 0, SAVE_OUT = False):
    """进行wilcox 差异蛋白分析, p value版本

    :param normal_run_path: normal file path.
    :param tumor_test_path: tumor file path.
    :param THRES: 假说检验的阈值
    :param SAVE_OUT: 是否保存输出结果
    :param save_dir: 输出目录
    :return: wilcox_up, wilcox_down(wilcoxon 得到的调控蛋白)

    """
    os.chdir(workdir)
    normal_run = pd.read_csv(normal_run_path, index_col=0)
    tumor_test = pd.read_csv(tumor_test_path, index_col=0)

    P_THRES = 0.05
    shapiro_list = []
    for ids in normal_run.index:
        _,p = stats.shapiro(normal_run.loc[ids,:])
        shapiro_list.append(p)

        _,p = stats.shapiro(tumor_test.loc[ids,:])
        shapiro_list.append(p)

    shapiro_rate = (np.array(shapiro_list) > P_THRES).sum() / len(shapiro_list)

    direct_list = (tumor_test.mean(axis=1) - normal_run.mean(axis=1)).values
    p_value_list = []
    for idx in tumor_test.index:
        p_v = stats.mannwhitneyu(tumor_test.loc[idx,:].values,normal_run.loc[idx,:].values).pvalue
        p_value_list.append(p_v)

    q_value_list = list(r_stats.p_adjust(FloatVector(p_value_list), method = 'BH'))
    tumor_wilcox_up = copy.deepcopy(tumor_test)
    tumor_wilcox_down = copy.deepcopy(tumor_test)
    tumor_wilcox_qvalues = copy.deepcopy(tumor_test)
    for idx in range(tumor_test.shape[0]):
        tumor_wilcox_qvalues.iloc[idx, :] = q_value_list[idx]

    fdr_down, fdr_up = [], []
    num_down, num_up = [], []

    tumor_wilcox_up.iloc[:,:] = 0
    tumor_wilcox_down.iloc[:,:] = 0

    tumor_wilcox_up[(np.array(p_value_list) < Q_THRES) & (direct_list > D_THRES)] = 1
    tumor_wilcox_down[(np.array(p_value_list) < Q_THRES) & (direct_list < -D_THRES)] = 1

    if SAVE_OUT:
        tumor_wilcox_up.to_csv('./wilcox_up.csv')
        tumor_wilcox_down.to_csv('./wilcox_down.csv')
        tumor_wilcox_qvalues.to_csv('./wilcox_qvalues.csv')

    return tumor_wilcox_up, tumor_wilcox_down


def rankc_j2(reoa_path, normal_data, tumor_data, workdir, cycle_rankc, fdr, max_exception_rankc):
    """
    Rankcomp找差异蛋白
    :param normal_data: normal data
    :param tumor_data: tumor data
    :param workdir: rankcomp j1 work sapace
    :param cycle_rankc: max cycles for filtering dysregulated genes
    :param fdr: FDR (False Discovery Rate) level
    :param max_exception_rankc: list of FDR level or exception number, which are used to select significantly stable pairs
    """
    n_n_sample = normal_data.shape[1]
    n_t_sample = tumor_data.shape[1]
    n_prot = normal_data.shape[0]
    # run rankcomp j1
    os.system('%s -s 1 -j 2 -m %d -f %f -a 2 2 %d %s/cohort_1.dat %s/cohort_2.dat %d %d %f %f'%
              (reoa_path, cycle_rankc, fdr, n_prot, workdir, workdir, n_n_sample, n_t_sample, 
               max_exception_rankc, max_exception_rankc))

def rankc_j2_v2(reoa_path, normal_data, tumor_data, workdir, cycle_rankc, fdr, max_exception_rankc):
    """
    Rankcomp找差异蛋白
    :param normal_data: normal data
    :param tumor_data: tumor data
    :param workdir: rankcomp j1 work sapace
    :param cycle_rankc: max cycles for filtering dysregulated genes
    :param fdr: FDR (False Discovery Rate) level
    :param max_exception_rankc: list of FDR level or exception number, which are used to select significantly stable pairs
    """
    n_n_sample = normal_data.shape[1]
    n_t_sample = tumor_data.shape[1]
    n_prot = normal_data.shape[0]
    # run rankcomp j1
    os.system('%s -s 1 -j 2 -m %d -f %f -a 0 2 %d %s/cohort_1.dat %s/cohort_2.dat %d %d %f %f'%
              (reoa_path, cycle_rankc, fdr, n_prot, workdir, workdir, n_n_sample, n_t_sample, 
               max_exception_rankc, max_exception_rankc))
    
def get_rankc_j2_results(workdir):
    """
    读取Rankcomp j2结果中的差异蛋白
    :param workdir: rankcomp j2 work space
    :return sp_concordant: concordant
    :return sp_reversed: reversed
    :return sp_normal: normal stable pairs
    """
    try:
        dep = pd.read_csv(workdir + '/gene_state_1.dat',sep='\t', header=None)    
    except:
        print('Error: DEP is empty!\n')
        
    return dep


def run_penda_fdr(rscript_path, workdir, fdr):
    os.chdir(workdir)
    os.system('Rscript %s ./cohort_1.csv ./cohort_2.csv %s ./penda_fdr_down.csv ./penda_fdr_up.csv'%(rscript_path, fdr))

    _penda_fdr_up = pd.read_csv('./penda_fdr_up.csv', index_col=0)
    _penda_fdr_down = pd.read_csv('./penda_fdr_down.csv', index_col=0)
    _penda_fdr_result = _penda_fdr_up | _penda_fdr_down

    return _penda_fdr_up, _penda_fdr_down, _penda_fdr_result

def run_penda(rscript_path, workdir):
    os.chdir(workdir)
    os.system('Rscript %s ./cohort_1.csv ./cohort_2.csv ./penda_down.csv ./penda_up.csv'%(rscript_path))

    _penda_up = pd.read_csv('./penda_up.csv', index_col=0)
    _penda_down = pd.read_csv('./penda_down.csv', index_col=0)
    _penda_result = _penda_up | _penda_down
    return _penda_up, _penda_down, _penda_result


def get_nd_param_sn(nd_result, nd_set, method, GROUP_DEP_THRES=0.5):
    sn = []
#     nd_data_label = []
    nd_param_type = []
    nd_method = []
    for k in range(len(nd_result)):
        _result = nd_result[k]
        _dep_prot = _result.index[(_result != 0).sum(axis=1) >= np.floor(_result.shape[1] * GROUP_DEP_THRES)].tolist()
        _undep_prot = _result.index[(_result != 0).sum(axis=1) < np.floor(_result.shape[1] * GROUP_DEP_THRES)].tolist()
        
        _dep_param = nd_set[k].loc[_dep_prot, :].mean(axis=1).values
        _undep_param = nd_set[k].loc[_undep_prot, :].mean(axis=1).values
        _sn = (np.mean(_dep_param) - np.mean(_undep_param)) / (np.std(_dep_param) + np.std(_undep_param))
        sn.append(_sn)
        nd_param_type.append('mean')
        nd_method.append(method)
        
        _dep_param = nd_set[k].loc[_dep_prot, :].var(axis=1).values
        _undep_param = nd_set[k].loc[_undep_prot, :].var(axis=1).values
        _sn = (np.mean(_dep_param) - np.mean(_undep_param)) / (np.std(_dep_param) + np.std(_undep_param))
        sn.append(_sn)
        nd_param_type.append('var')
        nd_method.append(method)
        
        _dep_param = (nd_set[k].loc[_dep_prot, :].std(axis=1) / nd_set[k].loc[_dep_prot, :].mean(axis=1)).values
        _undep_param = (nd_set[k].loc[_undep_prot, :].std(axis=1) / nd_set[k].loc[_undep_prot, :].mean(axis=1)).values
        _sn = (np.mean(_dep_param) - np.mean(_undep_param)) / (np.std(_dep_param) + np.std(_undep_param))
        sn.append(_sn)
        nd_param_type.append('cv')
        nd_method.append(method)
        
    
    nd_param_sn = pd.DataFrame({'sn': sn, 'param_type': nd_param_type, 'method': nd_method})
    
    return nd_param_sn





