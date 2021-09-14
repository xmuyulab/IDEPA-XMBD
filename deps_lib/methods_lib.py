import shutil
import random
import sys
import datetime
import copy 
import os

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from scipy import stats

from deps_lib import utils
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector

r_stats = importr('stats')


def quantile_dep(normal_data_path, tumor_data_path, workdir, quantile, factor, SAVE_OUT):
    
    os.chdir(workdir)
    normal_data = pd.read_csv(normal_data_path, index_col=0)
    tumor_data = pd.read_csv(tumor_data_path, index_col=0)

    normal_data = normal_data.loc[tumor_data.index.tolist(),:]

    quantile_up = copy.deepcopy(tumor_data)
    quantile_down = copy.deepcopy(tumor_data)
    quantile_up.iloc[:, :] = 0
    quantile_down.iloc[:, :] = 0
    for idx in normal_data.index:
        quantile_high = np.quantile(normal_data.loc[idx, :].values, 1-quantile)
        quantile_low = np.quantile(normal_data.loc[idx, :].values, quantile)
        if quantile_high > 0:
            up_thres = quantile_high * factor
        else:
            up_thres = quantile_high / factor

        if quantile_low > 0:
            down_thres = quantile_low / factor
        else:
            down_thres = quantile_low * factor
        quantile_up.loc[idx, :][tumor_data.loc[idx, :] > up_thres] = 1
        quantile_down.loc[idx, :][tumor_data.loc[idx, :] < down_thres] = 1
    if SAVE_OUT:
        quantile_up.to_csv('./quantile_up.csv')
        quantile_down.to_csv('./quantile_down.csv')
    
    return quantile_up, quantile_down 


def run_ttest(normal_run_path, tumor_test_path, workdir, Q_THRES = 0.05, D_THRES = 0, SAVE_OUT = False):
    """进行ttest 差异蛋白分析

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

    tumor_ttest_up[(np.array(q_value_list) < Q_THRES) & (direct_list > D_THRES)] = 1
    tumor_ttest_down[(np.array(q_value_list) < Q_THRES) & (direct_list < -D_THRES)] = 1

    if SAVE_OUT:

        tumor_ttest_up.to_csv('./ttest_up.csv')
        tumor_ttest_down.to_csv('./ttest_down.csv')
        tumor_ttest_qvalues.to_csv('./ttest_qvalues.csv')

    return tumor_ttest_up, tumor_ttest_down


def run_wilcox(normal_run_path, tumor_test_path, workdir, Q_THRES = 0.05, D_THRES = 0, SAVE_OUT = False):
    """进行wilcox 差异蛋白分析

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

    tumor_wilcox_up[(np.array(q_value_list) < Q_THRES) & (direct_list > D_THRES)] = 1
    tumor_wilcox_down[(np.array(q_value_list) < Q_THRES) & (direct_list < -D_THRES)] = 1

    if SAVE_OUT:
        tumor_wilcox_up.to_csv('./wilcox_up.csv')
        tumor_wilcox_down.to_csv('./wilcox_down.csv')
        tumor_wilcox_qvalues.to_csv('./wilcox_qvalues.csv')

    return tumor_wilcox_up, tumor_wilcox_down



def run_peng_method(normal_data, tumor_data, thres = 0.95):
    """Peng methods for DE protein
    
    @param normal_data
    @param tumor_data
    @return dep_up, dep_down
    """
    
    normal_rank = expression_2_rank(normal_data)
    normal_rank.columns = normal_data.columns

    tumor_rank = expression_2_rank(tumor_data)
    tumor_rank.columns = tumor_data.columns

    #### 1.确定normal中的稳定对
    sp_normal = get_stable_pairs(data = normal_rank,
                                thres = thres)
    # sp_tumor = get_stable_pairs(data = tumor_rank,
    #                             thres = thres)

    #### 确定逆转对
    # Fish exac test
    p_v = []
    for idx in sp_normal.index:
#         if (idx % 10000) == 0:
#             t = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
#             print("%s/%s: %s "%(idx, sp_normal.shape[0], t))
        idx_1 = sp_normal.loc[idx, 'idx_1']
        idx_2 = sp_normal.loc[idx, 'idx_2']
        _a = ((normal_rank.iloc[idx_1, :] - normal_rank.iloc[idx_2, :]) > 0).sum()
        _b = ((normal_rank.iloc[idx_1, :] - normal_rank.iloc[idx_2, :]) < 0).sum()
        _c = ((tumor_rank.iloc[idx_1, :] - tumor_rank.iloc[idx_2, :]) > 0).sum()
        _d = ((tumor_rank.iloc[idx_1, :] - tumor_rank.iloc[idx_2, :]) < 0).sum()
        p_v.append(stats.fisher_exact([[_a,_b],[_c,_d]])[1])

    q_v = list(r_stats.p_adjust(FloatVector(p_v), method = 'BH'))
    rev = sp_normal[np.array(q_v) < 0.1].reset_index(drop=True)
    p_direction = ((normal_rank.mean(axis=1) - tumor_rank.mean(axis=1)) > 0).values

    #### 3.过滤逆转对
    # 1. 变化方向与目标蛋白一致
    # 2. 变异系数低（Top 3）
    rev['idx_1_direction'] = p_direction[rev['idx_1'].values]
    rev['idx_2_direction'] = p_direction[rev['idx_2'].values]
    rev['concordant'] = (rev['idx_1_direction']&rev['idx_2_direction']) | (~rev['idx_1_direction'])&(~rev['idx_2_direction'])

    # 1. 变化方向一致
    rev = rev[rev['concordant']].reset_index(drop=True)

    # 2. 变异系数
    data_rank = pd.concat([normal_rank, tumor_rank], axis=1)
    cv_list = []
    for idx in data_rank.index:
        cv_list.append(cov(data_rank.loc[idx, :]))

    #### 4.确定参考组
    refer_p_set = {}
    for idx in range(normal_rank.shape[0]):
        n_1 = rev[rev['idx_1'] == idx].shape[0]
        n_2 = rev[rev['idx_2'] == idx].shape[0]
        refer_p = []
        if (n_2 + n_1) == 0:
            pass
        else:
            if (n_1 != 0):
                refer_p.extend(rev[rev['idx_1'] == idx]['idx_2'].values.tolist())
            if (n_2 != 0):
                refer_p.extend(rev[rev['idx_2'] == idx]['idx_1'].values.tolist())

        if len(refer_p) == 0:
            pass
        elif len(refer_p) <= 3:
            refer_p_set[idx] = refer_p
        elif len(refer_p) > 3:
            refer_p_f = []
            for i in np.partition(np.array(cv_list)[refer_p], 1)[:3]:
                refer_p_f.append(refer_p[np.where(np.array(cv_list)[refer_p] == i)[0][0]])
            refer_p_set[idx] = refer_p_f

    #### 5.确定差异蛋白
    dep_up = copy.deepcopy(tumor_rank)
    dep_down = copy.deepcopy(tumor_rank)

    dep_up.iloc[:, :] = False
    dep_down.iloc[:, :] = False

    for k in refer_p_set.keys():

        r_p = refer_p_set[k]
        p_d = p_direction[k]

        pd_tmp = np.zeros([len(r_p),tumor_rank.shape[1]])

        for _i, _p in enumerate(zip(r_p)):
            _p = _p[0]
            _l = ((normal_rank.iloc[k, :] - normal_rank.iloc[_p, :]) > 0).sum() > normal_rank.shape[1] * 0.8
            if _l:
                pd_tmp[_i, :] = ((tumor_rank.iloc[k, :] - tumor_rank.iloc[_p, :]) < 0).values
            else:
                pd_tmp[_i, :] = ((tumor_rank.iloc[k, :] - tumor_rank.iloc[_p, :]) > 0).values

        _dep = pd_tmp.sum(axis=0) > len(r_p) * 0.5

        if p_d:
            dep_down.iloc[k,:] = _dep
        else:
            dep_up.iloc[k, :] = _dep
            
    return dep_up, dep_down



# coefficient of variation
def cov(se):
    return np.std(se) / np.mean(se)

def expression_2_rank(expression_data):
    """
    将样本的丰度数据转换为次序数据
    :param expression_data: protein abundance
    :return rank_data: rank data
    """
    index_list = expression_data.index
    rank_data = pd.DataFrame(columns=['case'])
    
    for col in expression_data.columns:
        temp = pd.DataFrame(expression_data[col].sort_values())
        temp['%s_rank'%(col)] = np.arange(1, expression_data.shape[0]+1)
        temp.sort_index(inplace = True)
        
        rank_data['%s_rank'%(col)] = temp['%s_rank'%(col)]
    rank_data.drop(['case'], axis=1, inplace=True)
    rank_data = rank_data.loc[index_list,:]
    return rank_data


def get_stable_pairs(data, thres):
    idx_1 = [i for i in np.arange(data.shape[0]).tolist() for j in np.arange(data.shape[0])]
    idx_2 = np.arange(data.shape[0]).tolist() * data.shape[0]

    pairs = pd.DataFrame({'idx_1': idx_1, 'idx_2': idx_2})

    pairs = pairs[(pairs['idx_2'] - pairs['idx_1']) > 0].reset_index(drop=True)

    idx_data_set = {}
    for _idx in range(data.shape[0]):
        idx_data_set[_idx] = data.iloc[_idx,:].values

    tmp = []
    for _index in pairs.index:
        tmp.append(idx_data_set[pairs.loc[_index, 'idx_1']])
    pairs['values_1'] = tmp

    tmp = []
    for _index in pairs.index:
        tmp.append(idx_data_set[pairs.loc[_index, 'idx_2']])
    pairs['values_2'] = tmp

    pairs['diff'] = (pairs['values_1'] - pairs['values_2']) 

    label_list = []
    diff_label = []
    for _index in pairs.index:
        l_1 = ((pairs.loc[_index, 'diff'] > 0).sum() / data.shape[1]) > thres
        l_2 = ((pairs.loc[_index, 'diff'] < 0).sum() / data.shape[1]) > thres
        label_list.append(l_1 | l_2)
        
        _m_1 = data.iloc[pairs.loc[_index, 'idx_1'], :].mean()
        _m_2 = data.iloc[pairs.loc[_index, 'idx_2'], :].mean()
        if  _m_1 > _m_2:
            diff_label.append(True)
        else:
            diff_label.append(False)
        
    pairs['1>2'] = diff_label
    pairs['pairs_id'] = pairs['idx_1'].astype('str') + '_' + pairs['idx_2'].astype('str')
    returned_pairs = pairs[label_list].reset_index(drop=True).loc[:, ['idx_1', 'idx_2', '1>2', 'pairs_id']]
    return returned_pairs