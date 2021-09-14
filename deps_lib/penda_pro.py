import os
import copy
import datetime

import pandas as pd
import numpy as np
import scipy.stats as stats

from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector

r_stats = importr('stats')

def run_penda_pro(normal_path, tumor_path, FDR = 0.05, CONVERGENCE_THRESHOLD = 0.99, MAX_CYCLE = 48, K = 20, THRESHOLD_LH = 0.99):
    
    normal = pd.read_csv(normal_path, index_col=0)
    tumor = pd.read_csv(tumor_path, index_col=0)
    
    # 对normal队列进行从大到小排序
    normal['mean'] = normal.mean(axis=1)
    normal_sort = normal.sort_values(by='mean', ascending=False)
    normal_sort = normal_sort.drop(['mean'], axis=1)
    tumor_sort = tumor.loc[normal_sort.index.tolist(), :]
    tumor_rank = copy.deepcopy(tumor_sort)
    normal_rank = copy.deepcopy(normal_sort)

    for col in normal_sort.columns:
        _normal_sort = copy.deepcopy(normal_sort)
        _normal_sort = _normal_sort.sort_values(by=col, axis=0, ascending=False)
        _normal_sort['tmp'] = np.arange(1, _normal_sort.shape[0] + 1)
        normal_rank[col] = _normal_sort['tmp']

    for col in tumor_sort.columns:
        _tumor_sort = copy.deepcopy(tumor_sort)
        _tumor_sort = _tumor_sort.sort_values(by=col, axis=0, ascending=False)
        _tumor_sort['tmp'] = np.arange(1, _tumor_sort.shape[0] + 1)
        tumor_rank[col] = _tumor_sort['tmp']


    H = {}
    L = {}

    for idx in normal_rank.index:
        _normal_label_L = copy.deepcopy(normal_rank)
        _normal_label_H = copy.deepcopy(normal_rank)
        for col in normal_rank.columns:
            _normal_label_H[col] = normal_rank.loc[:,col] < normal_rank.loc[idx, col]
            _normal_label_L[col] = normal_rank.loc[:,col] > normal_rank.loc[idx, col]

        _H = _normal_label_H[_normal_label_H.sum(axis=1) > _normal_label_H.shape[1] * THRESHOLD_LH].index.tolist()
        _L = _normal_label_L[_normal_label_L.sum(axis=1) > _normal_label_L.shape[1] * THRESHOLD_LH].index.tolist()
        H[idx] = _H
        L[idx] = _L

        up_qvalues = copy.deepcopy(tumor_rank)
        down_qvalues = copy.deepcopy(tumor_rank)
    
    _tt = 0
    for col in tumor_rank.columns:
        _tt += 1
        t = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print('Columns : %s (%s / %s) %s'%(col, _tt, tumor_rank.shape[1], t))

        R = {}
        for idx in tumor_rank.index:
        # idx = tumor_col_rank.index[100]
            r_H = H[idx]
            r_L = L[idx]

            if len(r_H) >= K and len(r_L) >= K:
                R[idx] = r_L[:K] + r_H[-K:]

            elif len(r_H) < K and len(r_L) >= K:
                _n = K - len(r_H)
                if len(r_L) >= K + _n:
                    pass
                else:
                    r_L.extend('0' * (K + _n - len(r_L)))
                R[idx] = r_L[:(K+_n)] + r_H[-len(_H):]

            elif len(r_H) >= K and len(r_L) < K:
                _n = K - len(r_L)
                if len(r_H) >= K + _n:
                    pass
                else:
                    r_H.extend('0' * (K + _n - len(r_H)))
                R[idx] = r_L[:len(r_L)] + r_H[-(K+_n):]

            elif len(r_H) < K and len(r_L) < K:
                r_H.extend('0' * (K - len(r_H)))
                r_L.extend('0' * (K - len(r_L)))
                R[idx] = r_H[:K] + r_L[-K:]

        R_df = pd.DataFrame(R).T

        # 第一次计算P_value
        p_values = []
        p_values_less = []
        p_values_greater = []

        for idx in tumor_rank.index:
            r = R_df.loc[idx, :]
            
            if '0' in r.values:
                r = np.setdiff1d(r.values,'0')
            
            _normal_tmp = normal_rank.loc[r,:].mean(axis=1) - normal_rank.loc[idx, :].mean()
            num_normal_H = (_normal_tmp <= 0).sum()
            num_normal_L = (_normal_tmp > 0).sum()

            _tumor_tmp = tumor_rank.loc[r,col] - tumor_rank.loc[idx, col]
            num_tumor_H = (_tumor_tmp <= 0).sum()
            num_tumor_L = (_tumor_tmp > 0).sum()

            confusion_matrix = [[num_normal_H, num_normal_L], [num_tumor_H, num_tumor_L]]
            _, _p = stats.fisher_exact(confusion_matrix)
            p_values.append(_p)
            _, _p_less = stats.fisher_exact(confusion_matrix, alternative='less')
            p_values_less.append(_p_less)
            _, _p_greater = stats.fisher_exact(confusion_matrix, alternative='greater')
            p_values_greater.append(_p_greater)

            q_values = list(r_stats.p_adjust(FloatVector(p_values), method = 'BH'))
            q_values_less = list(r_stats.p_adjust(FloatVector(p_values_less), method = 'BH'))
            q_values_greater = list(r_stats.p_adjust(FloatVector(p_values_greater), method = 'BH'))

        p_values_df = pd.DataFrame({'two_side':p_values, 'down':p_values_less, 'up':p_values_greater}, index=tumor_rank.index)
        q_values_df = pd.DataFrame({'two_side':q_values, 'down':q_values_less, 'up':q_values_greater}, index=tumor_rank.index)

        # 迭代
        q_values_initial = q_values_df

        down_initial = q_values_df[q_values_initial['down'] < FDR].index.tolist()
        up_initial = q_values_df[q_values_initial['up'] < FDR].index.tolist()

        for cycle_id in range(MAX_CYCLE):
            t = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
#             print('# Cycle (%s/%s): %s'%(cycle_id, MAX_CYCLE, t))

            R_new = {}
            for idx in tumor_rank.index:

                _H = H[idx]
                _L = L[idx]

                _H_down = [x for x in _H if x not in down_initial]
                _L_up = [x for x in _L if x not in up_initial]

                if len(_H_down) >= K and len(_L_up) >= K:
                    R_new[idx] = _L_up[:K] + _H_down[-K:]

                elif len(_H_down) < K and len(_L_up) >= K:
                    _n = K - len(_H_down)
                    if len(_L_up) >= K + _n:
                        pass
                    else:
                        _L_up.extend('0' * (K + _n - len(_L_up)))
                    R_new[idx] = _L_up[:(K+_n)] + _H_down[-len(_H_down):]

                elif len(_H_down) >= K and len(_L_up) < K:
                    _n = K - len(_L_up)
                    if len(_H_down) >= K + _n:
                        pass
                    else:
                        _H_down.extend('0' * (K + _n - len(_H_down)))
                    R_new[idx] = _L_up[:len(_L_up)] + _H_down[-(K+_n):]

                elif len(_H_down) < K and len(_L_up) < K:
                    _H_down.extend('0' * (K - len(_H_down)))
                    _L_up.extend('0' * (K - len(_L_up)))
                    R_new[idx] = _H_down[:K] + _L_up[-K:]

            R_new_df = pd.DataFrame(R_new).T

            # 计算FDR值
            p_values_new = []
            p_values_less_new = []
            p_values_greater_new = []

            for idx in tumor_rank.index:
                r = R_new_df.loc[idx, :]
                if '0' in r.values:
                    r = np.setdiff1d(r.values,'0')
                _normal_tmp = normal_rank.loc[r,:].mean(axis=1) - normal_rank.loc[idx, :].mean()
                num_normal_H = (_normal_tmp <= 0).sum()
                num_normal_L = (_normal_tmp > 0).sum()

                _tumor_tmp = tumor_rank.loc[r,col] - tumor_rank.loc[idx, col]
                num_tumor_H = (_tumor_tmp <= 0).sum()
                num_tumor_L = (_tumor_tmp > 0).sum()

                confusion_matrix = [[num_normal_H, num_normal_L], [num_tumor_H, num_tumor_L]]
                _, _p = stats.fisher_exact(confusion_matrix)
                p_values_new.append(_p)
                _, _p_less = stats.fisher_exact(confusion_matrix, alternative='less')
                p_values_less_new.append(_p_less)
                _, _p_greater = stats.fisher_exact(confusion_matrix, alternative='greater')
                p_values_greater_new.append(_p_greater)

                q_values_new = list(r_stats.p_adjust(FloatVector(p_values_new), method = 'BH'))
                q_values_less_new = list(r_stats.p_adjust(FloatVector(p_values_less_new), method = 'BH'))
                q_values_greater_new = list(r_stats.p_adjust(FloatVector(p_values_greater_new), method = 'BH'))

            # p_values_df_new = pd.DataFrame({'two_side':p_values_new, 'down':p_values_less_new, 'up':p_values_greater_new}, index=tumor_rank.index)
            q_values_df_new = pd.DataFrame({'two_side':q_values_new, 'down':q_values_less_new, 'up':q_values_greater_new}, index=tumor_rank.index)

            down_new = q_values_df_new[q_values_df_new['down'] < FDR].index.tolist()
            up_new = q_values_df_new[q_values_df_new['up'] < FDR].index.tolist()

            down_intersection = set(down_new).intersection(set(down_initial))
            up_intersection = set(up_new).intersection(set(up_initial))

            if len(down_new) != 0 and len(up_new) != 0:
                down_inter_rate = len(down_intersection) / len(down_new)
                up_inter_rate = len(up_intersection) / len(up_new)
            else:
                down_inter_rate = 0
                up_inter_rate = 0
#             print('-- down_inter_rate: %.4f'%(down_inter_rate))
#             print('-- down_num: %s'%(len(down_new)))
#             print('-- up_inter_rate: %.4f'%(up_inter_rate))
#             print('-- up_num: %s'%(len(up_new)))
            if (down_inter_rate >= CONVERGENCE_THRESHOLD) and (up_inter_rate >= CONVERGENCE_THRESHOLD):
                break
            else:
                q_values_initial = q_values_df_new
                down_initial = down_new
                up_initial = up_new
                

        up_qvalues.loc[:, col] = q_values_df_new.loc[:, 'up']
        down_qvalues.loc[:, col] = q_values_df_new.loc[:, 'down']

    return up_qvalues, down_qvalues