import pandas as pd
import numpy as np
import random
import seaborn as sns
import matplotlib.pyplot as plt
import copy
import os



def get_simu_data_zscore(orig_data, DEREGULATED_RATE):
    """
    对经过z-score标准化数据进行处理，生成仿真数据
    :param orig_data
    :param DEREGULATED_RATE: 表达量数据改变比例
    :return simulation data
    """
    simu_data = copy.deepcopy(orig_data)
    for col in orig_data.columns:
        ind_list = random.sample(orig_data.index.tolist(), int(orig_data.shape[0] * DEREGULATED_RATE))
        for ind in ind_list:

            _direction = random.sample([1,-1], 1)[0]
            _range = random.sample([2,4,6], 1)[0]

            if _range == 2:
                _value = np.random.normal(loc =0.2 , scale= 0.05,size = (1))[0] 
            elif _range == 4:
                _value = np.random.normal(loc =0.4 , scale= 0.05,size = (1))[0] 
            else:
                _value = np.random.normal(loc =0.6 , scale= 0.05,size = (1))[0] 

            if _direction == 1:
                simu_data.loc[ind, col] = orig_data.loc[ind, col] + _value
            else:
                simu_data.loc[ind, col] = orig_data.loc[ind, col] - _value
    return simu_data



def create_simu_data(input_data_path, input_normal_path, output_workdir, DEP_RATE):

    if os.path.exists(output_workdir) == False:
        os.makedirs(output_workdir)
    os.chdir(output_workdir)

    data = pd.read_csv(input_data_path, index_col=0)
    normal_cohort = pd.read_csv(input_normal_path)

    tumor_cohort = copy.deepcopy(normal_cohort)
    tumor_cohort.columns = ['tumor']

    for i in range(tumor_cohort.shape[0]):
        tumor_cohort.iloc[i, 0] = normal_cohort.iloc[i, 0] + '_T'

    normal_data = data.loc[:, normal_cohort['normal']]
    simu_tumor_data = get_simu_data_zscore(orig_data = normal_data, DEREGULATED_RATE = DEP_RATE)
    simu_tumor_data.columns = tumor_cohort['tumor']

    simu_data = pd.concat([normal_data, simu_tumor_data], axis=1)

    simu_data.to_csv('./simu_data.csv')
    normal_cohort.to_csv('./normal.txt', sep='\t', index=False)
    tumor_cohort.to_csv('./tumor.txt', sep='\t', index=False)

    mc_data_workdir = output_workdir + '/methodsComp'
    if os.path.exists(mc_data_workdir) == False:
        os.makedirs(mc_data_workdir)
    os.chdir(mc_data_workdir)

    # 创建训练集与测试集
    cohort = pd.concat([normal_cohort, tumor_cohort], axis=1)
    paired_index = random.sample(cohort.index.tolist(), int(cohort.shape[0] * 0.3))
    paired_samples = cohort.iloc[paired_index,:].reset_index(drop=True)
    paired_data = simu_data.loc[:, paired_samples['normal'].tolist() + paired_samples['tumor'].tolist()]

    run_index = []
    for i in range(cohort.shape[0]):
        if i in paired_index:
            pass
        else:
            run_index.append(i)

    run_samples = cohort.iloc[run_index,:].reset_index(drop=True)
    mc_data = simu_data.loc[:, run_samples['normal'].tolist() + paired_samples['tumor'].tolist()]
    mc_normal = pd.DataFrame({'normal': run_samples['normal']})
    mc_tumor = pd.DataFrame({'tumor': paired_samples['tumor']})

    mc_data.to_csv('./data.csv')
    paired_data.to_csv('./paired_data.csv')
    paired_samples.to_csv('./paired_samples.txt', sep='\t', index=False)
    mc_normal.to_csv('./normal.txt', sep='\t', index=False)
    mc_tumor.to_csv('./tumor.txt', sep='\t', index=False)