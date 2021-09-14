#!/usr/bin/env python3
# -*- coding: utf-8 -*-


'''
Stable Pairs Analysis Dependent Library
'''

__author__ = 'Liu Yachen'


import os
import copy
import shutil
import random

import numpy as np
import pandas as pd

import seaborn as sns
import matplotlib.pyplot as plt

from deps_lib import utils


plt.rcParams['font.sans-serif'] = ['SimHei']  # 中文字体设置-黑体
plt.rcParams['axes.unicode_minus'] = False  # 解决保存图像是负号'-'显示为方块的问题
sns.set(font_scale=1.5)  # 解决Seaborn中文显示问题并调整字体大小

class StablePairs(object):
    
    def __init__(self, rd, sp_workdir, SP_THRES, RANDOM_VISUAL, NUM_VISUAL, CYCLE_RANKC, FDR, MAX_EXCEPTION_RANKC, fig_save_path='', reoa_path=None):
        
        self.rd = rd
        self.sp_workdir = sp_workdir
        self.SP_THRES = SP_THRES
        self.RANDOM_VISUAL = RANDOM_VISUAL
        self.NUM_VISUAL = NUM_VISUAL
        self.CYCLE_RANKC = CYCLE_RANKC
        self.FDR = FDR
        self.MAX_EXCEPTION_RANKC = MAX_EXCEPTION_RANKC
        self.fig_save_path = fig_save_path
        self.reoa_path = reoa_path
        
    def expression_2_rank(self, expression_data):
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

    def plot_pairs_rank(self, reversed_pairs, normal_data, tumor_data, pair_id, save_dir, if_reversed=True):
        """
        可视化一致对或逆转对的丰度
        :param reversed_pairs: reversed pairs
        :param normal_data: abundance of normal cohort
        :param tumor_data: abundance of tumor cohort
        :param pair_id: number of stable pairs that need to be visualized
        :param save_dir: figure path
        :param if_reversed: True, reversed pairs; Flase, concordant pairs
        """
        normal_data_rank = self.expression_2_rank(normal_data)
        tumor_data_rank = self.expression_2_rank(tumor_data)

        self.plot_pairs_expression(reversed_pairs = reversed_pairs, normal_data = normal_data_rank, 
                              tumor_data = tumor_data_rank, pair_id = pair_id, save_dir = save_dir, if_reversed = if_reversed)


    def plot_pairs_expression(self, reversed_pairs, normal_data, tumor_data, pair_id, save_dir, if_reversed=True):
        """
        可视化一致对或逆转对的丰度
        :param reversed_pairs: reversed pairs
        :param normal_data: abundance of normal cohort
        :param tumor_data: abundance of tumor cohort
        :param pair_id: number of stable pairs that need to be visualized
        :param save_dir: figure path
        :param if_reversed: True, reversed pairs; Flase, concordant pairs
        """

        protein_id = list(normal_data.index)
        id_1, id_2 = reversed_pairs.iloc[pair_id, :][0], reversed_pairs.iloc[pair_id, :][1]

        normal_data_1 = normal_data.iloc[id_1, :].values
        normal_data_2 = normal_data.iloc[id_2, :].values
        tumor_data_1 = tumor_data.iloc[id_1, :].values
        tumor_data_2 = tumor_data.iloc[id_2, :].values

        normal_x_list = np.arange(1, normal_data_1.shape[0]+1)
        tumor_x_list = np.arange(1, tumor_data_1.shape[0]+1)

        plt.figure(figsize=(18, 8))

        if if_reversed:
            plt.suptitle('The expression of protein in normal and tumor: Reversed', fontsize=25, color='black', fontweight='light')
        else:
            plt.suptitle('The expression of protein in normal and tumor: Concordant', fontsize=25, color='black', fontweight='light')
        plt.subplot(1,2,1)
        sns.lineplot(x=normal_x_list, y=(normal_data_1), color='red', style=True, dashes=[(1,1)], markers=True)
        sns.lineplot(x=normal_x_list, y=(normal_data_2), color='black', style=True, dashes=[(1,1)],markers=True)
        plt.legend(['(normal)protein: %s'%(protein_id[id_1]),'(normal)protein: %s'%(protein_id[id_2])])
        plt.xticks(fontsize = 12)
        plt.yticks(fontsize = 12)
        plt.xlabel('sample number', fontsize=15)

        plt.ylabel('protein rank', fontsize=15)
        plt.title('normal', fontsize=20)

        plt.subplot(1,2,2)
        sns.lineplot(x=tumor_x_list, y=(tumor_data_1), color='red', style=True, dashes=[(1,1)], markers=True)
        sns.lineplot(x=tumor_x_list, y=(tumor_data_2), color='black', style=True, dashes=[(1,1)],markers=True)
        plt.legend(['(tumor)protein: %s'%(protein_id[id_1]),'(tumor)protein: %s'%(protein_id[id_2])])
        plt.xticks(fontsize = 12)
        plt.yticks(fontsize = 12)
        plt.xlabel('sample number', fontsize=15)

        plt.ylabel('protein rank', fontsize=15)
        plt.title('tumor', fontsize=20)

        if if_reversed:
            plt.savefig(save_dir + '/%s-%s-reversed.pdf'%(protein_id[id_1], protein_id[id_2]), dpi=800, format='pdf', bbox_inches='tight')
        else:
            plt.savefig(save_dir + '/%s-%s-concordant.pdf'%(protein_id[id_1], protein_id[id_2]), dpi=800, format='pdf', bbox_inches='tight')


    def select_stable_pairs(self, specific_prot, express_data, stable_pairs_result):
        """
        包含特异性蛋白的稳定对
        :param specific_prot_path: specific protein path
        :param express_data: express data contain protein list
        :param stable_pairs_result: stable pairs result
        :return specific_p2: stable pair has two specific protein
        :return specific_p1: stable pair has a specific protein
        :return prot_new: protein list
        """

        specific_prot = specific_prot.loc[:,'specific_protein'].tolist()

        _data = express_data
        _result = stable_pairs_result

        prot_index = _data.index.tolist()
        prot_new = []

        for p in prot_index:
            prot_new.append(p.replace('.', '-'))

        specific_prot_filter = []
        specific_prot_filter_index = []
        for _sp in specific_prot:
            if _sp in prot_new:
                specific_prot_filter.append(_sp)
                specific_prot_filter_index.append(prot_new.index(_sp))

        _r_1 = _result.iloc[:,0].tolist()
        _r_2 = _result.iloc[:,1].tolist()
        t_1 = []
        t_2 = []
        for _t in specific_prot_filter_index:
            if _t in _r_1:
                t_1.append(_r_1.index(_t))
            if _t in _r_2:
                t_2.append(_r_2.index(_t))

        specific_p2 = _result.iloc[list(set(t_1).intersection(set(t_2))), :]
        specific_p1 = _result.iloc[list(set(t_1).union(set(t_2))), :]

        return specific_p2, specific_p1, prot_new

    def plot_result_random(self, normal_data, tumor_data, sp_reversed, sp_concordant, save_dir):
        """
        Visualize stable pairs by random
        :param normal_data: normal data
        :param tumor_data: tumor data
        :param sp_reversed: reversed pairs
        :param sp_concordant: concordant pairs
        :param save_dir: figure output directory
        """
        id_r = random.randint(0,sp_reversed.shape[0]-1)
        id_c = random.randint(0,sp_concordant.shape[0]-1)
        self.plot_pairs_rank(sp_reversed, normal_data, tumor_data, pair_id=id_r, save_dir=save_dir, if_reversed=True)
        self.plot_pairs_rank(sp_concordant, normal_data, tumor_data, pair_id=id_c, save_dir=save_dir, if_reversed=False)

    def concordant_result_select(self, specific_protein, normal_data, sp_concordant, NUM_VISUAL):
        """
        Select concordant pair to visualize
        :param specific_protein: specific protein
        :param normal_data: normal data
        :param sp_concordant: concordant pairs
        :param NUM_VISUAL: the number of pairs to visualize
        """
        sp_two, sp_one, sp_index = self.select_stable_pairs(specific_protein, normal_data, sp_concordant)
        print('\nDistribution of specific proteins in concordant pairs:')
        print('>>> Concordant pairs has two special proteins:  %s'%sp_two.shape[0])
        print('>>> Concordant pairs has one special proteins:  %s' %sp_one.shape[0])


        if sp_two.shape[0] > NUM_VISUAL:
            sp_two_select = sp_two.iloc[random.sample(range(0, sp_two.shape[0]), NUM_VISUAL), :]
        else:
            sp_two_select = sp_two

        if sp_one.shape[0] > NUM_VISUAL:
            sp_one_select = sp_one.iloc[random.sample(range(0, sp_one.shape[0]), NUM_VISUAL), :]
        else:
            sp_one_select = sp_one

        sp_two_select_prot = copy.deepcopy(sp_two_select)
        sp_one_select_prot = copy.deepcopy(sp_one_select)
        for idx in range(sp_two_select.shape[0]):
            for col in range(sp_two_select.shape[1]):
                sp_two_select_prot.iloc[idx, col] = sp_index[sp_two_select.iloc[idx, col]]

        for idx in range(sp_one_select.shape[0]):
            for col in range(sp_one_select.shape[1]):
                sp_one_select_prot.iloc[idx, col] = sp_index[sp_one_select.iloc[idx, col]]

        col_name = ['Protein-1', 'Protein-2']
        sp_two_select_prot.columns = col_name
        sp_one_select_prot.columns = col_name

        return sp_one_select_prot, sp_two_select_prot

    def reversed_result_select(self, specific_protein, normal_data, sp_reversed, NUM_VISUAL):
        """
        Select reversed pair to visualize
        :param specific_protein: specific protein path
        :param normal_data: normal data
        :param sp_concordant: concordant pairs
        :param NUM_VISUAL: the number of pairs to visualize
        """
        sp_two, sp_one, sp_index = self.select_stable_pairs(specific_protein, normal_data, sp_reversed)
        print('\nDistribution of specific proteins in reversed pairs:')
        print('>>> Reversal pairs has two special protein:  %s'%sp_two.shape[0])
        print('>>> Reversal pairs has one special protein:  %s' %sp_one.shape[0])

        if sp_two.shape[0] > NUM_VISUAL:
            sp_two_select = sp_two.iloc[random.sample(range(0, sp_two.shape[0]), NUM_VISUAL), :]
        else:
            sp_two_select = sp_two

        if sp_one.shape[0] > NUM_VISUAL:
            sp_one_select = sp_one.iloc[random.sample(range(0, sp_one.shape[0]), NUM_VISUAL), :]
        else:
            sp_one_select = sp_one

        sp_two_select_prot = copy.deepcopy(sp_two_select)
        sp_one_select_prot = copy.deepcopy(sp_one_select)
        for idx in range(sp_two_select.shape[0]):
            for col in range(sp_two_select.shape[1]):
                sp_two_select_prot.iloc[idx, col] = sp_index[sp_two_select.iloc[idx, col]]

        for idx in range(sp_one_select.shape[0]):
            for col in range(sp_one_select.shape[1]):
                sp_one_select_prot.iloc[idx, col] = sp_index[sp_one_select.iloc[idx, col]]

        col_name = ['Protein-1', 'Protein-2']
        sp_two_select_prot.columns = col_name
        sp_one_select_prot.columns = col_name

        return sp_one_select_prot, sp_two_select_prot


    def plot_result_select(self, normal_data, tumor_data, sp_reversed, sp_concordant, save_dir, specific_protein, NUM_VISUAL):
        """
        Select the result to be visualized
        :param normal_data
        :param tumor_data
        :param sp_reversed
        :param sp_concordant
        :param save_dir
        :param specific_protein
        :param NUM_VISUAL
        """
        mode = input('What do you want to visualize？ \n1.concordant pairs\n2.reversed pairs\n')

        if mode == '1':
            sp_one_select, sp_two_select = self.concordant_result_select(specific_protein, normal_data, sp_concordant, NUM_VISUAL)

            if sp_one_select.shape[0] == 0:
                print('\nWarning！There are no concordant pairs containing specific proteins.')
            else:
                if sp_one_select.shape[0] != 0 and sp_two_select.shape[0] == 0:
                    print('\nWarning！There are no concordant pairs containing two specific proteins.')               
                id_con = input("\nWhich concordant pairs do you want to visualize?\nPlease enter the number：\n\nEither\n%s\n\nBoth\n%s\n"%
                               (sp_one_select, sp_two_select))
                self.plot_pairs_rank(sp_concordant, normal_data, tumor_data, pair_id=int(id_con), save_dir=save_dir, if_reversed=False)

        elif mode == '2':
            sp_one_select, sp_two_select = self.reversed_result_select(specific_protein, normal_data, sp_reversed, NUM_VISUAL)

            if sp_one_select.shape[0] == 0:
                print('\nWarning！There are no reversal pairs containing specific proteins.')
            else:
                if sp_one_select.shape[0] != 0 and sp_two_select.shape[0] == 0:
                    print('\nWarning！There are no reversal pairs containing two specific proteins.')
                id_rev = input("\nWhich reversed pairs do you want to visualize?\nPlease enter the number：\n\nEither\n%s\n\nBoth\n%s\n"%
                               (sp_one_select, sp_two_select))
                self.plot_pairs_rank(sp_reversed, normal_data, tumor_data, pair_id=int(id_rev), save_dir=save_dir, if_reversed=True)
                
    def run_stablePairs(self):
        # 创建工作目录，若不存在则创建
        if os.path.exists(self.sp_workdir) == False:
            os.makedirs(self.sp_workdir)
        os.chdir(self.sp_workdir)

        # 数据读入
        data = self.rd.data_imput
        normal_cohort = self.rd.normal_cohort
        tumor_cohort = self.rd.tumor_cohort

        normal_data = self.rd.data_normal_imput
        tumor_data = self.rd.data_tumor_imput

        # 将数据写入硬盘
        # --------------------此处放入utils
        utils.write_normal_dat(normal_data, self.sp_workdir)
        utils.write_tumor_dat(tumor_data, self.sp_workdir)

        # run rankcomp j1
        utils.rankc_j1(self.reoa_path, normal_data, tumor_data, self.sp_workdir, self.CYCLE_RANKC, self.FDR, self.MAX_EXCEPTION_RANKC)

        # 读取rankcomp j1的结果
        sp_concordant, sp_reversed, sp_normal = utils.get_rankc_j1_results(self.sp_workdir)
        
        if_return = False
        print('>>> Normal stable pairs: %s'%sp_normal.shape[0])
        print('>>> Concordant pairs: %s'%sp_concordant.shape[0])
        print('>>> Reversed pairs: %s'%sp_reversed.shape[0])
        print('\nThe stable pairs of each type should be not less than %s!'%self.SP_THRES)
        if sp_normal.shape[0] < self.SP_THRES:
            print('\n    Error: The number of normal stable pairs less than %s!'%self.SP_THRES)
            print('\nYou should adjust the input parameters')
            if_return = True
        if sp_concordant.shape[0] < self.SP_THRES:
            print('\n    Error: The number of concordant pairs less than %s!'%self.SP_THRES)
            print('\nYou should adjust the input parameters')
            if_return = True
        if sp_reversed.shape[0] < self.SP_THRES:
            print('\n    Error: The number of reversed pairs less than %s!'%self.SP_THRES)
            print('\nYou should adjust the input parameters')
            if_return = True
            
        if self.fig_save_path == '':
            save_dir = self.sp_workdir
        else:
            save_dir = self.fig_save_path

        # 若找到的normal stable pairs, reversed, conocordant不少于SP_THRES。则将结果保存至本地
        if if_return:
            pass
        else:
            if self.RANDOM_VISUAL:
                self.plot_result_random(normal_data, tumor_data, sp_reversed, sp_concordant, save_dir)
            else:
                self.plot_result_select(normal_data, tumor_data, sp_reversed, sp_concordant, save_dir, self.rd.specific_protein_adapt, self.NUM_VISUAL)


