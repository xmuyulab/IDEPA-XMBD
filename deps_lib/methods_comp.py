
import pandas as pd
import numpy as np
import shutil
import random
import os
import copy

import sys
import datetime

import seaborn as sns
import matplotlib.pyplot as plt

from scipy import stats

from deps_lib import utils
from deps_lib import penda_pro
from deps_lib import methods_lib

from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector

r_stats = importr('stats')
        
class methodsComp(object):
    
    def __init__(self, rd, r_penda_path, r_penda_fdr_path,  reoa_path, 
                 CYCLE_RANKC, FDR, MAX_EXCEPTION_RANKC, PLOT_METHODS_COMPARE_RESULT, METHODS_LIST):
        
        self.rd = rd
        self.r_penda_path = r_penda_path
        self.r_penda_fdr_path = r_penda_fdr_path
        self.reoa_path = reoa_path
        self.CYCLE_RANKC = CYCLE_RANKC
        self.FDR = FDR
        self.MAX_EXCEPTION_RANKC = MAX_EXCEPTION_RANKC
        self.PLOT_METHODS_COMPARE_RESULT = PLOT_METHODS_COMPARE_RESULT
        self.METHODS_LIST = METHODS_LIST
        
 
    def method_comp(self, rankc_v1_up = None, rankc_v1_down = None, 
                    rankc_v2_up = None, rankc_v2_down = None, 
                    penda_up = None, penda_down = None, 
                    penda_up_fdr = None, penda_down_fdr = None, 
                    penda_up_pro = None, penda_down_pro = None,
                    peng_up = None, peng_down = None,
                    ttest_up = None, ttest_down = None, 
                    wilcox_up = None, wilcox_down = None, 
                    quantile_up = None, quantile_down = None,
                    val_normal = None, val_tumor = None, 
                    save_path = None, PLOT_METHODS_COMPARE_RESULT = False):
        """Show analysis results
        :param rankc_v1_up
        :param rankc_v1_down
        :param rankc_v2_up
        :param rankc_v2_down
        :param penda_up: Penda up deregulated protein.
        :param penda_down: Penda down deregulated protein.
        :param penda_up_fdr
        :param penda_down_fdr
        :param peng_up
        :param peng_down
        :param ttest_up
        :param ttest_down
        :param wilcox_up
        :param wilcox_down
        :param val_normal: Gold standard data(Tumor)
        :param val_tumor: Gold standard data(Normal)
        :param save_path: fig save path
        :param PLOT_METHODS_COMPARE_RESULT: Whether to show results.
        :return: Visualize the results.
        """

        # 数据读取，整理


        val_normal.columns = val_tumor.columns

        # 创建金标准数据
        gold_S = val_tumor - val_normal
        gold_S[gold_S>0] = 2
        gold_S[gold_S<0] = 1

        gold_S = gold_S.astype("int")
        
        precision_mess = {}
        positive_mess = {}

        ############## rankcomp v1,v2 ###############################
        if 'RankComp' in self.METHODS_LIST:

            # 计算统计值
            u_v1, u_v1_T, d_v1, d_v1_T , r_v1, r_v1_T = [], [], [], [], [], []
            u_v2, u_v2_T, d_v2, d_v2_T , r_v2, r_v2_T = [], [], [], [], [], []

            for col in rankc_v1_up.columns:
                u_v1.extend([rankc_v1_up[col].sum()])
                d_v1.extend([rankc_v1_down[col].sum()])
                r_v1.extend([rankc_v1_up[col].sum() + rankc_v1_down[col].sum()])

                u_v1_T.extend([(gold_S[rankc_v1_up] == 2).sum()[col]])
                d_v1_T.extend([(gold_S[rankc_v1_down] == 1).sum()[col]])
                r_v1_T.extend([(gold_S[rankc_v1_down] == 1).sum()[col] + (gold_S[rankc_v1_up] == 2).sum()[col]])
                
                u_v2.extend([rankc_v2_up[col].sum()])
                d_v2.extend([rankc_v2_down[col].sum()])
                r_v2.extend([rankc_v2_up[col].sum() + rankc_v2_down[col].sum()])

                u_v2_T.extend([(gold_S[rankc_v2_up] == 2).sum()[col]])
                d_v2_T.extend([(gold_S[rankc_v2_down] == 1).sum()[col]])
                r_v2_T.extend([(gold_S[rankc_v2_down] == 1).sum()[col] + (gold_S[rankc_v2_up] == 2).sum()[col]])


            # 计算precision
            prec_u_v1 = np.array(u_v1_T) / np.array(u_v1)
            prec_d_v1 = np.array(d_v1_T) / np.array(d_v1)

            prec_u_v2 = np.array(u_v2_T) / np.array(u_v2)
            prec_d_v2 = np.array(d_v2_T) / np.array(d_v2)

            prec_r_v1 = np.array(r_v1_T) / np.array(r_v1)
            prec_r_v2 = np.array(r_v2_T) / np.array(r_v2)
            
            precision_mess['rankc_v1'] = prec_r_v1
            positive_mess['rankc_v1'] = r_v1
            precision_mess['rankc_v2'] = prec_r_v2
            positive_mess['rankc_v2'] = r_v2

        ################## penda ########################################
        if 'Penda' in self.METHODS_LIST:
            u_penda, d_penda, u_penda_T, d_penda_T, r_penda, r_penda_T = [], [], [], [], [], []

            for col in penda_down.columns:
                u_penda.extend([penda_up[col].sum()])
                d_penda.extend([penda_down[col].sum()])
                r_penda.extend([penda_up[col].sum() + penda_down[col].sum()])

                u_penda_T.extend([(gold_S[penda_up] == 2).sum()[col]])
                d_penda_T.extend([(gold_S[penda_down] == 1).sum()[col]])
                r_penda_T.extend([(gold_S[penda_down] == 1).sum()[col] + (gold_S[penda_up] == 2).sum()[col]])

            # 计算precision
            prec_u_penda = np.array(u_penda_T) / np.array(u_penda)
            prec_d_penda = np.array(d_penda_T) / np.array(d_penda)
            prec_r_penda = np.array(r_penda_T) / np.array(r_penda)   
            
            precision_mess['penda'] = prec_r_penda
            positive_mess['penda'] = r_penda

        ################## penda_fdr ########################################
        if 'Penda fdr' in self.METHODS_LIST:
            u_penda_fdr, d_penda_fdr, u_penda_fdr_T, d_penda_fdr_T, r_penda_fdr, r_penda_fdr_T = [], [], [], [], [], []

            for col in penda_down_fdr.columns:
                u_penda_fdr.extend([penda_up_fdr[col].sum()])
                d_penda_fdr.extend([penda_down_fdr[col].sum()])
                r_penda_fdr.extend([penda_up_fdr[col].sum() + penda_down_fdr[col].sum()])

                u_penda_fdr_T.extend([(gold_S[penda_up_fdr] == 2).sum()[col]])
                d_penda_fdr_T.extend([(gold_S[penda_down_fdr] == 1).sum()[col]])
                r_penda_fdr_T.extend([(gold_S[penda_down_fdr] == 1).sum()[col] + (gold_S[penda_up_fdr] == 2).sum()[col]])

            # 计算precision
            prec_u_penda_fdr = np.array(u_penda_fdr_T) / np.array(u_penda_fdr)
            prec_d_penda_fdr = np.array(d_penda_fdr_T) / np.array(d_penda_fdr)
            prec_r_penda_fdr = np.array(r_penda_fdr_T) / np.array(r_penda_fdr) 

            precision_mess['penda_fdr'] = prec_r_penda_fdr
            positive_mess['penda_fdr'] = r_penda_fdr

        ################## penda_pro ########################################
        if 'Penda pro' in self.METHODS_LIST:
            u_penda_pro, d_penda_pro, u_penda_pro_T, d_penda_pro_T, r_penda_pro, r_penda_pro_T = [], [], [], [], [], []

            for col in penda_down_pro.columns:
                u_penda_pro.extend([penda_up_pro[col].sum()])
                d_penda_pro.extend([penda_down_pro[col].sum()])
                r_penda_pro.extend([penda_up_pro[col].sum() + penda_down_pro[col].sum()])

                u_penda_pro_T.extend([(gold_S[penda_up_pro] == 2).sum()[col]])
                d_penda_pro_T.extend([(gold_S[penda_down_pro] == 1).sum()[col]])
                r_penda_pro_T.extend([(gold_S[penda_down_pro] == 1).sum()[col] + (gold_S[penda_up_pro] == 2).sum()[col]])

            # 计算precision
            prec_u_penda_pro = np.array(u_penda_pro_T) / np.array(u_penda_pro)
            prec_d_penda_pro = np.array(d_penda_pro_T) / np.array(d_penda_pro)
            prec_r_penda_pro = np.array(r_penda_pro_T) / np.array(r_penda_pro) 
            
            precision_mess['penda_pro'] = prec_r_penda_pro
            positive_mess['penda_pro'] = r_penda_pro
            
            
        ################## peng method ########################################
        if 'Peng method' in self.METHODS_LIST:
            u_peng, d_peng, u_peng_T, d_peng_T, r_peng, r_peng_T = [], [], [], [], [], []

            for col in peng_down.columns:
                u_peng.extend([peng_up[col].sum()])
                d_peng.extend([peng_down[col].sum()])
                r_peng.extend([peng_up[col].sum() + peng_down[col].sum()])

                u_peng_T.extend([(gold_S[peng_up] == 2).sum()[col]])
                d_peng_T.extend([(gold_S[peng_down] == 1).sum()[col]])
                r_peng_T.extend([(gold_S[peng_down] == 1).sum()[col] + (gold_S[peng_up] == 2).sum()[col]])

            # 计算precision
            prec_u_peng = np.array(u_peng_T) / np.array(u_peng)
            prec_d_peng = np.array(d_peng_T) / np.array(d_peng)
            prec_r_peng = np.array(r_peng_T) / np.array(r_peng)   
            
            precision_mess['peng_method'] = prec_r_peng
            positive_mess['peng_method'] = r_peng
            

        ################# ttest #######################################
        if 'T-test' in self.METHODS_LIST:
            u_tt, d_tt, u_tt_T, d_tt_T, r_tt, r_tt_T = [], [], [], [], [], []

            for col in ttest_down.columns:

                u_tt.extend([ttest_up[col].sum()])
                d_tt.extend([ttest_down[col].sum()])
                r_tt.extend([ttest_up[col].sum() + ttest_down[col].sum()])

                u_tt_T.extend([(gold_S[ttest_up == 1] == 2).sum()[col]])
                d_tt_T.extend([(gold_S[ttest_down == 1] == 1).sum()[col]])
                r_tt_T.extend([(gold_S[ttest_down == 1] == 1).sum()[col] + (gold_S[ttest_up == 1] == 2).sum()[col]])

            # 计算precision
            prec_u_tt = np.array(u_tt_T) / np.array(u_tt)
            prec_d_tt = np.array(d_tt_T) / np.array(d_tt)
            prec_r_tt = np.array(r_tt_T) / np.array(r_tt)
            
            precision_mess['ttest'] = prec_r_tt
            positive_mess['ttest'] = r_tt

        ################# wilcox #######################################
        if 'Wilcoxon' in self.METHODS_LIST:
            u_wilc, d_wilc, u_wilc_T, d_wilc_T, r_wilc, r_wilc_T = [], [], [], [], [], []

            for col in wilcox_down.columns:

                u_wilc.extend([wilcox_up[col].sum()])
                d_wilc.extend([wilcox_down[col].sum()])
                r_wilc.extend([wilcox_up[col].sum() + wilcox_down[col].sum()])

                u_wilc_T.extend([(gold_S[wilcox_up == 1] == 2).sum()[col]])
                d_wilc_T.extend([(gold_S[wilcox_down == 1] == 1).sum()[col]])
                r_wilc_T.extend([(gold_S[wilcox_down == 1] == 1).sum()[col] + (gold_S[wilcox_up == 1] == 2).sum()[col]])

            # 计算precision
            prec_u_wilc = np.array(u_wilc_T) / np.array(u_wilc)
            prec_d_wilc = np.array(d_wilc_T) / np.array(d_wilc)
            prec_r_wilc = np.array(r_wilc_T) / np.array(r_wilc)
            
            precision_mess['wilcox'] = prec_r_wilc
            positive_mess['wilcox'] = r_wilc

        ################# quantile #######################################
        if 'Quantile' in self.METHODS_LIST:
            u_quant, d_quant, u_quant_T, d_quant_T, r_quant, r_quant_T = [], [], [], [], [], []

            for col in quantile_down.columns:

                u_quant.extend([quantile_up[col].sum()])
                d_quant.extend([quantile_down[col].sum()])
                r_quant.extend([quantile_up[col].sum() + quantile_down[col].sum()])

                u_quant_T.extend([(gold_S[quantile_up == 1] == 2).sum()[col]])
                d_quant_T.extend([(gold_S[quantile_down == 1] == 1).sum()[col]])
                r_quant_T.extend([(gold_S[quantile_down == 1] == 1).sum()[col] + (gold_S[quantile_up == 1] == 2).sum()[col]])

            # 计算precision
            prec_u_quant = np.array(u_quant_T) / np.array(u_quant)
            prec_d_quant = np.array(d_quant_T) / np.array(d_quant)
            prec_r_quant = np.array(r_quant_T) / np.array(r_quant)
            
            precision_mess['quantile'] = prec_r_quant
            positive_mess['quantile'] = r_quant


        ################ plot result ################################

        plt.rcParams['font.sans-serif'] = ['SimHei']  # 中文字体设置-黑体
        sns.set(font_scale=1.5)  # 解决Seaborn中文显示问题并调整字体大小
        plt.rcParams['savefig.dpi'] = 600 #图片像素

        plt.rcParams['figure.dpi'] = 300 #分辨率
        plt.rcParams['axes.unicode_minus']=False # 用来显示负号



        # 画图
        x_list = np.arange(1,val_tumor.shape[1]+1)
        plt.figure(figsize=(24,18))
        plt.subplot(2,1,1)
        if 'RankComp' in self.METHODS_LIST:
            plt.scatter(x_list,prec_u_v1, marker='^', color = 'g', s=150, label='U_v1: %.4f'%np.mean(prec_u_v1))
            plt.scatter(x_list,prec_u_v2, marker='^', color = 'r', s=150, label='U_v2: %.4f'%np.mean(prec_u_v2))
        if 'Penda' in self.METHODS_LIST:
            plt.scatter(x_list,prec_u_penda, marker='^', color = 'y', s=150, label='U_penda: %.4f'%np.mean(prec_u_penda))
        if 'Penda fdr' in self.METHODS_LIST:
            plt.scatter(x_list,prec_u_penda_fdr, marker='^', color = 'm', s=150, label='U_penda_fdr: %.4f'%np.mean(prec_u_penda_fdr))
        if 'Penda pro' in self.METHODS_LIST:
            plt.scatter(x_list,prec_u_penda_pro, marker='^', color = 'c', s=150, label='U_penda_pro: %.4f'%np.mean(prec_u_penda_pro))
        if 'Peng method' in self.METHODS_LIST:
            plt.scatter(x_list,prec_u_peng, marker='^', color = 'lime', s=150, label='U_peng: %.4f'%np.mean(prec_u_peng))
        if 'T-test' in self.METHODS_LIST:
            plt.scatter(x_list,prec_u_tt, marker='^', color = 'k', s=150, label='U_ttest: %.4f'%np.mean(prec_u_tt))
        if 'Wilcoxon' in self.METHODS_LIST:
            plt.scatter(x_list,prec_u_wilc, color = 'b', marker='^',s=150, label='U_wilc: %.4f'%np.mean(prec_u_wilc))
        if 'Quantile' in self.METHODS_LIST:
            plt.scatter(x_list,prec_u_quant, color='gold', marker='^',s=150, label='U_quant: %.4f'%np.mean(prec_u_quant))

        for i in x_list:
            plt.vlines(i, 0.65, 1,color="darkgrey",linestyles='dotted')#竖线

        plt.xlabel('samples', fontsize=30)
        plt.ylabel('precision', fontsize=30)
        plt.title("The precisions of DE protein: UP", fontsize=35)
        # plt.title("The precisions of DE protein by RankComp(v1-v2)-Penda-Wilcox: %s"%title, fontsize=25)
        plt.legend(loc=[1.05, 0], fontsize=25)
        plt.xticks(fontsize = 20)
        plt.yticks(fontsize = 20)


        plt.subplot(2,1,2)
        if 'RankComp' in self.METHODS_LIST:
            plt.plot(x_list, u_v1, linestyle = '-', color = 'g', marker = '^', markersize = 10, linewidth = 4, label='U_v1: %d'%np.mean(u_v1))
            plt.plot(x_list, u_v2, linestyle = '-', color = 'r', marker = '^', markersize = 10, linewidth = 4, label='U_v2: %d'%np.mean(u_v2))
        if 'Penda' in self.METHODS_LIST:
            plt.plot(x_list, u_penda, linestyle = '-', color = 'y', marker = '^', markersize = 10, linewidth = 4, label='U_penda: %d'%np.mean(u_penda))
        if 'Penda fdr' in self.METHODS_LIST:
            plt.plot(x_list, u_penda_fdr, linestyle = '-', color = 'm', marker = '^', markersize = 10, linewidth = 4, label='U_penda_fdr: %d'%np.mean(u_penda_fdr))
        if 'Penda pro' in self.METHODS_LIST:
            plt.plot(x_list, u_penda_pro, linestyle = '-', color = 'c', marker = '^', markersize = 10, linewidth = 4, label='U_penda_pro: %d'%np.mean(u_penda_pro))
        if 'Peng method' in self.METHODS_LIST:
            plt.plot(x_list, u_peng, linestyle = '-', color = 'lime', marker = '^', markersize = 10, linewidth = 4, label='U_peng: %d'%np.mean(u_peng))
        if 'T-test' in self.METHODS_LIST:
            plt.plot(x_list, u_tt, linestyle = '-', color = 'k', marker = '^', markersize = 10, linewidth = 4, label='U_ttest: %d'%np.mean(u_tt))
        if 'Wilcoxon' in self.METHODS_LIST:
            plt.plot(x_list, u_wilc, linestyle = '-', color = 'b', marker = '^', markersize = 10, linewidth = 4, label='U_wilc: %d'%np.mean(u_wilc))
        if 'Quantile' in self.METHODS_LIST:
            plt.plot(x_list, u_quant, linestyle = '-', color='gold', marker = '^', markersize = 10, linewidth = 4, label='U_quant: %d'%np.mean(u_quant))


        plt.xlabel('samples', fontsize=30)
        plt.ylabel('number', fontsize=30)
        plt.title("The number of DE protein: UP", fontsize=35)
        plt.legend(loc=[1.05, 0], fontsize=25)
        plt.xticks(fontsize = 20)
        plt.yticks(fontsize = 20)
        plt.tight_layout()

        plt.savefig(save_path+'/up_deregulated.pdf', dpi=800, format='pdf', bbox_inches='tight') #指定分辨
        if PLOT_METHODS_COMPARE_RESULT:
            plt.show()

        plt.figure(figsize=(24,18))
        plt.subplot(2,1,1)
        if 'RankComp' in self.METHODS_LIST:
            plt.scatter(x_list,prec_d_v1, marker='v', color = 'g', s=150, label='D_v1: %.4f'%np.mean(prec_d_v1))
            plt.scatter(x_list,prec_d_v2, marker='v', color = 'r', s=150, label='D_v2: %.4f'%np.mean(prec_d_v2))
        if 'Penda' in self.METHODS_LIST:
            plt.scatter(x_list,prec_d_penda, marker='v', color = 'y', s=150, label='D_penda: %.4f'%np.mean(prec_d_penda))
        if 'Penda fdr' in self.METHODS_LIST:
            plt.scatter(x_list,prec_d_penda_fdr, marker='v', color = 'm', s=150, label='D_penda_fdr: %.4f'%np.mean(prec_d_penda_fdr))
        if 'Penda pro' in self.METHODS_LIST:
            plt.scatter(x_list,prec_d_penda_pro, marker='v', color = 'c', s=150, label='D_penda_pro: %.4f'%np.mean(prec_d_penda_pro))
        if 'Peng method' in self.METHODS_LIST:
            plt.scatter(x_list,prec_d_peng, marker='v', color = 'lime', s=150, label='D_peng: %.4f'%np.mean(prec_d_peng))
        if 'T-test' in self.METHODS_LIST:
            plt.scatter(x_list,prec_d_tt, marker='v', color = 'k', s=150, label='D_ttest: %.4f'%np.mean(prec_d_tt))
        if 'Wilcoxon' in self.METHODS_LIST:
            plt.scatter(x_list,prec_d_wilc, marker='v', color = 'b', s=150, label='D_wilc: %.4f'%np.mean(prec_d_wilc))
        if 'Quantile' in self.METHODS_LIST:
            plt.scatter(x_list,prec_d_quant, marker='v', color='gold', s=150, label='D_quant: %.4f'%np.mean(prec_d_quant))

        for i in x_list:
            plt.vlines(i, 0.65, 1,color="darkgrey",linestyles='dotted')#竖线

        plt.xlabel('samples', fontsize=30)
        plt.ylabel('precision', fontsize=30)
        plt.title("The precisions of DE protein: DOWN", fontsize=35)
        # plt.title("The precisions of DE protein by RankComp(v1-v2)-Penda-Wilcox: %s"%title, fontsize=25)
        plt.legend(loc=[1.05, 0], fontsize=25)
        plt.xticks(fontsize = 20)
        plt.yticks(fontsize = 20)

        plt.subplot(2,1,2)
        if 'RankComp' in self.METHODS_LIST:
            plt.plot(x_list, d_v1, linestyle = '-', color = 'g', marker = 'v', markersize = 10, linewidth = 4, label='D_v1: %d'%np.mean(d_v1))
            plt.plot(x_list, d_v2, linestyle = '-', color = 'r', marker = 'v', markersize = 10,  linewidth = 4, label='D_v2: %d'%np.mean(d_v2))
        if 'Penda' in self.METHODS_LIST:
            plt.plot(x_list, d_penda, linestyle = '-', color = 'y', marker = 'v', markersize = 10, linewidth = 4, label='D_penda: %d'%np.mean(d_penda))
        if 'Penda fdr' in self.METHODS_LIST:
            plt.plot(x_list, d_penda_fdr, linestyle = '-', color = 'm', marker = 'v', markersize = 10, linewidth = 4, label='D_penda_fdr: %d'%np.mean(d_penda_fdr))
        if 'Penda pro' in self.METHODS_LIST:
            plt.plot(x_list, d_penda_pro, linestyle = '-', color = 'c', marker = 'v', markersize = 10, linewidth = 4, label='D_penda_pro: %d'%np.mean(d_penda_pro))
        if 'Peng method' in self.METHODS_LIST:
            plt.plot(x_list, d_peng, linestyle = '-', color = 'lime', marker = 'v', markersize = 10, linewidth = 4, label='D_peng: %d'%np.mean(d_peng))
        if 'T-test' in self.METHODS_LIST:
            plt.plot(x_list, d_tt, linestyle = '-', color = 'k', marker = 'v', markersize = 10, linewidth = 4, label='D_ttest: %d'%np.mean(d_tt))
        if 'Wilcoxon' in self.METHODS_LIST:
            plt.plot(x_list, d_wilc, linestyle = '-', color = 'b', marker = 'v', markersize = 10, linewidth = 4, label='D_wilc: %d'%np.mean(d_wilc))
        if 'Quantile' in self.METHODS_LIST:
            plt.plot(x_list, d_quant, linestyle = '-', color='gold', marker = 'v', markersize = 10, linewidth = 4, label='D_quant: %d'%np.mean(d_quant))
            


        plt.xlabel('samples', fontsize=30)
        plt.ylabel('number', fontsize=30)
        plt.title("The number of DE protein: DOWN", fontsize=35)
        plt.legend(loc=[1.05, 0], fontsize=25)
        plt.xticks(fontsize = 20)
        plt.yticks(fontsize = 20)
        plt.tight_layout()

        plt.savefig(save_path+'/down_deregulated.pdf', dpi=800, format='pdf', bbox_inches='tight') #指定分辨
        if PLOT_METHODS_COMPARE_RESULT:
            plt.show()

        plt.figure(figsize=(24,18))
        plt.subplot(2,1,1)
        if 'RankComp' in self.METHODS_LIST:
            plt.scatter(x_list,prec_r_v1, marker='o', color = 'g', s=150, label='R_v1: %.4f'%np.mean(prec_r_v1))
            plt.scatter(x_list,prec_r_v2, marker='o', color = 'r', s=150, label='R_v2: %.4f'%np.mean(prec_r_v2))
        if 'Penda' in self.METHODS_LIST:
            plt.scatter(x_list,prec_r_penda, marker='o', color = 'y', s=150, label='R_penda: %.4f'%np.mean(prec_r_penda))
        if 'Penda fdr' in self.METHODS_LIST:
            plt.scatter(x_list,prec_r_penda_fdr, marker='o', color = 'm', s=150, label='R_penda_fdr: %.4f'%np.mean(prec_r_penda_fdr))
        if 'Penda pro' in self.METHODS_LIST:
            plt.scatter(x_list,prec_r_penda_pro, marker='o', color = 'c', s=150, label='R_penda_pro: %.4f'%np.mean(prec_r_penda_pro))
        if 'Peng method' in self.METHODS_LIST:
            plt.scatter(x_list,prec_r_peng, marker='o', color = 'lime', s=150, label='R_peng: %.4f'%np.mean(prec_r_peng))
        if 'T-test' in self.METHODS_LIST:
            plt.scatter(x_list,prec_r_tt, marker='o', color = 'k', s=150, label='R_ttest: %.4f'%np.mean(prec_r_tt))
        if 'Wilcoxon' in self.METHODS_LIST:
            plt.scatter(x_list,prec_r_wilc, marker='o', color = 'b', s=150, label='R_wilc: %.4f'%np.mean(prec_r_wilc))
        if 'Quantile' in self.METHODS_LIST:
            plt.scatter(x_list,prec_r_quant, marker='o', color='gold', s=150, label='R_quant: %.4f'%np.mean(prec_r_quant))

        for i in x_list:
            plt.vlines(i, 0.65, 1,color="darkgrey",linestyles='dotted')#竖线

        plt.xlabel('samples', fontsize=30)
        plt.ylabel('precision', fontsize=30)
        plt.title("The precisions of DE protein: Deregulation", fontsize=35)
        # plt.title("The precisions of DE protein by RankComp(v1-v2)-Penda-Wilcox: %s"%title, fontsize=25)
        plt.legend(loc=[1.05, 0], fontsize=25)
        plt.xticks(fontsize = 20)
        plt.yticks(fontsize = 20)

        plt.subplot(2,1,2)
        if 'RankComp' in self.METHODS_LIST:
            plt.plot(x_list, r_v1, linestyle = '-', color = 'g', marker = 'o', markersize = 10, linewidth = 4, label='R_v1: %d'%np.mean(r_v1))
            plt.plot(x_list, r_v2, linestyle = '-', color = 'r', marker = 'o', markersize = 10,  linewidth = 4, label='R_v2: %d'%np.mean(r_v2))
        if 'Penda' in self.METHODS_LIST:
            plt.plot(x_list, r_penda, linestyle = '-', color = 'y', marker = 'o', markersize = 10, linewidth = 4, label='R_penda: %d'%np.mean(r_penda))
        if 'Penda fdr' in self.METHODS_LIST:
            plt.plot(x_list, r_penda_fdr, linestyle = '-', color = 'm', marker = 'o', markersize = 10, linewidth = 4, label='R_penda_fdr: %d'%np.mean(r_penda_fdr))
        if 'Penda pro' in self.METHODS_LIST:
            plt.plot(x_list, r_penda_pro, linestyle = '-', color = 'c', marker = 'o', markersize = 10, linewidth = 4, label='R_penda_pro: %d'%np.mean(r_penda_pro))
        if 'Peng method' in self.METHODS_LIST:
            plt.plot(x_list, r_peng, linestyle = '-', color = 'lime', marker = 'o', markersize = 10, linewidth = 4, label='R_peng: %d'%np.mean(r_peng))
        if 'T-test' in self.METHODS_LIST:
            plt.plot(x_list, r_tt, linestyle = '-', color = 'k', marker = 'o', markersize = 10, linewidth = 4, label='R_ttest: %d'%np.mean(r_tt))
        if 'Wilcoxon' in self.METHODS_LIST:
            plt.plot(x_list, r_wilc, linestyle = '-', color = 'b', marker = 'o', markersize = 10, linewidth = 4, label='R_wilc: %d'%np.mean(r_wilc))
        if 'Quantile' in self.METHODS_LIST:
            plt.plot(x_list, r_quant, linestyle = '-', color = 'gold', marker = 'o', markersize = 10, linewidth = 4, label='R_quant: %d'%np.mean(r_quant))
            


        plt.xlabel('samples', fontsize=30)
        plt.ylabel('number', fontsize=30)
        plt.title("The number of DE protein: Deregulated", fontsize=35)
        plt.legend(loc=[1.05, 0], fontsize=25)
        plt.xticks(fontsize = 20)
        plt.yticks(fontsize = 20)
        plt.tight_layout()

        plt.savefig(save_path+'/deregulated.pdf', dpi=800, format='pdf', bbox_inches='tight') #指定分辨
        if PLOT_METHODS_COMPARE_RESULT:
            plt.show()
            
        return precision_mess, positive_mess

    def run_methodsComp(self, mc_workdir, method_comp_label=True):
        
        print('########## Comparison of methods ################')
        if os.path.exists(mc_workdir) == False:
            os.makedirs(mc_workdir)
        os.chdir(mc_workdir)

        # 数据读入
        data = self.rd.data_imput
        normal_cohort = self.rd.normal_cohort
        tumor_cohort = self.rd.tumor_cohort

        normal_run_data = self.rd.data_normal_imput
        tumor_run_data = self.rd.data_tumor_imput

        normal_val_data = self.rd.paired_data_normal_imput
        tumor_val_data = self.rd.paired_data_tumor_imput

        ################# Rankcomp ####################
#         if 'RankComp' in self.METHODS_LIST:
#             t = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
#             print('# Rankcom %s'%(t))
#             rankc_workspace = mc_workdir + '/rankcomp'
#             if os.path.exists(rankc_workspace) == True:
#                 shutil.rmtree(rankc_workspace)
#             os.makedirs(rankc_workspace)
#             os.chdir(rankc_workspace)

#             # 将rankcom输入数据写入硬盘
#             utils.write_normal_dat(normal_run_data, rankc_workspace)
#             utils.write_tumor_dat(tumor_run_data, rankc_workspace)

#             # run rankcomp v1 and v2
#             utils.rankc_j2(self.reoa_path, normal_run_data, tumor_run_data, rankc_workspace, 
#                            self.CYCLE_RANKC, self.FDR, self.MAX_EXCEPTION_RANKC)
#             # get rankcomp result
#             rankc_result = utils.get_rankc_j2_results(rankc_workspace)
#             rankc_result = rankc_result.iloc[:, :-1]
#             rankc_result.index = tumor_run_data.index
#             rankc_result.columns = tumor_run_data.columns
#         else:
#             rankc_result = None
            
        ################# Rankcomp ####################
        if 'RankComp' in self.METHODS_LIST:
            t = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print('# Rankcom %s'%(t))
            rankc_workspace = mc_workdir + '/rankcomp'
            if os.path.exists(rankc_workspace) == True:
                shutil.rmtree(rankc_workspace)
            os.makedirs(rankc_workspace)
            os.chdir(rankc_workspace)

            # 将rankcom输入数据写入硬盘
            utils.write_normal_dat(normal_run_data, rankc_workspace)
            utils.write_tumor_dat(tumor_run_data, rankc_workspace)

            # run rankcomp v1 and v2
            utils.rankc_j2_qvalues(self.reoa_path, normal_run_data, tumor_run_data, rankc_workspace, 
                                   self.CYCLE_RANKC, self.FDR, self.MAX_EXCEPTION_RANKC)
            # get rankcomp result
            rankc_v1_up, rankc_v1_down, rankc_v2_up, rankc_v2_down = utils.get_rankc_j2_qvalues_result(rankc_workspace, tumor_run_data)

        else:
            rankc_result = None

        ################ Penda #########################
        if 'Penda' in self.METHODS_LIST:
            t = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print('# Penda %s'%(t))
            penda_workspace = mc_workdir + '/penda'
            if os.path.exists(penda_workspace) == True:
                shutil.rmtree(penda_workspace)
            os.makedirs(penda_workspace)
            os.chdir(penda_workspace)

            # run penda
            utils.write_penda_data(normal_run_data, tumor_run_data, penda_workspace)
            utils.run_penda(self.r_penda_path)
            penda_up, penda_down, penda_result = utils.get_penda_result(penda_workspace)
            penda_up = penda_up.loc[tumor_run_data.index.tolist(),tumor_run_data.columns.tolist()]
            penda_down = penda_down.loc[tumor_run_data.index.tolist(), tumor_run_data.columns.tolist()]
            penda_result = penda_result.loc[tumor_run_data.index.tolist(),tumor_run_data.columns.tolist()]
        else:
            penda_result = None
            penda_down = None
            penda_up = None
        
        ################ Penda_fdr #########################
        if 'Penda fdr' in self.METHODS_LIST:
            t = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print('# Penda fdr %s'%(t))
            penda_fdr_workspace = mc_workdir + '/penda_fdr'
            if os.path.exists(penda_fdr_workspace) == True:
                shutil.rmtree(penda_fdr_workspace)
            os.makedirs(penda_fdr_workspace)
            os.chdir(penda_fdr_workspace)

            # run penda_fdr
            utils.write_penda_data(normal_run_data, tumor_run_data, penda_fdr_workspace)
            utils.run_penda_fdr(self.r_penda_fdr_path, self.FDR)
            penda_up_fdr, penda_down_fdr, penda_result_fdr = utils.get_penda_fdr_result(penda_fdr_workspace)
            penda_up_fdr = penda_up_fdr.loc[tumor_run_data.index.tolist(),tumor_run_data.columns.tolist()]
            penda_down_fdr = penda_down_fdr.loc[tumor_run_data.index.tolist(), tumor_run_data.columns.tolist()]
            penda_result_fdr = penda_result_fdr.loc[tumor_run_data.index.tolist(),tumor_run_data.columns.tolist()]
        else:
            penda_up_fdr = None
            penda_down_fdr = None
            penda_result_fdr = None

        ################ Penda pro #########################
        if 'Penda pro' in self.METHODS_LIST:
            t = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print('# Penda pro %s'%(t))
            penda_pro_workspace = mc_workdir + '/penda_pro'
            if os.path.exists(penda_pro_workspace) == True:
                shutil.rmtree(penda_pro_workspace)
            os.makedirs(penda_pro_workspace)
            os.chdir(penda_pro_workspace)

            # run penda pro
            utils.write_penda_data(normal_run_data, tumor_run_data, penda_pro_workspace)
            penda_up_pro_qv, penda_down_pro_qv = penda_pro.run_penda_pro(normal_path = './normal_run.csv', 
                                                                         tumor_path = './tumor_run.csv', 
                                                                         FDR = 0.05, 
                                                                         CONVERGENCE_THRESHOLD = 0.95, 
                                                                         MAX_CYCLE = 48, 
                                                                         K = 20, 
                                                                         THRESHOLD_LH = 0.99)
            penda_up_pro_qv.to_csv('./penda_up_pro_qv.csv')
            penda_down_pro_qv.to_csv('./penda_down_pro_qv.csv')
            penda_up_pro = penda_up_pro_qv < self.FDR
            penda_down_pro = penda_down_pro_qv < self.FDR

            penda_up_pro = penda_up_pro.loc[tumor_run_data.index.tolist(),tumor_run_data.columns.tolist()]
            penda_down_pro = penda_down_pro.loc[tumor_run_data.index.tolist(), tumor_run_data.columns.tolist()]
            # penda_result_pro = penda_up_pro + penda_down_pro
        else:
            penda_up_pro = None
            penda_down_pro = None
        
        
        ################ Peng_method #########################
        if 'Peng method' in self.METHODS_LIST:
            t = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print('# Peng method %s'%(t))
            peng_method_workspace = mc_workdir + '/peng_method'
            if os.path.exists(peng_method_workspace) == True:
                shutil.rmtree(peng_method_workspace)
            os.makedirs(peng_method_workspace)
            os.chdir(peng_method_workspace)

            # run penda
            utils.write_penda_data(normal_run_data, tumor_run_data, penda_workspace)
            peng_up, peng_down = methods_lib.run_peng_method(normal_run_data, tumor_run_data)
            penda_up = penda_up.loc[tumor_run_data.index.tolist(),tumor_run_data.columns.tolist()]
            penda_down = penda_down.loc[tumor_run_data.index.tolist(), tumor_run_data.columns.tolist()]
            peng_up.to_csv('./peng_up.csv')
            peng_down.to_csv('./peng_down.csv')
        else:
            peng_down = None
            peng_up = None
        
        
        ################ T-test #######################
        if 'T-test' in self.METHODS_LIST:
            t = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print('# T-test %s'%(t))
            ttest_workspace = mc_workdir + '/ttest'
            if os.path.exists(ttest_workspace) == True:
                shutil.rmtree(ttest_workspace)
            os.makedirs(ttest_workspace)
            os.chdir(ttest_workspace)

            utils.write_penda_data(normal_run_data, tumor_run_data, ttest_workspace)

            ttest_up, ttest_down = methods_lib.run_ttest(normal_run_path = './normal_run.csv',
                                                         tumor_test_path = './tumor_run.csv',
                                                         workdir = ttest_workspace,
                                                         Q_THRES = self.FDR,
                                                         D_THRES = 0,
                                                         SAVE_OUT = True)
        else:
            ttest_up = None
            ttest_down = None

        ############## wilcoxon ######################
        if 'Wilcoxon' in self.METHODS_LIST:
            t = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print('# Wilcoxon %s'%(t))
            wilcox_workspace = mc_workdir + '/wilcoxon'
            if os.path.exists(wilcox_workspace) == True:
                shutil.rmtree(wilcox_workspace)
            os.makedirs(wilcox_workspace)
            os.chdir(wilcox_workspace)

            utils.write_penda_data(normal_run_data, tumor_run_data, wilcox_workspace)

            wilcox_up, wilcox_down = methods_lib.run_wilcox(normal_run_path = './normal_run.csv', 
                                                            tumor_test_path = './tumor_run.csv',
                                                            workdir = wilcox_workspace, 
                                                            Q_THRES = self.FDR, D_THRES = 0,
                                                            SAVE_OUT = True)
        else:
            wilcox_up = None
            wilcox_down = None
            
        ############## quantile ######################
        if 'Quantile' in self.METHODS_LIST:
            t = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print('# quantile %s'%(t))
            quantile_workspace = mc_workdir + '/quantile'
            if os.path.exists(quantile_workspace) == True:
                shutil.rmtree(quantile_workspace)
            os.makedirs(quantile_workspace)
            os.chdir(quantile_workspace)

            utils.write_penda_data(normal_run_data, tumor_run_data, quantile_workspace)

            quantile_up, quantile_down = methods_lib.quantile_dep(normal_data_path = './normal_run.csv',
                                                                  tumor_data_path = './tumor_run.csv',
                                                                  workdir = quantile_workspace, 
                                                                  quantile = 0.05, 
                                                                  factor = 1.2,
                                                                  SAVE_OUT = True)
        else:
            quantile_up = None
            quantile_down = None

            
        ############## Method comparison ######################
        if method_comp_label:
            t = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print('\n# method comparison %s'%(t))

            precision_mess, positive_mess = self.method_comp(rankc_v1_up = rankc_v1_up, rankc_v1_down = rankc_v1_down, 
                                                             rankc_v2_up = rankc_v2_up, rankc_v2_down = rankc_v2_down,
                                                             penda_up = penda_up, penda_down = penda_down, 
                                                             penda_up_fdr = penda_up_fdr, penda_down_fdr = penda_down_fdr, 
                                                             penda_up_pro = penda_up_pro, penda_down_pro = penda_down_pro,
                                                             peng_up = peng_up, peng_down = peng_down,
                                                             ttest_up = ttest_up, ttest_down = ttest_down,
                                                             wilcox_up = wilcox_up, wilcox_down = wilcox_down,
                                                             quantile_up = quantile_up, quantile_down = quantile_down,
                                                             val_normal = normal_val_data, val_tumor = tumor_val_data,
                                                             save_path = mc_workdir, 
                                                             PLOT_METHODS_COMPARE_RESULT = self.PLOT_METHODS_COMPARE_RESULT)
        
            return precision_mess, positive_mess