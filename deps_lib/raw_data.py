#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np
import pandas as pd
import os
import copy

from scipy import stats

from deps_lib import utils


class InputData(object):
    """
    用于储存输入数据，并且对数据进行预处理
    """
    def __init__(self, data=None, normal_cohort=None, tumor_cohort=None, 
                 specific_protein=None, HAS_SPECIFIC_PROTEIN=None,
                 paired_data=None, paired_samples=None, HAS_PAIRED_DATA=False,  
                 INDEX_OLD_CHAR=None, INDEX_NEW_CHAR=None, COLUMNS_OLD_CHAR=None, COLUMNS_NEW_CHAR=None,
                 NORMALIZATION=False, LOG_LABEL=False, IMPUT_LABEL=False, NA_LABEL=None, NA_RATE=0.3):
        
        # unpaired data
        self.data = data
        self.normal_cohort = normal_cohort
        self.tumor_cohort = tumor_cohort
        self.data_imput = None
        self.data_normal_imput = None
        self.data_tumor_imput = None
        self.normal_cohort_adapt = None
        self.tumor_cohort_adapt = None
        
        # stable pairs
        self.HAS_SPECIFIC_PROTEIN = HAS_SPECIFIC_PROTEIN
        self.specific_protein = specific_protein
        self.specific_protein_adapt = None
        
        # normal cohort and tumor cohort is paired
        self.HAS_PAIRED_DATA = HAS_PAIRED_DATA
        self.paired_data = paired_data
        self.paired_samples = paired_samples
        self.paired_data_imput = None
        self.paired_data_normal_cohort = None
        self.paired_data_tumor_cohort = None
        self.paired_normal_cohort_adapt = None
        self.paired_tumor_cohort_adapt = None
        self.paired_data_normal_imput = None
        self.paired_data_tumor_imput = None
        
        # preprocess 
        self.NORMALIZATION = NORMALIZATION
        self.LOG_LABEL = LOG_LABEL
        self.IMPUT_LABEL = IMPUT_LABEL
        self.NA_LABEL = NA_LABEL
        self.NA_RATE = NA_RATE 
        self.INDEX_OLD_CHAR = INDEX_OLD_CHAR
        self.INDEX_NEW_CHAR = INDEX_NEW_CHAR
        self.COLUMNS_OLD_CHAR = COLUMNS_OLD_CHAR
        self.COLUMNS_NEW_CHAR = COLUMNS_NEW_CHAR
        
        
    def data_col_adapt(self, data, characters = [' ', '-'], new_char = '.'):
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

    def data_index_adapt(self, data, characters = ['-'], new_char = '.'):
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

    def data_col_index_adapt(self, data, col_characters = [' ', '-'], index_characters = ['-'], col_new_char = '.',
                            index_new_char = '.'):
        """
        由于后续需要用到R语言工具，对data行列名中不合法的字符进行替换。
        :param data: protein abundance
        :param col_characters
        :param index_characters
        :param new_char: new character
        :return data: renamed data
        """
        data_tmp1 = self.data_col_adapt(data=data, characters = col_characters, new_char = col_new_char)
        data_tmp2 = self.data_index_adapt(data=data_tmp1, characters = index_characters, new_char = index_new_char)
        return data_tmp2

    def data_index_recover(self, data, index_characters = ['-'], new_char = '.'):
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

    def col_adapt(self, samples_col, characters = [' ', '-'], new_char = '.'):
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
    
    def adapt_specific_protein(self, specific_protein, characters = ['-'], new_char = '.'):
        """
        由于后续需要用到R语言工具，对specific_protein中不合法的字符进行替换。
        :param specific_protein
        :param character: characters that need to be replaced
        :param new_char: new character
        :return specific_protein: renamed protein
        """
        for ch in characters:
            for idx in range(specific_protein.shape[0]):
                for col in range(specific_protein.shape[1]):
                    specific_protein.iloc[idx, col] = specific_protein.iloc[idx, col].replace(ch, new_char)

        return specific_protein

        
    def data_preprocess(self, preprocess_workdir, r_bpca_path):
        """
        Data preprocess
        主要包含这几部分：控制缺失值比例，对数化，z-score normalization, 调整列名与基因名，缺失值填充
        """
        
        preprocess_workdir = os.path.abspath(preprocess_workdir)
        if os.path.exists(preprocess_workdir) == False:
            os.makedirs(preprocess_workdir)
        os.chdir(preprocess_workdir)
        
        if self.HAS_SPECIFIC_PROTEIN:
            self.specific_protein_adapt = self.adapt_specific_protein(specific_protein = self.specific_protein,
                                                                      characters = self.INDEX_OLD_CHAR, 
                                                                      new_char = self.INDEX_NEW_CHAR)
            self.specific_protein_adapt.to_csv('./specific_protein_adapt.csv', index=False)
        
        if not self.HAS_PAIRED_DATA:

            # 将代表NA的标记修改为np.nan
            if self.IMPUT_LABEL:
                self.data = self.data.replace(self.NA_LABEL, np.nan)
            else:
                pass
            
            data_tumor = self.data.loc[:, self.tumor_cohort['tumor']]
            data_normal = self.data.loc[:, self.normal_cohort['normal']]
            # filter NA
            data_tumor_fna = data_tumor[data_tumor.isna().sum(axis=1) <= self.NA_RATE * data_tumor.shape[1]]
            data_normal_fna = data_normal[data_normal.isna().sum(axis=1) <= self.NA_RATE * data_normal.shape[1]]
            data_fna_index = set(data_tumor_fna.index).intersection(set(data_normal_fna.index))
            data_fna = self.data.loc[data_fna_index,:]

        
            # 2. 对数化
            if self.LOG_LABEL:
                data_fna_lt = copy.deepcopy(data_fna)
                for col in data_fna.columns:
                    data_fna_lt[col] = np.log2(data_fna[col].values)
            else:
                data_fna_lt = copy.deepcopy(data_fna)

            # 3. z-score normalization
            if self.NORMALIZATION:
                if self.NORMALIZATION == 'z-score':
                    data_fna_lt_zn = pd.DataFrame(stats.zscore(data_fna_lt, axis=0, nan_policy='omit'), 
                                                  index=data_fna_lt.index, columns=data_fna_lt.columns)
            else:
                data_fna_lt_zn = copy.deepcopy(data_fna_lt)

            # 4. Adapt columns and index 
            data_fna_lt_zn_da = self.data_col_index_adapt(data = data_fna_lt_zn, col_characters = self.COLUMNS_OLD_CHAR,
                                                         index_characters = self.INDEX_OLD_CHAR, col_new_char = self.COLUMNS_NEW_CHAR,
                                                         index_new_char = self.INDEX_NEW_CHAR)


            self.normal_cohort_adapt = self.col_adapt(samples_col = self.normal_cohort, 
                                                      characters = self.COLUMNS_OLD_CHAR, 
                                                      new_char = self.COLUMNS_NEW_CHAR)
            self.tumor_cohort_adapt = self.col_adapt(samples_col = self.tumor_cohort,
                                                    characters = self.COLUMNS_OLD_CHAR,
                                                    new_char = self.COLUMNS_NEW_CHAR)
            
            # 5. Imputation   
            self.tumor_cohort_adapt.to_csv('./tumor_cohort.txt', sep='\t', index=False)
            self.normal_cohort_adapt.to_csv('./normal_cohort.txt', sep='\t', index=False)
            data_fna_lt_zn_da.to_csv('./data.csv')

            utils.imputation_bpca_r(data_path = './data.csv', normal_cohort_path = './normal_cohort.txt', 
                                   tumor_cohort_path = './tumor_cohort.txt', r_bpca_path = r_bpca_path, 
                                   imput_result_path = './data_imput.csv')
            
            self.data_imput = pd.read_csv('./data_imput.csv', index_col=0)
            self.data_tumor_imput = self.data_imput.loc[:, self.tumor_cohort_adapt['tumor']]
            self.data_normal_imput = self.data_imput.loc[:, self.normal_cohort_adapt['normal']]
            
        else:
            self.paired_normal_cohort = pd.DataFrame({'normal':self.paired_samples.loc[:, 'normal']})
            self.paired_tumor_cohort = pd.DataFrame({'tumor':self.paired_samples.loc[:, 'tumor']})
            

            # 将代表NA的标记修改为np.nan
            if self.IMPUT_LABEL:
                data_fix_na = self.data.replace(self.NA_LABEL, np.nan)
                paired_data_fix_na = self.paired_data.replace(self.NA_LABEL, np.nan)

            else:
                data_fix_na = self.data
                paired_data_fix_na = self.paired_data

            
            data_tumor = data_fix_na.loc[:, self.tumor_cohort['tumor']]
            data_normal = data_fix_na.loc[:, self.normal_cohort['normal']]
            paired_data_tumor = paired_data_fix_na.loc[:, self.paired_tumor_cohort['tumor']]
            paired_data_normal = paired_data_fix_na.loc[:, self.paired_normal_cohort['normal']]
            
            # filter NA
            data_tumor_fna = data_tumor[data_tumor.isna().sum(axis=1) <= self.NA_RATE * data_tumor.shape[1]]
            data_normal_fna = data_normal[data_normal.isna().sum(axis=1) <= self.NA_RATE * data_normal.shape[1]]
            data_fna_index = set(data_tumor_fna.index).intersection(set(data_normal_fna.index))
            
            paired_data_tumor_fna = paired_data_tumor[paired_data_tumor.isna().sum(axis=1) <= self.NA_RATE * paired_data_tumor.shape[1]]
            paired_data_normal_fna = paired_data_normal[paired_data_normal.isna().sum(axis=1) <= self.NA_RATE * paired_data_normal.shape[1]]
            paired_data_fna_index = set(paired_data_tumor_fna.index).intersection(set(paired_data_normal_fna.index))
            
            _data_index = data_fna_index.intersection(paired_data_fna_index)
            
            data_fna = data_fix_na.loc[_data_index,:]
            paired_data_fna = paired_data_fix_na.loc[_data_index,:]

        
            # 2. 对数化
            if self.LOG_LABEL:
                data_fna_lt = copy.deepcopy(data_fna)
                for col in data_fna.columns:
                    data_fna_lt[col] = np.log2(data_fna[col].values)
                paired_data_fna_lt = copy.deepcopy(paired_data_fna)
                for col in paired_data_fna.columns:
                    paired_data_fna_lt[col] = np.log2(paired_data_fna[col].values)
            else:
                data_fna_lt = copy.deepcopy(data_fna)
                paired_data_fna_lt = copy.deepcopy(paired_data_fna)

            # 3. z-score normalization
            if self.NORMALIZATION:
                if self.NORMALIZATION == 'z-score':
                    data_fna_lt_zn = pd.DataFrame(stats.zscore(data_fna_lt, axis=0, nan_policy='omit'), 
                                                  index=data_fna_lt.index, 
                                                  columns=data_fna_lt.columns)
                    paired_data_fna_lt_zn = pd.DataFrame(stats.zscore(paired_data_fna_lt, axis=0, nan_policy='omit'),
                                                         index=paired_data_fna_lt.index,
                                                         columns=paired_data_fna_lt.columns)
            else:
                data_fna_lt_zn = copy.deepcopy(data_fna_lt)
                paired_data_fna_lt_zn = copy.deepcopy(paired_data_fna_lt)

            # 4. Adapt columns and index 
            data_fna_lt_zn_da = self.data_col_index_adapt(data = data_fna_lt_zn, 
                                                          col_characters = self.COLUMNS_OLD_CHAR,
                                                          index_characters = self.INDEX_OLD_CHAR, 
                                                          col_new_char = self.COLUMNS_NEW_CHAR,
                                                          index_new_char = self.INDEX_NEW_CHAR)

            paired_data_fna_lt_zn_da = self.data_col_index_adapt(data = paired_data_fna_lt_zn, 
                                                                 col_characters = self.COLUMNS_OLD_CHAR, 
                                                                 index_characters = self.INDEX_OLD_CHAR, 
                                                                 col_new_char = self.COLUMNS_NEW_CHAR,
                                                                 index_new_char = self.INDEX_NEW_CHAR)

            self.normal_cohort_adapt = self.col_adapt(samples_col = self.normal_cohort,
                                                      characters = self.COLUMNS_OLD_CHAR,
                                                      new_char = self.COLUMNS_NEW_CHAR)
            self.tumor_cohort_adapt = self.col_adapt(samples_col = self.tumor_cohort,
                                                     characters = self.COLUMNS_OLD_CHAR,
                                                     new_char = self.COLUMNS_NEW_CHAR)
            self.paired_normal_cohort_adapt = self.col_adapt(samples_col = self.paired_normal_cohort,
                                                             characters = self.COLUMNS_OLD_CHAR,
                                                             new_char = self.COLUMNS_NEW_CHAR)
            self.paired_tumor_cohort_adapt = self.col_adapt(samples_col = self.paired_tumor_cohort,
                                                            characters = self.COLUMNS_OLD_CHAR,
                                                            new_char = self.COLUMNS_NEW_CHAR)
            
            # 5. Imputation
            self.tumor_cohort_adapt.to_csv('./tumor_cohort.txt', sep='\t', index=False)
            self.normal_cohort_adapt.to_csv('./normal_cohort.txt', sep='\t', index=False)
            data_fna_lt_zn_da.to_csv('./data.csv')

            utils.imputation_bpca_r(data_path = './data.csv', 
                                    normal_cohort_path = './normal_cohort.txt', 
                                    tumor_cohort_path = './tumor_cohort.txt', 
                                    r_bpca_path = r_bpca_path, 
                                    imput_result_path = './data_imput.csv')
            
            self.data_imput = pd.read_csv('./data_imput.csv', index_col=0)
            self.data_tumor_imput = self.data_imput.loc[:, self.tumor_cohort_adapt['tumor']]
            self.data_normal_imput = self.data_imput.loc[:, self.normal_cohort_adapt['normal']]
            
   
            self.paired_tumor_cohort_adapt.to_csv('./paired_tumor_cohort.txt', sep='\t', index=False)
            self.paired_normal_cohort_adapt.to_csv('./paired_normal_cohort.txt', sep='\t', index=False)
            paired_data_fna_lt_zn_da.to_csv('./paired_data.csv')

            utils.imputation_bpca_r(data_path = './paired_data.csv', 
                                    normal_cohort_path = './paired_normal_cohort.txt', 
                                    tumor_cohort_path = './paired_tumor_cohort.txt', 
                                    r_bpca_path = r_bpca_path,
                                    imput_result_path = './paired_data_imput.csv')
            
            self.paired_data_imput = pd.read_csv('./paired_data_imput.csv', index_col=0)
            self.paired_data_tumor_imput = self.paired_data_imput.loc[:, self.paired_tumor_cohort_adapt['tumor']]
            self.paired_data_normal_imput = self.paired_data_imput.loc[:, self.paired_normal_cohort_adapt['normal']]