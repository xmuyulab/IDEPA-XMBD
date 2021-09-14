import pandas as pd
import os
import copy

def read_data(data_path):
    """
    从硬盘上读入蛋白质丰度数据，文件格式为.csv，第一列为蛋白质名。
    """
    data = pd.read_csv(data_path, sep=',', index_col=0)
    return data

def read_paired_data(data_path):
    """
    从硬盘上读入蛋白质丰度数据，文件格式为.csv，第一列为蛋白质名。
    """
    data = pd.read_csv(data_path, sep=',', index_col=0)
    return data

def read_paired_samples(paired_samples_path):
    """
    从硬盘上读入蛋白质样本信息，文件格式为.csv，文件包含两列信息：tumor队列样品名，normal队列样品名。
    """
    paired_samples = pd.read_csv(paired_samples_path, sep='\t')
    return paired_samples

def read_normal_cohort(normal_cohort_path):
    """
    读取normal cohort samples message
    """
    normal_cohort = pd.read_csv(normal_cohort_path, sep='\t')
    return normal_cohort

def read_tumor_cohort(tumor_cohort_path):
    """
    读取tumor cohort samples message
    """
    tumor_cohort = pd.read_csv(tumor_cohort_path, sep='\t')
    return tumor_cohort
    
def read_specific_protein(specific_protein_path):
    """
    读取specific protein message
    """
    specific_protein = pd.read_csv(specific_protein_path)
    return specific_protein

def imputation_bpca_r(data_path, normal_cohort_path, tumor_cohort_path, r_bpca_path, imput_result_path):
    """
    对data中不同队列数据分别填充
    """
    os.system('Rscript %s %s %s %s %s'%(r_bpca_path, data_path, tumor_cohort_path, normal_cohort_path, imput_result_path))
    
    
def write_normal_dat(normal_data, dir_path):
    """
    将normal_data以dat的格式写入硬盘
    :param normal_data: normal data
    :param dir_path: target directory
    """
    normal_data.to_csv(dir_path + '/normal.dat', sep='\t', index=0, header=0)
    
def write_tumor_dat(tumor_data, dir_path):
    """
    将tumor_data以dat的格式写入硬盘
    :param tumor_data: tumor data
    :param dir_path: target directory
    """
    tumor_data.to_csv(dir_path + '/tumor.dat', sep='\t', index=0, header=0)
    
    
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
    
    
def rankc_j1(reoa_path, normal_data, tumor_data, workdir, cycle_rankc, fdr, max_exception_rankc):
    """
    Rankcomp找normal 稳定对，逆转对，一致对
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
    os.system('%s -s 0 -j 1 -m %d -f %f -a 0 2 %d %s/normal.dat %s/tumor.dat %d %d %f %f'%
              (reoa_path, cycle_rankc, fdr, n_prot, workdir, workdir, n_n_sample, n_t_sample, 
               max_exception_rankc, max_exception_rankc))
    
def get_rankc_j1_results(workdir):
    """
    读取Rankcomp j1结果中的normal 稳定对，逆转对，一致对
    :param workdir: rankcomp j1 work space
    :return sp_concordant: concordant
    :return sp_reversed: reversed
    :return sp_normal: normal stable pairs
    """
    try:
        sp_concordant = pd.read_csv(workdir + '/concordant_pairs_1_0.dat',sep='\t', header=None)    
    except:
        print('Error: Concordant pair is empty!\n')

    try:
        sp_reversed = pd.read_csv(workdir + '/reversed_pairs_1_0.dat', sep='\t', header=None)
    except:
        print('Error: Reversal pair is empty!\n')
        
    try:
        sp_normal = pd.read_csv(workdir + '/stable_pairs_0.dat', sep='\t', header=None)
    except:
        print('Error: Normal stable pairs is empty!\n')
        
    return sp_concordant, sp_reversed, sp_normal

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
    os.system('%s -s 1 -j 2 -m %d -f %f -a 0 2 %d %s/normal.dat %s/tumor.dat %d %d %f %f'%
              (reoa_path, cycle_rankc, fdr, n_prot, workdir, workdir, n_n_sample, n_t_sample, 
               max_exception_rankc, max_exception_rankc))
    
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
    os.system('%s -s 1 -j 2 -m %d -f %f -a 2 2 %d %s/normal.dat %s/tumor.dat %d %d %f %f'%
              (reoa_path, cycle_rankc, fdr, n_prot, workdir, workdir, n_n_sample, n_t_sample, 
               max_exception_rankc, max_exception_rankc))
    
def rankc_j2_qvalues(reoa_path, normal_data, tumor_data, workdir, cycle_rankc, fdr, max_exception_rankc):
    """
    Rankcomp找差异蛋白，并且包含q values信息
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
    os.system('%s -s 1 -j 2 -m %d -f %f -a 2 -v 2 %d %s/normal.dat %s/tumor.dat %d %d %f %f'%
              (reoa_path, cycle_rankc, fdr, n_prot, workdir, workdir, n_n_sample, n_t_sample, 
               max_exception_rankc, max_exception_rankc))

    
    

def get_rankc_j2_qvalues_result(workdir, data_run_tumor):
    """
    读取Rankcomp j2 包含q values的结果
    :param work
    :param data_run_tumor
    """
    
    os.chdir(workdir)
    
    rankc_v1_up = copy.deepcopy(data_run_tumor)
    rankc_v1_up.iloc[:,:] = False
    rankc_v1_down = copy.deepcopy(data_run_tumor)
    rankc_v1_down.iloc[:, :] = False
    rankc_v1_up_qvalues = copy.deepcopy(data_run_tumor)
    rankc_v1_up_qvalues.iloc[:, :] = 1
    rankc_v1_down_qvalues = copy.deepcopy(data_run_tumor)
    rankc_v1_down_qvalues.iloc[:, :] = 1

    rankc_v2_up = copy.deepcopy(data_run_tumor)
    rankc_v2_up.iloc[:,:] = False
    rankc_v2_down = copy.deepcopy(data_run_tumor)
    rankc_v2_down.iloc[:, :] = False
    rankc_v2_up_qvalues = copy.deepcopy(data_run_tumor)
    rankc_v2_up_qvalues.iloc[:, :] = 1
    rankc_v2_down_qvalues = copy.deepcopy(data_run_tumor)
    rankc_v2_down_qvalues.iloc[:, :] = 1

    for idx in range(data_run_tumor.shape[1]):
        _down_v1_name = 'down_regulated_1c%sv0.dat'%idx
        _up_v1_name = 'up_regulated_1c%sv0.dat'%idx
        if os.path.getsize(_down_v1_name) != 0:
            _down_v1 = pd.read_csv(_down_v1_name, sep='\t', header=None)
            rankc_v1_down.iloc[_down_v1[0].tolist(), idx] = True
            rankc_v1_down_qvalues.iloc[_down_v1[0].tolist(), idx] = _down_v1[1].tolist()
        if os.path.getsize(_up_v1_name) != 0:
            _up_v1 = pd.read_csv(_up_v1_name, sep='\t', header=None)
            rankc_v1_up.iloc[_up_v1[0].tolist(), idx] = True
            rankc_v1_up_qvalues.iloc[_up_v1[0].tolist(), idx] = _up_v1[1].tolist()
        
        _down_v2_name = 'down_regulated_1c%s_0.dat'%idx
        _up_v2_name = 'up_regulated_1c%s_0.dat'%idx
        if os.path.getsize(_down_v2_name) != 0:
            _down_v2 = pd.read_csv(_down_v2_name, sep='\t', header=None, engine='python')
            rankc_v2_down.iloc[_down_v2[0].tolist(), idx] = True
            rankc_v2_down_qvalues.iloc[_down_v2[0].tolist(), idx] = _down_v2[1].tolist()
        if os.path.getsize(_up_v2_name) != 0:
            _up_v2 = pd.read_csv(_up_v2_name, sep='\t', header=None, engine='python')
            rankc_v2_up.iloc[_up_v2[0].tolist(), idx] = True
            rankc_v2_up_qvalues.iloc[_up_v2[0].tolist(), idx] = _up_v2[1].tolist()
            

    rankc_v1_up.to_csv('./rankc_v1_up.csv')
    rankc_v1_down.to_csv('./rankc_v1_down.csv')
    rankc_v1_up_qvalues.to_csv('./rankc_v1_up_qvalues.csv')
    rankc_v1_down_qvalues.to_csv('./rankc_v1_down_qvalues.csv')

    rankc_v2_up.to_csv('./rankc_v2_up.csv')
    rankc_v2_down.to_csv('./rankc_v2_down.csv')
    rankc_v2_up_qvalues.to_csv('./rankc_v2_up_qvalues.csv')
    rankc_v2_down_qvalues.to_csv('./rankc_v2_down_qvalues.csv')

    return rankc_v1_up, rankc_v1_down, rankc_v2_up, rankc_v2_down

    
def get_rankc_j2_results(workdir):
    """
    读取Rankcomp j2结果中的差异蛋白
    :param workdir: rankcomp j2 work space
    :return dep: difference expression result
    """
    try:
        dep = pd.read_csv(workdir + '/gene_state_1.dat',sep='\t', header=None)    
    except:
        print('Error: DEP is empty!\n')
        
    return dep

def write_penda_data(normal_run, tumor_run, dir_path):
    """
    将penda所需数据写入硬盘
    :param normal_run
    :param tumor_run
    """
    normal_run.to_csv(dir_path + '/normal_run.csv')
    tumor_run.to_csv(dir_path + '/tumor_run.csv')
    

def run_penda_fdr(r_penda_path, fdr):
    
    os.system('Rscript %s ./normal_run.csv ./tumor_run.csv %s ./penda_fdr_down.csv ./penda_fdr_up.csv'%(r_penda_path, fdr))
    
def get_penda_fdr_result(workdir):
 
    _penda_up = pd.read_csv(workdir + '/penda_fdr_up.csv', index_col=0)
    _penda_down = pd.read_csv(workdir + '/penda_fdr_down.csv', index_col=0)
    _penda_result = _penda_up | _penda_down
    return _penda_up, _penda_down, _penda_result

def run_penda(r_penda_path):
    
    os.system('Rscript %s ./normal_run.csv ./tumor_run.csv ./penda_down.csv ./penda_up.csv'%(r_penda_path))
    
def get_penda_result(workdir):
  
    _penda_up = pd.read_csv(workdir + '/penda_up.csv', index_col=0)
    _penda_down = pd.read_csv(workdir + '/penda_down.csv', index_col=0)
    _penda_result = _penda_up | _penda_down
    return _penda_up, _penda_down, _penda_result


    
