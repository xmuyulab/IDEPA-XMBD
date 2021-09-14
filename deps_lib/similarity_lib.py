import pandas as pd
import numpy as np
import copy
import os


def get_concordance_curve(qv_1, qv_2, N):

    concordance_score = []
    for col in qv_1.columns:
        _cs = []
        for n in range(1, N+1):
            _cs.append(len(set(qv_1[col].sort_values(ascending=True)[:n].index).intersection(set(qv_2[col].sort_values(ascending=True)[:n].index))))
        concordance_score.append(_cs)
    concordance_score = np.array(concordance_score)
    cs_mean = concordance_score.mean(axis=0)
    return cs_mean

def get_concordance_curve_ar(qv_1, qv_2, N):
    qv_1_mean = qv_1.mean(axis=1).sort_values(ascending=True)
    qv_2_mean = qv_2.mean(axis=1).sort_values(ascending=True)
    _cs = []
    for n in range(1, N+1):
        _cs.append(len(set(qv_1_mean[:n].index).intersection(set(qv_2_mean[:n].index))))
    return _cs

def get_concordance_score(y):
    x = range(len(y))
    area = np.trapz(y=y, x=x)
    score = area / (len(y)**2 / 2)
    return score


def get_algorithm_similarity(result_dir, workdir, N):
    """
    get algorithm similarity result
    :param result_dir 
    :param workdir
    :param N
    """
    if os.path.exists(workdir) == False:
        os.makedirs(workdir)
    os.chdir(workdir)

    penda_pro_up_qv_path = result_dir + '/penda_pro/penda_up_pro_qv.csv'
    penda_pro_down_qv_path = result_dir + '/penda_pro/penda_down_pro_qv.csv'

    rankc_v1_up_qv_path = result_dir + '/rankcomp/rankc_v1_up_qvalues.csv'
    rankc_v1_down_qv_path = result_dir + '/rankcomp/rankc_v1_down_qvalues.csv'

    rankc_v2_up_qv_path = result_dir + '/rankcomp/rankc_v2_up_qvalues.csv'
    rankc_v2_down_qv_path = result_dir + '/rankcomp/rankc_v2_down_qvalues.csv'

    ttest_up_path = result_dir + '/ttest/ttest_up.csv'
    ttest_down_path = result_dir + '/ttest/ttest_down.csv'
    ttest_qvalues_path = result_dir + '/ttest/ttest_qvalues.csv'

    wilcox_up_path = result_dir + '/wilcoxon/wilcox_up.csv'
    wilcox_down_path = result_dir + '/wilcoxon/wilcox_down.csv'
    wilcox_qvalues_path = result_dir + '/wilcoxon/wilcox_qvalues.csv'

    penda_pro_up_qv = pd.read_csv(penda_pro_up_qv_path, index_col=0)
    penda_pro_down_qv = pd.read_csv(penda_pro_down_qv_path, index_col=0)

    rankc_v1_up_qv = pd.read_csv(rankc_v1_up_qv_path, index_col=0)
    rankc_v1_down_qv = pd.read_csv(rankc_v1_down_qv_path, index_col=0)

    rankc_v2_up_qv = pd.read_csv(rankc_v2_up_qv_path, index_col=0)
    rankc_v2_down_qv = pd.read_csv(rankc_v2_down_qv_path, index_col=0)

    ttest_up = pd.read_csv(ttest_up_path, index_col=0)
    ttest_down = pd.read_csv(ttest_down_path, index_col=0)
    ttest_qvalues = pd.read_csv(ttest_qvalues_path, index_col=0)

    wilcox_up = pd.read_csv(wilcox_up_path, index_col=0)
    wilcox_down = pd.read_csv(wilcox_down_path, index_col=0)
    wilcox_qvalues = pd.read_csv(wilcox_qvalues_path, index_col=0)

    # 得到ttest_up 与 ttest_down 的 qvalues
    ttest_up_qv = copy.deepcopy(ttest_up)
    ttest_down_qv = copy.deepcopy(ttest_down)

    ttest_up_qv.iloc[:,:] = 1
    ttest_down_qv.iloc[:,:] = 1

    for idx in ttest_up.index:
        if ttest_up.loc[idx, :].sum() == ttest_up.shape[1]:
            ttest_up_qv.loc[idx, :] = ttest_qvalues.loc[idx, :]
        elif ttest_down.loc[idx, :].sum() == ttest_down.shape[1]:
            ttest_down_qv.loc[idx, :] = ttest_qvalues.loc[idx, :]

    wilcox_up_qv = copy.deepcopy(wilcox_up)
    wilcox_down_qv = copy.deepcopy(wilcox_down)

    wilcox_up_qv.iloc[:,:] = 1
    wilcox_down_qv.iloc[:,:] = 1

    for idx in wilcox_up.index:
        if wilcox_up.loc[idx, :].sum() == wilcox_up.shape[1]:
            wilcox_up_qv.loc[idx, :] = wilcox_qvalues.loc[idx, :]
        elif wilcox_down.loc[idx, :].sum() == wilcox_down.shape[1]:
            wilcox_down_qv.loc[idx, :] = wilcox_qvalues.loc[idx, :]

    up_v1_v2_cc = get_concordance_curve(rankc_v1_up_qv, rankc_v2_up_qv, N)
    up_v1_pp_cc = get_concordance_curve(rankc_v1_up_qv, penda_pro_up_qv, N)
    up_v1_tt_cc = get_concordance_curve(rankc_v1_up_qv, ttest_up_qv, N)
    up_v1_wil_cc = get_concordance_curve(rankc_v1_up_qv, wilcox_up_qv, N)

    up_v2_pp_cc = get_concordance_curve(rankc_v2_up_qv, penda_pro_up_qv, N)
    up_v2_tt_cc = get_concordance_curve(rankc_v2_up_qv, ttest_up_qv, N)
    up_v2_wil_cc = get_concordance_curve(rankc_v2_up_qv, wilcox_up_qv, N)

    up_pp_tt_cc = get_concordance_curve(penda_pro_up_qv, ttest_up_qv, N)
    up_pp_wil_cc = get_concordance_curve(penda_pro_up_qv, wilcox_up_qv, N)
    up_tt_wil_cc = get_concordance_curve(ttest_up_qv, wilcox_up_qv, N)


    down_v1_v2_cc = get_concordance_curve(rankc_v1_down_qv, rankc_v2_down_qv, N)
    down_v1_pp_cc = get_concordance_curve(rankc_v1_down_qv, penda_pro_down_qv, N)
    down_v1_tt_cc = get_concordance_curve(rankc_v1_down_qv, ttest_down_qv, N)
    down_v1_wil_cc = get_concordance_curve(rankc_v1_down_qv, wilcox_down_qv, N)

    down_v2_pp_cc = get_concordance_curve(rankc_v2_down_qv, penda_pro_down_qv, N)
    down_v2_tt_cc = get_concordance_curve(rankc_v2_down_qv, ttest_down_qv, N)
    down_v2_wil_cc = get_concordance_curve(rankc_v2_down_qv, wilcox_down_qv, N)

    down_pp_tt_cc = get_concordance_curve(penda_pro_down_qv, ttest_down_qv, N)
    down_pp_wil_cc = get_concordance_curve(penda_pro_down_qv, wilcox_down_qv, N)
    down_tt_wil_cc = get_concordance_curve(ttest_down_qv, wilcox_down_qv, N)

    _x = range(len(up_tt_wil_cc))

    import matplotlib.pyplot as plt

    plt.figure(figsize=(12,8))
    plt.plot(up_v1_v2_cc , linewidth=4, label='up_v1_v2: %.4f'%(get_concordance_score(up_v1_v2_cc)))
    plt.plot(up_v1_pp_cc , linewidth=4, label='up_v1_pp: %.4f'%(get_concordance_score(up_v1_pp_cc)))
    plt.plot(up_v1_tt_cc , linewidth=4, label='up_v1_tt: %.4f'%(get_concordance_score(up_v1_tt_cc)))
    plt.plot(up_v1_wil_cc , linewidth=4, label='up_v1_wil: %.4f'%(get_concordance_score(up_v1_wil_cc)))
    plt.plot(up_v2_pp_cc , linewidth=4, label='up_v2_pp: %.4f'%(get_concordance_score(up_v2_pp_cc)))
    plt.plot(up_v2_tt_cc , linewidth=4, label='up_v2_tt: %.4f'%(get_concordance_score(up_v2_tt_cc)))
    plt.plot(up_v2_wil_cc , linewidth=4, label='up_v2_wil: %.4f'%(get_concordance_score(up_v2_wil_cc)))
    plt.plot(up_pp_tt_cc , linewidth=4, label='up_pp_tt: %.4f'%(get_concordance_score(up_pp_tt_cc)))
    plt.plot(up_pp_wil_cc , linewidth=4, label='up_pp_wil: %.4f'%(get_concordance_score(up_pp_wil_cc)))
    plt.plot(up_tt_wil_cc , linewidth=4, label='up_tt_wil: %.4f'%(get_concordance_score(up_tt_wil_cc)))

    plt.legend(fontsize=15)
    plt.xlabel("Top N protein", fontsize=15)
    plt.ylabel("Concordance", fontsize=15)
    plt.title("Concordance curve", fontsize=20)
    plt.savefig(workdir+'/similarity_up.pdf', dpi=800, format='pdf')


    plt.figure(figsize=(12,8))
    plt.plot(down_v1_v2_cc , linewidth=4, label='down_v1_v2: %.4f'%(get_concordance_score(down_v1_v2_cc)))
    plt.plot(down_v1_pp_cc , linewidth=4, label='down_v1_pp: %.4f'%(get_concordance_score(down_v1_pp_cc)))
    plt.plot(down_v1_tt_cc , linewidth=4, label='down_v1_tt: %.4f'%(get_concordance_score(down_v1_tt_cc)))
    plt.plot(down_v1_wil_cc , linewidth=4, label='down_v1_wil: %.4f'%(get_concordance_score(down_v1_wil_cc)))
    plt.plot(down_v2_pp_cc , linewidth=4, label='down_v2_pp: %.4f'%(get_concordance_score(down_v2_pp_cc)))
    plt.plot(down_v2_tt_cc , linewidth=4, label='down_v2_tt: %.4f'%(get_concordance_score(down_v2_tt_cc)))
    plt.plot(down_v2_wil_cc , linewidth=4, label='down_v2_wil: %.4f'%(get_concordance_score(down_v2_wil_cc)))
    plt.plot(down_pp_tt_cc , linewidth=4, label='down_pp_tt: %.4f'%(get_concordance_score(down_pp_tt_cc)))
    plt.plot(down_pp_wil_cc , linewidth=4, label='down_pp_wil: %.4f'%(get_concordance_score(down_pp_wil_cc)))
    plt.plot(down_tt_wil_cc , linewidth=4, label='down_tt_wil: %.4f'%(get_concordance_score(down_tt_wil_cc)))

    plt.legend(fontsize=15)
    plt.xlabel("Top N protein", fontsize=15)
    plt.ylabel("Concordance", fontsize=15)
    plt.title("Concordance curve", fontsize=20)
    plt.savefig(workdir+'/similarity_down.pdf', dpi=800, format='pdf')