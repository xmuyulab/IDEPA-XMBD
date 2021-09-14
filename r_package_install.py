#!/usr/bin/
# -*- coding: utf-8 -*-

'''
Install r package

'''

import rpy2
from rpy2.robjects.packages import importr
import rpy2.robjects.packages as rpackages
import rpy2.robjects as ro
import os

# R vector of strings
from rpy2.robjects.vectors import StrVector



# R package names(CRAN)
cran_packages = ('devtools', 'remotes', 'ggplot2', 'ggpubr', 'survival', 'survminer', 'BiocManager', 
                 'SuperExactTest', 'mixtools', 'dplyr', 'data.table', 'RSQLite')

# R package names(BIO)
bio_packages = ('clusterProfiler', 'org.Hs.eg.db', 'pcaMethods')


# import R utils package
utils = rpackages.importr('utils')

if not rpackages.isinstalled('foreign'):
    # 安装老版本的foreign
    rpy2.robjects.r("""
        packageurl <- "https://cran.r-project.org/src/contrib/Archive/foreign/foreign_0.8-71.tar.gz"
        install.packages(packageurl, repos=NULL, type="source")
        """)


# select a mirror for R packages
utils.chooseCRANmirror(ind=1) # select the first mirror in the list

# Selectively install what needs to be install.
cran_to_install = [x for x in cran_packages if not rpackages.isinstalled(x)]
if len(cran_to_install):
    utils.install_packages(StrVector(cran_to_install))

if not rpackages.isinstalled('penda'):

    home_dir = os.getcwd()

    penda_dir = home_dir + '/deps_lib'
    os.chdir(penda_dir)

    os.system('git clone https://github.com/bcm-uga/penda.git')

    os.chdir(penda_dir + '/penda')

    rpy2.robjects.r('''
            # install.packages("devtools")
            # mixtools htmltools scales yaml lazyeval plyr rlang ggplot2 gtools caTools KernSmooth
            devtools::load_all(); devtools::document(); devtools::install()
            ''')
    
    
utils.chooseBioCmirror(ind=1) # select the first mirror in the list
biocinstaller = importr("BiocManager")
# Selectively install what needs to be install.
bio_to_install = [x for x in bio_packages if not rpackages.isinstalled(x)]
if len(bio_to_install):
    biocinstaller.install(StrVector(bio_to_install))
    







