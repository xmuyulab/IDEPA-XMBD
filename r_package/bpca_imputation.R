#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("pcaMethods")

library(pcaMethods)

options(warn =-1)

args <- commandArgs()

data_path = args[6]
sample_path = args[7]
output_path = args[8]

sample = read.csv(sample_path, header=TRUE, sep=',')
sample_11 = lapply(sample, function(x) as.character(x))
sample_1 = sample_11$tumor
sample_2 = sample_11$normal

data = read.csv(data_path, header = TRUE, sep=',')

data = data.matrix(data)

data_1 = data[,as.character(sample_1)]
data_2 = data[,as.character(sample_2)]

data_1_filter = scale(data_1)
data_2_filter = scale(data_2)

bpca_1 <- pca(data_1_filter, method="bpca", nPcs=5, scale = "none", center = FALSE)
bpca_2 <- pca(data_2_filter, method="bpca", nPcs=5, scale = "none", center = FALSE)

cObs_1 <- completeObs(bpca_1)
cObs_2 <- completeObs(bpca_2)

result = cbind(cObs_1, cObs_2)

write.csv(result, output_path, sep=',', col.names = TRUE)
