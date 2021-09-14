library(pcaMethods)

options(warn =-1)
args <- commandArgs()

data_path = args[6]
tumor_col_path = args[7]
normal_col_path = args[8]
output_path = args[9]

data = read.csv(data_path, header=TRUE, sep=',')
tumor_col = read.table(tumor_col_path, header=TRUE)
normal_col = read.table(normal_col_path, header=TRUE)

row.names(data) = data[,1]
data = data[,-1]
data = data.matrix(data)


data_t = data[,as.character(tumor_col$tumor)]
data_n = data[,as.character(normal_col$normal)]


bpca_t <- pca(data_t, method="bpca", nPcs=5, scale = "none", center = FALSE)
bpca_n <- pca(data_n, method="bpca", nPcs=5, scale = "none", center = FALSE)

cObs_t <- completeObs(bpca_t)
cObs_n <- completeObs(bpca_n)

result = cbind(cObs_t, cObs_n)

write.csv(result, output_path, sep=',', col.names = TRUE, row.names = TRUE)
