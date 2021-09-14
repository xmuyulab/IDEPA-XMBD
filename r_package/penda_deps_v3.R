#!/usr/bin/Rscript

library(penda)
library(data.table)
require(dplyr)

options(warn =-1)
args <- commandArgs()

normal_data_path = args[6]

tumor_data_path = args[7]


normal_data = read.csv(normal_data_path, header=TRUE, sep=',')

tumor_data = read.csv(tumor_data_path, header=TRUE, sep=',')

row.names(normal_data) = normal_data[,1]
normal_data = normal_data[,-1]
row.names(tumor_data) = tumor_data[,1]
tumor_data = tumor_data[,-1]

normal_data = data.matrix(normal_data)
tumor_data = data.matrix(tumor_data)

threshold_dataset = 1.0

dataset_penda = penda::make_dataset(normal_data,
                                  tumor_data,
                                  detectlowvalue = FALSE,
                                  detectNA = FALSE,
                                  threshold = threshold_dataset)

normal_data = dataset_penda$data_ctrl
tumor_data = dataset_penda$data_case

threshold_LH = 0.99
s_max = 30
L_H_list = penda::compute_lower_and_higher_lists(normal_data, threshold = threshold_LH,
                                                 s_max = s_max)
L = L_H_list$L
H = L_H_list$H


threshold = 0.4
iterations = 30
quant_test = 0.05
factor_test = 1.2


penda_res = penda::penda_test(samples = tumor_data,
                              controls = normal_data,
                              threshold = threshold,
                              iterations = iterations,
                              L_H_list = L_H_list,
                              quant_test = quant_test,
                              factor_test = factor_test)

write.csv(penda_res$down_genes, args[8], sep=',', col.names = T, row.names = T)
write.csv(penda_res$up_genes, args[9], sep=',', col.names = T, row.names = T)
