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

TEST_RATE = 0.3
num_test = floor(dim(normal_data)[2] * TEST_RATE)

normal_test = normal_data[, (1:num_test)]
normal_run = normal_data[, -(1:num_test)]
tumor_test = tumor_data[, (1:num_test)]
tumor_run = tumor_data[, -(1:num_test)]


simulation = list(initial_data = normal_test, simulated_data = tumor_test)
colnames(simulation$initial_data) = colnames(simulation$simulated_data)


threshold_LH = 0.99
s_max = 30
L_H_list = penda::compute_lower_and_higher_lists(normal_run, threshold = threshold_LH,
                                                 s_max = s_max)

L = L_H_list$L
H = L_H_list$H


threshold_values = c(0, 0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95, 0.99)

iterations = 30
quant_test = 0.05
factor_test = 1.2

which_threshold = penda::choose_threshold(normal_run, L_H_list = L_H_list,
                                          iterations = iterations,
                                          simulation = simulation,
                                          threshold_values = threshold_values,
                                          quant_test = quant_test,
                                          factor_test = factor_test)

thres = data.frame(which_threshold)
thres = data.frame(thres['patient'],thres['threshold'], thres['FDR'])
thres <- dplyr::mutate_all(thres,as.character)
thres <- dplyr::mutate_all(thres,as.numeric)

thres_temp = thres[which(thres$FDR < args[8]),]

thresh_min = 0
for(i in unique(c(thres_temp$patient))){
  temp = min(thres_temp[which(thres_temp$patient == i),]['threshold'])
  if(thresh_min < temp){
    thresh_min = temp
  }
}


penda_res = penda::penda_test(samples = tumor_data,
                              controls = normal_data,
                              threshold = thresh_min,
                              iterations = iterations,
                              L_H_list = L_H_list,
                              quant_test = quant_test,
                              factor_test = factor_test)



write.csv(penda_res$down_genes, args[9], sep=',', col.names = T, row.names = T)
write.csv(penda_res$up_genes, args[10], sep=',', col.names = T, row.names = T)
