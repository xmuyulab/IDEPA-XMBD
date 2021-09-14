#!/usr/bin/Rscript

library(penda)
library(data.table)
require(dplyr)

options(warn =-1)

args <- commandArgs()

normal_run_path = args[6]
print(normal_run_path)

normal_test_path = args[7]
print(args[7])
print(normal_test_path)
tumor_run_path = args[8]
tumor_test_path = args[9]

normal_run = read.csv(normal_run_path, header=TRUE, sep=',')
normal_test = read.csv(normal_test_path, header=TRUE, sep=',')
tumor_run = read.csv(tumor_run_path, header=TRUE, sep=',')
tumor_test = read.csv(tumor_test_path, header=TRUE, sep=',')

normal_run = data.matrix(normal_run)
normal_test = data.matrix(normal_test)
tumor_run = data.matrix(tumor_run)
tumor_test = data.matrix(tumor_test)

row.names(normal_run) = normal_run[,1]
normal_run = normal_run[,-1]
row.names(tumor_run) = tumor_run[,1]
tumor_run = tumor_run[,-1]

row.names(normal_test) = normal_test[,1]
normal_test = normal_test[,-1]
row.names(tumor_test) = tumor_test[,1]
tumor_test = tumor_test[,-1]

# penda 表达量不能有小于0的值
MIN = min(min(normal_test), min(normal_run), min(tumor_run), min(tumor_test))

normal_test = normal_test + abs(1.1 * MIN)
normal_run = normal_run + abs(1.1 * MIN)
tumor_run = tumor_run + abs(1.1 * MIN)
tumor_test = tumor_test + abs(1.1 * MIN)

simulation = list(initial_data = normal_test, simulated_data = tumor_test)
colnames(simulation$initial_data) = colnames(simulation$simulated_data)

#DEREGULARATED_PROP = 0.3
#simulation = penda::simplified_simulation(data = normal_test, 
#                                          proportion = DEREGULARATED_PROP, 
#                                          threshold = 10,
#                                          modifier = 4,
#                                          factor = 4)



threshold_LH = 0.99
s_max = 30
L_H_list = penda::compute_lower_and_higher_lists(normal_run, threshold = threshold_LH,
                                                 s_max = s_max)

L = L_H_list$L
H = L_H_list$H

threshold_values = c(0, 0.05, 0.1, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99)
# threshold = 0.05
iterations = 30
quant_test = 0.05
factor_test = 1.2

which_threshold = penda::choose_threshold(normal_run, L_H_list = L_H_list,
                                          iterations = iterations,
                                          simulation = simulation,
                                          threshold_values = threshold_values,
                                          quant_test = 0,
                                          factor_test = 1)

thres = data.frame(which_threshold)
thres = data.frame(thres['patient'],thres['threshold'], thres['FDR'])
thres <- dplyr::mutate_all(thres,as.character)
thres <- dplyr::mutate_all(thres,as.numeric)

thres_temp = thres[which(thres$FDR < args[12]),]

thresh_min = 0
for(i in unique(c(thres_temp$patient))){
  temp = min(thres_temp[which(thres_temp$patient == i),]['threshold'])
  if(thresh_min < temp){
    thresh_min = temp
  }
}


penda_res = penda::penda_test(samples = simulation$simulated_data,
                                   controls = simulation$initial_data,
                                   threshold = thresh_min,
                                   iterations = iterations,
                                   L_H_list = L_H_list,
                                   quant_test = quant_test,
                                   factor_test = factor_test)
  


write.csv(penda_res$down_genes, args[10], sep=',', col.names = T, row.names = F)
write.csv(penda_res$up_genes, args[11], sep=',', col.names = T, row.names = F)
