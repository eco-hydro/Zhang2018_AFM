rm(list = ls())
source('R/mainfunc/main.R', encoding = 'UTF-8')
source('R/mainfunc/main_PLS.R', encoding = 'UTF-8')
library(purrr)

# load("data/phenology_list20170105.rda")
# load("../TP_phenoAgr/agrmet25_TP.rda")
load("../data/00basemap.rda")

file = 'E:/1.论文代码/04论文-再分析数据/mete_month.mat'
dataIn <- R.matlab::readMat(file)

Id_del <- c(16035, 18443, 20925, 21407, 21618, 21732, 21827)
dataIn.trim <- lapply(dataIn, function(x) x[-Id_del, ])

x <- dataIn.trim$TminMonth


infos <- llply(dataIn.trim, function(x) aggregate(x, list(cluster = Id_cluster), mean, na.rm = T))
# infos %<>% llply(set_colnames)

infos %<>% map(~ set_colnames(.x, c("clustet", paste0("month", 1:12))))
writelist_ToXlsx(infos, "mete_Avg_monthly.xlsx")
