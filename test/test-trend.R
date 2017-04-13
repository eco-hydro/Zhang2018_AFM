rm(list = ls())
source("R/mainfunc/main_TREND.R", encoding = "utf-8")
load("E:/GitHub/Vegetation/TP_phenology/a1.rda")

library(microbenchmark)
library(compiler)

mkTrend_cmf <- cmpfun(mkTrend)

microbenchmark(mkTrend(a1[1:200]))
