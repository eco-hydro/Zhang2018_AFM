library(stringr)
library(ncdf4)
rm(list = ls())
source('R/mainfunc/main.R', encoding = 'UTF-8')
files <- dir("TRS/", full.names=T, pattern = "*.txt") %>% {.[order(as.numeric(str_extract(., "\\d{1,5}")))]}

TRS <- llply(files, fread, sep = "\t", header = T, .progress = "text")
x <- fread(files[1])

trsnames <- c("TRS2", "TRS5", "TRS6")
SOS <- llply(TRS, function(x) x[, c(2, 5, 8)], .progress = "text")
SOS.list <- llply(1:3, function(j) ldply(SOS, `[`, ,j = j, .progress = "text")[-Id_del, ]) %>% 
  set_names(paste0(trsnames, ".SOS"))

EOS <- llply(TRS, function(x) x[, c(3, 6, 9)], .progress = "text")
EOS.list <- llply(1:3, function(j) ldply(EOS, `[`, ,j = j, .progress = "text")[-Id_del, ]) %>% 
  set_names(paste0(trsnames, ".EOS"))
# save(EOS.list, SOS.list, file = "TRS_suppy20170105.rda")

## test if eos if [gt] than sos
#  perfect results. Unique(losn) is zero
# losn <- laply(TRS, function(x) {
#   los <- x[, c(4, 7, 10)]
#   length(which(los < 0))#return
# }, .progress = "text")#

load("../data/AVHRR_pheno_origin.rda")
load("data/TRS_suppy20170105.rda")
#1:2, 4:5; "trs2", "trs5", 
method.group <- list(7:8, 11:14, 15:18, 19:22, 23:26, 27:30, 31:34) %>% 
  set_names(c("der", paste(rep(c("zhang", "gu"), c(3, 3)), 
                           rep(c("spl", "beck", "elmr"), 2), sep = ".")))

## 对数据格式进行优化
library(foreach)
library(iterators)
library(parallel)
# progress <- create_progress_bar("text")
# progress$init(length(pheno))
# pheno.df <- foreach(x = pheno, i = icount()) %do% {
#   progress$step()
#   x$year <- 1:32
#   x$id <- i
#   x
#   # data.table(id = i, year = 1:32, x)
# }

pheno_year <- llply(1:32, function(year) ldply(pheno, function(x) x[year, 7:34], .progress = "text", .id = NULL))
## directly use data.table to manipulate
# pheno.df <- llply(pheno, function(x) as.data.frame(x[, 7:34]), .progress = "text")

# file_nc <- "../../003Paper_TP_Phenology/data_004paper/AVHRR_NDVI_HantSmooth.nc"
# fid <- nc_open(file_nc)
# Id_grid <- ncvar_get(fid, "Id")
# Id <- fid$dim$Id$vals
cl <- makeCluster(7, outfile = "log.txt")
tmp <-clusterEvalQ(cl, {
  library(magrittr)
  library(data.table)
  #close all warnings, due to optimx check initial start point, can ignore it
  options(warn =-1)
})
clusterExport(cl, c("method.group", "filter_rational.dt"), envir = globalenv())

pheno.df <- parLapplyLB(cl, pheno, function(x) as.data.frame(x))
pheno.trim <- parLapplyLB(cl, pheno.df, filter_rational)

vars_all <- colnames(pheno.trim[[1]])
pheno.list <- llply(seq_along(vars_all), function(i) ldply(pheno.trim, function(x) x[, i], .id = NULL), 
                    .progress = "text") %>% set_names(vars_all)
NAMES <- read.table("data/names.txt", header = F, stringsAsFactors = F)$V1
pheno_list <- c(SOS.list, EOS.list, pheno.list) %>% set_names(NAMES)
nainfo <- ldply(pheno_list, function(x_df) apply(x_df, 1, function(x) length(which(is.na(x)))) %>% summary, .progress = "text")

## 数据划分为两部分，SOS, EOS

SOS_list <- pheno_list[c(1:3, 7, seq( 9, 32, 4), seq(10, 32, 4)) %>% sort]
EOS_list <- pheno_list[c(4:6, 8, seq(11, 32, 4), seq(12, 32, 4)) %>% sort]

EOS.mean  <- llply(EOS_list, function(x_df) apply(x_df, 1, mean, na.rm = TRUE), .progress = "text") %>% do.call(cbind.data.frame, .)
EOS.sd  <- llply(EOS_list, function(x_df) apply(x_df, 1, sd, na.rm = TRUE), .progress = "text") %>% do.call(cbind.data.frame, .)

SOS.mean  <- llply(SOS_list, function(x_df) apply(x_df, 1, mean, na.rm = TRUE), .progress = "text") %>% do.call(cbind.data.frame, .)
SOS.sd  <- llply(SOS_list, function(x_df) apply(x_df, 1, sd, na.rm = TRUE), .progress = "text") %>% do.call(cbind.data.frame, .)

# save(SOS_list, SOS.mean, SOS.sd, EOS_list, EOS.mean, EOS.sd, file = "data/phenology_list20170105.rda")
save(SOS_list, EOS_list, file = "data/phenology_list20170105.rda")



## only EOS PHENOLOGICAL METRICS EXTRACT
# Id_eos <- rep(list(2, 3:4), c(3, 6))

# Id_sos <- rep(list(1, 1:2), c(3, 6))
# pheno.sos <- llply(pheno.trim, function(x) mapply(`[`, x, Id_sos) %>% do.call(cbind, .), 
#                    .progress = "text")
# ## transfer the structure of data stored
# # varsEOS.origin <- colnames(pheno.eos[[1]])
# # phenoEOS.mat <- llply(seq_along(vars), function(i) ldply(pheno.eos, function(x) x[, i+1], .id = NULL), 
# #                       .progress = "text") %>% set_names(vars)
# phenoEOS.mat <- llply(seq_along(vars_origin), function(i) ldply(pheno.sos, function(x) x[, i], .id = NULL), 
#                       .progress = "text") %>% set_names(vars_origin)

## ignore trs2 because too mach calculated failed
# trs6 in Step01__.R
# phenoEOS <- c(phenoEOS.mat[1], list(trs6 = trs6.df[-Id_del, ]), phenoEOS.mat[2:14])
# ## -------------- GET SPATIAL DISTRIBUTION OF END OF SEASON PHENOLOGY
# 
# gridClip@data <- EOS.mean
# spplot(gridClip, as.table = TRUE)
# save(phenoEOS, EOS.sd, EOS.mean, file = "data/AVHRR_EOS_TP.rda")

# system.time({
#   for (i in seq_along(phenoEOS)){
#     fwrite(phenoEOS[[i]], file = sprintf("../AVHRR_%s.txt", vars[i]))
#   }
# })