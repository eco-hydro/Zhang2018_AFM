library(maptools)
library(ncdf4)
library(lubridate)
library(floodmap)
library(magrittr)
library(zoo)
library(plyr)
library(optimx)
rm(list = ls())
# source('R_004paper/PhenologyMetrics_new.R', encoding = 'UTF-8')
source('R/物候数据计算_avhrr/optimDoubleLostics.R', encoding = 'UTF-8')
source('R/物候数据计算_avhrr/PhenoExtract.R', encoding = 'UTF-8')
source('R/parallel/parallel_remote_functions.R', encoding = 'UTF-8')

Id_del <- c(16035, 18443, 20925, 21407, 21618, 21732, 21827)
outdir <- "E:/GitHub/Vegetation/TRS6/"
files_finish <- get_files(outdir, pattern="*.txt")
Id <- 1:22211
if (length(files_finish) > 0){
  Id_finish <- stringr::str_extract(files_finish, "\\d{1,}") %>% as.numeric()
  # Id_cal <- setdiff(Id, Id_finish)
}else{
  Id_cal <- Id
}

## ------------------- match functions -----------------------
ls.isfun <- function(x) eval(parse(text = sprintf("is.function(%s)", x)))
funcs <- c(ls()[sapply(ls(), ls.isfun)], ".error")
# OPTIMX METHODS
methods <- c('BFGS','CG','Nelder-Mead','L-BFGS-B','nlm','nlminb','spg','ucminf','Rcgmin','Rvmmin','newuoa','bobyqa','nmkb','hjkb')

file_ahvrr <- "../data/AVHRR_NDVI_HantSmooth.nc"
fid <- nc_open(file_ahvrr)
NDVI.avhrr <- ncvar_get(fid, "NDVI")
NDVI.avhrr %<>% {array(., dim = c(nrow(.), 24, 32))} %>% aperm(c(3, 2, 1))
# apply(NDVI.avhrr, 2, mean, na.rm = F) %T>% plot() 
NDVI_list <- alply(NDVI.avhrr, 3, function(x) x) 
for (i in seq_along(NDVI_list)) NDVI_list[[i]] <- list(xx=NDVI_list[[i]], i=i)
# get AVHRR NDVI components date information
date_up <- seq.Date(as.Date("1982-01-01"), as.Date("2013-12-31"), by = "month")
get_middays <- function(date) days_in_month(date) %>% {date + floor(quantile(seq(.), c(1/4, 3/4))) - 1}
dates_avhrr <- lapply(date_up, get_middays) %>% do.call(c, .)
dates_avhrr.doy <- as.numeric(format(dates_avhrr, "%j")) %>% matrix(ncol = 24, byrow = TRUE)

nyear <- 32; ngrid <- dim(NDVI.avhrr)[3]
xx <- NDVI_list[[1000]]

Id <- seq_along(NDVI_list)

# j <- 1
# t <- dates_avhrr.doy[j, ]
# tout <- t[1]:t[length(t)]
# xs <- zoo(xx[j, ], t)
# # xl <- spline(doy, xx[i, ], xout = tout)$y %>% zoo(., tout)
# a <- FitDoubleLogElmore(xs, t,tout)
# b <- FitDoubleLogBeck(xs, t, tout)

# a <- getYears_phenology(xx, dates_avhrr.doy)
options(warn =-1)#close all warnings, due to optimx check initial start point, can ignore it
# options(warn = 0)#open warnings

## ------------- THEN GO TO parallel_initial.R ---------------------------

# phenoogyMetrics <- result
# ##-------------------------------------
# get_Metrics <- function(data) sapply(data, function(x) unlist(x[1:4]))
# ##-------------------------------------
# file_metrics <- "data_out/TP_MetricsOnly.rda"
# if (!file.exists(file_metrics)){
#   clusterExport(cl, c("get_Metrics"))
#   ## 首先分析生长季的长度, 综合分析应该截去多个生长季的区域，生长季跨年的区域仍然保留
#   sm <- snow.time(metrics <- parLapply(cl, phenoogyMetrics, get_Metrics))
#   plot(sm)#show the time parallel use
#   save(metrics, file = file_metrics)
# }else{
#   load(file_metrics)
# }
# # save(phenoogyMetrics, file = "data_out/TP_phenoogyMetrics.rda")

# ## ---------------------- MULTI-COMPUTER PARALLEL ----------------------

# sm <- snow.time(metrics <- parLapply(cl, phenoogyMetrics, get_Metrics))
# plot(sm)#show the time parallel use

## 单独计算TRS6指标
library(parallel)
cl <- makeCluster(7, outfile = "log.txt")

tmp <-clusterEvalQ(cl, {
  library(magrittr)
  library(zoo)
  library(plyr)
  library(optimx)
  library(zoo)
  #close all warnings, due to optimx check initial start point, can ignore it
  options(warn =-1)
})

outdir <- "E:/GitHub/Vegetation/TP_phenology/TRS/"
clusterExport(cl, c("dates_avhrr.doy", "fprintf", funcs), envir = globalenv())
Id_cal <- seq_along(NDVI_list)
RESULT <- parLapplyLB(cl, NDVI_list[Id_cal], getYears_phenology, dates_avhrr.doy, outdir, PhenoFUN=get_Metric_TRS)

# llply(xlist, getYears_phenology, dates_avhrr.doy, outdir, fun=get_Metric_TRS)
# 
# for (i in seq_along(Id_cal)){
#   tmp <- getYears_phenology(NDVI_list[[Id_cal[i]]], dates_avhrr.doy, outdir, FUN=get_Metric_TRS)
# }
# 
# xlist = NDVI_list[1000:1200]
files <- get_files(outdir, full.names = T)
xlist <- llply(files[order(Id_finish)][-Id_del],
               fread, select = 5:6,
               .progress = "text")#, select = c(2, 3)
# trs6 <- ldply(xlist, `[[`, 6)
trs6 <- llply(xlist, function(x_df){
        # Id_na <- which((x_df[, 2] - x_df[, 1] ) <= 0)
        Id_na <- apply(x_df, 1, order) %>% apply(2, function(x) {
          x_unique <- unique(diff(x))
          (length(x_unique) != 1) || (x_unique != 1)#should be attention that `|` different from `||`
        }) %>% which
        x_df[Id_na, ] <- NA#if Id_na is NULL have no any negtive effects
        x_df#需要对trs2结果进行改进
      }, .progress = "text")
trs6.df <- ldply(trs6, `[[`, "eos")

load("data/00basemap.rda")

a <- llply(xlist, function(x) apply(x[, c(3, 6)], 2, mean), .progress = "text")
x <- do.call(rbind.data.frame, a) %>% set_colnames(c("trs5", "trs6"))
gridClip@data <- x
spplot(gridClip)