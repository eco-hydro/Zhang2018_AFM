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

Id <- seq_along(NDVI_list); #Id <- sample(Id, 10)

outdir <- "E:/GitHub/Vegetation/TP_phenology/TRS/"
xlist <- NDVI_list[Id]
RESULT <- llply(xlist, getYears_phenology, dates_avhrr.doy, outdir, PhenoFUN=get_Metric_TRS, .progress = "text")

## 修改为调试模式
result <- list()
for (i in seq_along(Id)){
  fprintf("[%d]---\n", i)
  
  ## 标准化之前
  xx <- melt(xlist[[i]]$xx) %>% set_names(c("year", "x", "value"))
  ggplot(xx, aes(x, value, year, color = year, fill = year, group = year)) + geom_line()
  
  ## 标准化之后
  x <- apply(xlist[[i]]$xx, 1, function(x) {
    minval <- min(x, na.rm = T)
    maxval <- max(x, na.rm = T)
    (x - minval)/(maxval - minval)
  }) %>% t
  xx <- melt(x) %>% set_names(c("year", "x", "value"))
  ggplot(xx, aes(x, value, year, color = year, fill = year, group = year)) + geom_line()
  
  result[[i]] <- getYears_phenology(xlist[[i]], dates_avhrr.doy, outdir, PhenoFUN=get_Metric_TRS)
}


ggplot(a, aes(x, value, year, fill = year, group = year)) + geom_line()
## 读取数据




