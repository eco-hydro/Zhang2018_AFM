library(ncdf4)
library(optimx)
library(ppcor)
rm(list = ls())

source('R/mainfunc/main.R', encoding = 'UTF-8')
source('R/mainfunc/main_VSI.R', encoding = 'UTF-8')
load("../data/00basemap.rda")

## 测试是否进行Hant-smooth处理对结果造成的差异
## ------------------- match functions -----------------------
# ls.isfun <- function(x) eval(parse(text = sprintf("is.function(%s)", x)))
# funcs <- c(ls()[sapply(ls(), ls.isfun)], ".error")

file_mete <- "data/Mete_analysis_dataIn_origin.rda"
if (!file.exists(file_mete)){
  # file_ahvrr <- "../data/AVHRR_NDVI_HantSmooth.nc"
  file_ahvrr <- "../data/AVHRR_NDVI_TP_GRID.nc"
  fid <- nc_open(file_ahvrr)
  NDVI.avhrr <- ncvar_get(fid, "NDVI")
  NDVI.avhrr %<>% apply(1, function(x) matrix(x, nrow=2) %>% apply(2, max, na.rm = T)) %>% t
  
  metes <- R.matlab::readMat("../data/plsMeteInputData.mat") %>% llply(function(x) x[, 13:396])
  # varnames <- names(vars)
  
  xlist <- list()
  for (i in 1:nrow(NDVI.avhrr)){
    if (i %% 100 == 0) cat(sprintf("[%05dth]\n", i))
    NDVI <- NDVI.avhrr[i, ] %>% matrix(nrow = 12)
    NDVI_t1 <- matrix(c(NDVI.avhrr[i, -1], NA), nrow = 12)
    mete <- llply(metes, function(x) matrix(x[i, ], nrow = 12))
    
    dataIn <- c(list(NDVI = NDVI, t1 = NDVI_t1), mete)
    # xlist <- lapply(1:12, function(i) lapply(dataIn, `[`, i = i, )) %>% set_names(1:12)
    xlist[[i]] <- lapply(1:12, function(i) lapply(dataIn, `[`, i = i, ) %>% 
                           {na.omit(do.call(cbind.data.frame, .))}) %>% set_names(1:12)
  }
  xlist %<>% set_names(1:nrow(NDVI.avhrr)); xlist <- xlist[-Id_del]
  save(xlist, file = file_mete)
}else{
  load(file_mete)
}

library(parallel)
cl <- makeCluster(7, outfile = "log.txt")
tmp <- clusterEvalQ(cl, library(magrittr))
clusterExport(cl, c("pcor_month", "vsi_month", "listk"), envir = globalenv())

## vsi
# no t1
system.time(VSI <- parLapplyLB(cl, xlist, vsi_year))
save(VSI, file = "data/Vegeteation_Sensitivity_20170220_origin.rda")

# have t1
system.time(VSI <- parLapplyLB(cl, xlist, vsi_year, t1 = TRUE))
save(VSI, file = "data/Vegeteation_Sensitivity_20170220_origin.rda")

vsi_coef <- VSI_data(VSI)

brks <- c(0.03, 0.15, 10) %>% {c(-rev(.), 0, .)}

source('R/mainfunc/main_VSI.R', encoding = 'UTF-8')
VSI_show_coef(vsi_coef$coef[, 1:30], brks, filename = "VSI_coef_origin_Tavg.png")

VSI_show_dominate(vsi_coef$opt, filename = "VSI_dominate_origin_Tavg.png")

## debug error
for (i in 22004:22204){
  print(i)
  tmp <- vsi_year(xlist[[i]])
}

## pcor
system.time(PCOR_t1 <- parLapplyLB(cl, xlist, pcor_year, t1 = TRUE))
system.time(PCOR <- parLapplyLB(cl, xlist, pcor_year, t1 = FALSE))

meteVars <- c("t1", "Tmax", "Tavg", "Tmin", "Prec", "Srad")[-1]
Id <- match(meteVars, colnames(PCOR_t1[[1]]))

monthId <- 4:10 %>% {set_names(as.list(.), .)}
varCoef.list <- llply(monthId, function(i) ldply(PCOR_t1, function(x) x[i, Id], .id = NULL), .progress = "text")
x <- list.cbind(varCoef.list) %>% set_colnames(paste0("M", colnames(.)))

VSI_show_coef(x[, seq(1, length(meteVars)*6)], brks = c(Rc, 1), filename = "pcor_coef_origin-t1.png", sign = T)


# sign <- sign(x)
# "t1"   "Prec" "Srad" "Tmin" "Tavg" "Tmax"

save(PCOR, file = "Vegetation_PCOR.rda")
# pcor(na.omit(full))$estimate[1, 2:7]
