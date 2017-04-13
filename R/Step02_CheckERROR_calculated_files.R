library(lubridate)
library(floodmap)
library(magrittr)
library(zoo)
library(plyr)
library(data.table)
library(stringr)
rm(list = ls())

outdir <- "E:/Phenology/AVHRR/"
indir <- "E:/Phenology/tmp/AVHRR/"
indir2 <- "E:/Phenology/tmp/AVHRR2/"

files <- dir(indir, full.names = T, pattern = "*.txt")
# data <- llply(files, fread, header = T, sep = "\t", .progress = "text")
data <- list()
for (i in seq_along(files)){
  if (i %% 1000 == 0) cat(sprintf("[%d] %s\n", i, basename(files[i])))
  data[[i]] <- tryCatch(fread(files[i], header = T, sep="\t"), 
                        warning = function(w){
                          message(sprintf("[%d]WARNING %s %s", i, basename(files[i]), w))
                          suppressWarnings(fread(files[i], header = T, sep="\t"))#still return data
                        },
                        error = function(e) message(sprintf("[%d]ERROR %s %s", i, basename(files[i]), e)))
}
data %<>% set_names(basename(files))
# tmp <- data[1:1000]
# data %<>% llply(., function(x) {
#   class(x) <- "data.frame"
#   # x <- `[<-.data.frame`(x, x < 0, value = NA)
#   # x[x < 0] <- NA
#   rownames(x) <- 1982:2013
#   class(x) <- c("data.table", "data.frame")
#   return(x)
# }, .progress = "text")

Id_null <- which(sapply(data, is.null))
if (length(Id_null)) data <- data[-Id_null]
nyears <- sapply(data, nrow) %T>% {print(table(.))}

## --------- 第二次计算的数据 -----------
files <- dir(indir2, full.names = T, pattern = "*.txt")
data2 <- list()
for (i in seq_along(files)){
  if (i %% 200 == 0) cat(sprintf("[%d] %s\n", i, basename(files[i])))
  data2[[i]] <- tryCatch(fread(files[i], header = T, sep="\t"), 
                        warning = function(w){
                          message(sprintf("[%d]WARNING %s %s", i, basename(files[i]), w))
                          suppressWarnings(fread(files[i], header = T, sep="\t"))#still return data
                        },
                        error = function(e) message(sprintf("[%d]ERROR %s %s", i, basename(files[i]), e)))
}
data2 %<>% set_names(basename(files))
nyears <- sapply(data2, nrow) %T>% {print(table(.))}
Id_null <- which(sapply(data2, is.null) | nyears < 32)
if (length(Id_null)) data2 <- data2[-Id_null]
## delete year column var in order to keep same form
data2 <- lapply(data2, function(x) x[, year:= NULL])#delete year column
# for (i in seq_along(data2))

pheno <- c(data, data2)
Id_finish <- names(pheno) %>% str_extract("\\d{1,}") %>% as.numeric()
pheno %<>% set_names(Id_finish) %>% {.[order(Id_finish)]}

Id <- 1:22211#just for GIMMS NDVI3g
Id_cal <- setdiff(Id, Id_finish)
Id_del <- Id_cal

## SAVE AVHRR NDVI PHENOLOGY RESULT
save(pheno, Id_del, file = "AVHRR_pheno_origin.rda")