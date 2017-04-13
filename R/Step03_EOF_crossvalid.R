library(lubridate)
library(maptools)
library(floodmap)
library(magrittr)
library(zoo)
library(plyr)
library(data.table)
library(stringr)
library(ncdf4)
# rm(list = ls())

source('R/mainfunctions.R', encoding = 'UTF-8')
load("data/AVHRR_pheno_origin.rda")
file_nc <- "E:/GitHub/003Paper_TP_Phenology/data_004paper/AVHRR_NDVI_HantSmooth.nc"
# fid <- nc_open(file_nc)
# Id_grid <- ncvar_get(fid, "Id")
# Id <- fid$dim$Id$vals
pheno.df <- llply(pheno, function(x) as.data.frame(x))
method.group <- list(1:2, 4:5, 7:8, 11:14, 15:18, 19:22, 23:26, 27:30, 31:34) %>% 
  set_names(c("trs2", "trs5", "der", 
              paste(rep(c("zhang", "gu"), c(3, 3)), 
                    rep(c("spl", "beck", "elmr"), 2), sep = ".")))
x <- pheno.df[[1]]

pheno.trim <- list()
for (i in seq_along(pheno.df)){#1:300
  if (i %% 100 == 0) cat(sprintf("[%d]----------\n", i))
  pheno.trim[[i]] <- tryCatch(filter_rational(pheno.df[[i]]), 
                              error = function(e) message(sprintf("[%d]ERROR %s", i, e)))
}
# pheno.trim <- llply(pheno.df[1:100], filter_rational, .progress = "text")
names(pheno.trim) <- names(pheno.df)

Id <- 1:22211#just for GIMMS NDVI3g
Id_remain <- setdiff(Id, Id_del)
# Id_del <- Id_cal
save(pheno.trim, Id, Id_remain, Id_del, file = "AVHRR_pheno_trim.rda")

load("data/AVHRR_pheno_trim.rda")
## only EOS PHENOLOGICAL METRICS EXTRACT
Id_eos <- rep(list(2, 3:4), c(3, 6))
pheno.eos <- llply(pheno.trim, function(x) mapply(`[`, x, Id_eos) %>% do.call(cbind, .), 
                   .progress = "text")
## transfer the structure of data stored
varsEOS.origin <- colnames(pheno.eos[[1]])

phenoEOS.mat <- llply(seq_along(vars), function(i) ldply(pheno.eos, function(x) x[, i+1], .id = NULL), 
           .progress = "text") %>% set_names(vars)

nainfo <- ldply(phenoEOS.mat, function(x_df) apply(x_df, 1, function(x) length(which(is.na(x)))) %>% summary)
## ignore trs2 because too mach calculated failed
# trs6 in Step01__.R
phenoEOS <- c(phenoEOS.mat[1], list(trs6 = trs6.df[-Id_del, ]), phenoEOS.mat[2:14])
## -------------- GET SPATIAL DISTRIBUTION OF END OF SEASON PHENOLOGY
EOS.mean  <- llply(phenoEOS, function(x_df) apply(x_df, 1, mean, na.rm = TRUE), .progress = "text") %>% do.call(cbind.data.frame, .)
EOS.sd  <- llply(phenoEOS, function(x_df) apply(x_df, 1, sd, na.rm = TRUE), .progress = "text") %>% do.call(cbind.data.frame, .)

gridClip@data <- EOS.mean
spplot(gridClip, as.table = TRUE)

save(phenoEOS, EOS.sd, EOS.mean, file = "data/AVHRR_EOS_TP.rda")

# system.time({
#   for (i in seq_along(phenoEOS)){
#     fwrite(phenoEOS[[i]], file = sprintf("../AVHRR_%s.txt", vars[i]))
#   }
# })

