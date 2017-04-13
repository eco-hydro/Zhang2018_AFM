# library(snow)
library(parallel)
library(magrittr)
# library(stringr)#stringr package need to be installed
source('R/parallel/parallel_remote_functions.R', encoding = "utf-8")
source('R/parallel/parallel_package.R', encoding = "utf-8")
## ------------------- debug model -----------------------
# debugSource('R/snowSOCK.R')
# source('R/snow.R', encoding = "utf-8")
# source('R/RngStream.R', encoding = "utf-8")
# source('R/detectCores.R', encoding = "utf-8")
# 
# libpath <- dirname(system.file(package = "parallel"))
# initDefaultClusterOptions(libpath)

#' @param host remote computer host name or ip
#' @param user remote computer login username
#' @param rscript remote computer R software bin path
#' @param rshcmd plink login command with password
localhost <- list(host = "localhost", rshcmd = "plink")
zhang <- list(
  host = "zhang",
  user = "kongdd",
  rscript = "C:/Program Files/R/R-3.3.2/bin/Rscript.exe",
  rshcmd = "plink -pw xiuyuan156"
)
kongdd <- list(
  host = "kongdd-pc",
  user = "kongdd",
  rscript = "C:/Program Files/R/R-3.3.2/bin/Rscript.exe",
  rshcmd = "plink -pw genie156"
)
kongdd2 <- list(
  host = "kongdd-pc",
  user = "kongdd",
  rscript = "C:/Program Files/R/R-3.3.2/bin/Rscript.exe",
  rshcmd = "ssh"
)
kong <- list(
  host = "kong",
  user = "kongdd",
  rscript = "C:/Program Files/R/R-3.3.1/bin/Rscript.exe",
  rshcmd = "plink -pw genie156"
)

wang <- list(
  host = "wang-pc",
  user = "wang",
  port = 11044,
  rscript = "C:/Program Files/R/R-3.3.2/bin/Rscript.exe",
  rshcmd = "plink -pw 18314514006h"
)
# plink -pw 18314514006h wang@wang-pc taskkill /f /IM Rscript.exe
# plink -pw xiuyuan156 kongdd@zhang taskkill /f /IM Rscript.exe
# PsExec \\wang-pc -u wang -p 18314514006h ipconfig
# psservice \\wang-pc -u wang -p 18314514006h start sshd
# psservice \\kong -u kongdd -p genie156 start sshd
kizon <- list(
  host = "kizon",
  user = "kizon",
  rscript = "D:/program file/R/R-3.3.1/bin",
  rshcmd = "plink -pw 0312"
)
fkk <- list(
  host = "WIN-FANKKDESK",
  user = "tufcc",
  rscript = "C:/Program Files/R/R-3.3.0/bin/Rscript.exe",
  rshcmd = "plink -pw happy"
)
gu <- list(
  host = "GU",
  user = "xihui",
  rscript = "C:/Program Files/R/R-3.2.2/bin/Rscript.exe",
  rshcmd = "plink -pw zxq910221"
)
works <- list(zhang, wang, kong, kongdd, localhost) %>% set_names(sapply(., `[[`, "host"))
# cl <- makePSOCKcluster(names = rep(works, c(0, 0, 0, 1, 0)))
# cl <- makeCluster(6)
# cl <- makePSOCKcluster(names = works[2])

## test for kongdd-pc
cl <- makePSOCKcluster(names = rep(works, c(32, 0, 0, 0, 6)), outfile = "log.txt")
## test for wang-pc
# cl <- makePSOCKcluster(names = rep(works, c(0, 5, 0, 7, 0)))
# cl <- makePSOCKcluster(names = list(gu))