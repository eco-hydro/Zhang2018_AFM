tmp <-clusterEvalQ(cl, {
  library(magrittr)
  library(zoo)
  library(plyr)
  library(optimx)
  library(zoo)
  #close all warnings, due to optimx check initial start point, can ignore it
  options(warn =-1)
})

Id_host <- selectHost(cl)
hosts <- names(Id_host)
indir <- "E:/Phenology/"

dir.status <- clusterSCall(cl, check_dir, indir) %>% print()
# works.sysinfo <- clusterLCall(cl, sysinfo, localhost)
## find which Id grid have been finished, and then skip it in parallel
Id <- read.table("Id_cal.txt", header = T)$x
files_finish <- clusterSCall(cl, get_files, indir=indir, pattern="*.txt", name=F)
if (length(files_finish) > 0){
  Id_finish <- stringr::str_extract(files_finish, "\\d{1,}") %>% as.numeric()
  Id_cal <- setdiff(Id, Id_finish)
}else{
  Id_cal <- Id
}
write.table(Id, file = "Id_cal.txt", sep = "\n", row.names = F)
# Id_cal <- read.table("Id_cal.txt", header = T)$x
# write.table(Id_cal[])
## save phenological data into txt files in case of script crash and lost everything
clusterExport(cl, c("dates_avhrr.doy", "fprintf", funcs), envir = globalenv())
RESULT <- parLapplyLB(cl, NDVI_list[Id_cal], getYears_phenology2, dates_avhrr.doy, indir)
# invoid to use plink and stopcluster frequently, because waste a lot of time
# stopCluster(cl)

# shell("plink -pw genie156 kongdd@pc taskkill /f /IM Rscript.exe", intern = T)
# shell("plink -pw genie156 kongdd@pc taskkill /f /IM Rscript.exe", intern = T)
