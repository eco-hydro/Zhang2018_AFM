rm(list = ls())
source('R/mainfunc/main.R', encoding = 'UTF-8')
source('R/mainfunc/main_PLS.R', encoding = 'UTF-8')

load("data/phenology_list20170105.rda")
load("../TP_phenoAgr/agrmet25_TP.rda")
load("../data/00basemap.rda")
load("../data/Fig5_plsMete.rda")

## global functions
get_month <- function(doys) month(parse_date_time(paste0("2010", as.integer(doys)), "yj"))
## 准备气象数据, 上年7月至平均MOY
getPlsr_varId_monthly <- function(phen) lapply(1:32, function(i) beginMonth:(phen+12) + 12*(i - 1))

## ------------------ PLS模型物候数据准备 ---------------------
SOS_sel <- SOS_list[c(1, 11, 13, 15)]
EOS_sel <- EOS_list[c(3, 5, 11, 13)]

## only mean phenology data remained
getMean_pheno <- function(pheno){
  # mean <- mean(pheno)
  pheno_mean <- array(NA, dim=c(22204, 32, 4))
  for (i in 1:4){
    pheno_mean[,, i] <- as.matrix(pheno[[i]])
  }
  apply(pheno_mean, c(1,2), mean, na.rm = T) %>% round()
  # pheno$Mean <- apply(pheno_mean, c(1,2), mean, na.rm = T) %>% as.data.frame()
  # pheno %<>% llply(function(x) set_colnames(as.matrix(x), NULL))
  # return(pheno)
}
pheno <- list(SOS = getMean_pheno(SOS_sel), EOS = getMean_pheno(EOS_sel))

# pheno$EOS <- as.matrix(EOS_sel$Gu.beck.DD)#test previous error

clusters <- c(levels(Id_cluster)[-7], "TP")
info <- llply(pheno, function(x) {
  df <- aggregate(x, list(cluster = Id_cluster), mean, na.rm = T)[-7, -1] %>% 
  {set_colnames(data.frame(year=1982:2013, t(.)), c("year", clusters[-9]))} %>% 
    set_rownames(NULL)
  df$TP <- apply(x, 2, mean, na.rm =T)
  return(df)#return
}, .progress = "text")
#EOS_1 present previous annual EOS phenology metric
info$EOS_1 <- cbind(year = 1982:2013, rbind(NA, info$EOS[-1, -1]))

x <- melt(info, "year") %>% set_colnames(c("year", "cluster", "doys", "method"))
## plot to check
ggplot(subset(x, method == "EOS"), aes(year, doys, color = cluster)) + 
  geom_point() + geom_line() + 
  facet_wrap( ~ cluster) + 
  stat_smooth(method ="loess")

info_doys <- lapply(info, function(x) x[, -1])#cluster mean doys, used to construct pls model

## ---------------- PLS模型数据预处理 -------------------------
## 根据EOS挑选pls变量个数
Month.mean <- ldply(info_doys, function(x) apply(x, 2, mean, na.rm=T) %>% sapply(get_month), .id = NULL) %>% 
  set_rownames(names(info_doys))
methods <- rownames(Month.mean)
clusters <- colnames(Month.mean)
# [1] "常绿阔叶林"     "高寒草甸"       "高寒草原"       "高寒灌木、草甸" "荒漠"           "落叶阔叶林"    
# [7] "温带草原"       "温性草原"       "TP"   

## 首先转化mete格式
mete.list <- lapply(clusters, function(cluster) llply(mete, `[[`, cluster) %>% do.call(cbind.data.frame,.)) %>% set_names(clusters)

dataIn <- foreach(x = mete.list, doys_Mean = Month.mean) %do%{
  # doys_Mean
  Id <- llply(doys_Mean, getPlsr_varId_monthly) %>% set_names(methods)
  llply(Id, function(Idi) lapply(Idi, function(i)x[i,]) %>% set_names(1:32))
} %>% set_names(clusters)

dataIn.pls <- llply(dataIn, function(xlist) {
  llply(xlist, function(x) {
    nvar <- nrow(x[[1]])
    vars <- paste0(rep(colnames(x[[1]]), rep(nvar, 5)), c(paste0("_", 7:12), 1:(nvar - 6)))
    xx <- ldply(x, function(x) as.numeric(as.matrix(x)), .id = NULL) %>% set_rownames(1982:2013) %>% set_colnames(vars)
    xx
  })
})

doys <- llply(clusters, function(cluster) {
  x <- ldply(info_doys, `[[`, cluster, .id=NULL)
  as.data.frame(t(x)) %>% set_colnames(methods) %>% set_rownames(1982:2013)
}) %>% set_names(clusters)

## -----------------------------------------------------------------
## 采用新方法进行改进
source('R/mainfunc/main_PLS.R', encoding = 'UTF-8')

## debug version
# for (i in seq_along(clusters)){
#   X = cbind(doy_cluster[, c("SOS", "EOS_1")],  metes$EOS)
#   Y = doy_cluster$EOS %>% matrix(ncol = 1)
#   print(i)
#   Id_nona <- which(complete.cases(X) & complete.cases(Y))
#   suppressMessages(valid_pls(X[Id_nona, ], Y[Id_nona], minVar = 10, method = "VIP", plsModel = F))
# }

## minval = 10
# with sos
v <- list()
# no sos
v[[1]] <- PLS_main.SOS(doys, dataIn.pls, Idetrend = TRUE, hSOS = FALSE)
v[[2]] <- PLS_main.SOS(doys, dataIn.pls, Idetrend = FALSE, hSOS = FALSE)

names(v) <- c("noDetrend", "Detrend")
validInfo <- llply(v, `[[`, "validInfo")
coef <- llply(v, `[[`, "coef")
coef_out <- llply(coef, function(x) ldply(x, function(x) x[1, ]))

writelist_ToXlsx(validInfo, "PLS_validInfo.xlsx")
writelist_ToXlsx(coef_out, "SOS_Info.xlsx")



str(v, 1)

lapply(v[[1]]$coef, `[`, i = 1, j = ) %>% melt_list("cluster")
lapply(v[[2]]$coef, `[`, i = 1, j = ) %>% melt_list("cluster")

xlist <- EVAP$x
Id <- sample(seq_along(xlist), 10)
x <- xlist[Id]


## traditional PLS method test
hSOS = TRUE
Idetrend = TRUE
methods <- "Mean"
# Idetrend <- TRUE#FALSE
# hSOS <- TRUE
filename <-  c("Fig7_plsmean", if (hSOS) "SOS" else NULL, 
               if (Idetrend) "detrend" else NULL) %>% 
  paste(collapse = "_", sep = "") %>% paste0(., ".pdf")

progress <- create_progress_bar("text")
progress$init(length(clusters) * length(methods));

PLS_list <- foreach(doy_cluster=doys, metes=dataIn.pls, .final = .%>% set_names(clusterEN), i = icount()) %do% {
  if (hSOS){
    X <- cbind(doy_cluster[, c("SOS"), drop = F],  metes$EOS)
    # X = cbind(doy_cluster[, c("SOS", "EOS_1")],  metes$EOS)
  }else{
    X <- metes$EOS
  }
  Y = doy_cluster$EOS %>% matrix(ncol = 1)
  
  if (Idetrend){
    X %<>% {pracma::detrend(as.matrix(.))}
    Y %>% pracma::detrend()
  }
  Id_nona <- which(complete.cases(X) & complete.cases(Y))
  
  progress$step()
  suppressMessages(valid_pls_tradition(X[Id_nona, ], Y[Id_nona], comps = 2))
}

validInfo <- lapply(PLS_list, `[[`, "validInfo") %>%
  {melt_list(., "cluster")[, c(9, 1:5, 7:8)]}

openxlsx::write.xlsx(subset(validInfo, predI == 2), 
                     "validInfo_tradition.xlsx")



