rm(list = ls())
source('R/mainfunc/mainfunctions.R', encoding = 'UTF-8')
source('R/mainfunc/Trends_func.R', encoding = 'UTF-8')
source('R/mainfunc/mainfunc_pls.R', encoding = 'UTF-8')
# load("data/AVHRR_EOS_TP.rda")
# load("data/Fig1_agrdata2.rda")
load("../data/00basemap.rda")
# load("data/Fig1_infoValid.rda")
load("../data/Fig4_infoDoys.rda")
if (!file.exists("../data/Fig5_plsMete.rda")){
  dataIn <- R.matlab::readMat("../data/plsMeteInputData.mat")
  Id_del <- c(16035, 18443, 20925, 21407, 21618, 21732, 21827)
  dataIn.trim <- lapply(dataIn, function(x) x[-Id_del, ])
  
  mete <- llply(dataIn.trim, function(x){
    result <- aggregate(x, list(cluster = Id_cluster), mean, na.rm = T)
    df <- data.frame(t(result[-7, -1])) %>% set_rownames(NULL) %>% set_colnames(levels(Id_cluster)[-7])
    df$TP <- apply(x, 2, mean, na.rm = T)
    return(df)
  }, .progress = "text")
  time_month <- seq(ymd("1981-01-01"), ymd("2013-12-01"), by = "month")
  save(mete, time_month, file = "../data/Fig5_plsMete.rda")
}else{
  load("../data/Fig5_plsMete.rda")
}

## global functions
beginMonth = 7
get_month <- function(doys) month(parse_date_time(paste0("2010", as.integer(doys)), "yj"))
## 准备气象数据, 上年7月至平均MOY
getPlsr_varId_monthly <- function(phen) lapply(1:32, function(i) beginMonth:(phen+12) + 12*(i - 1))

# main scripts ------------------------------------------------------------
doys <- info_doys[[1]]$TP
Month.mean <- ldply(info_doys, function(x) apply(x, 2, mean, na.rm=T) %>% 
                      sapply(get_month), .id = NULL) %>% set_rownames(names(info_doys))
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
## 同时需要把doy加进去

comps <- 10
plsr_fits <- foreach(doy_cluster=doys, metes=dataIn.pls, i = icount, .final = .%>% set_names(clusters)) %:%
  foreach(Y=doy_cluster, X=metes, .final = .%>% set_names(methods)) %do%{#
    # Y
    # print(str(X))
    plsr_fit = plsreg1(X, scale(Y), comps=comps, crosval=T)
    plsr_fit
  }
Q2 <- llply(plsr_fits, llply, `[[`, 'Q2')
info <- llply(Q2, ldply, function(x) summary(x[, "Q2"]))
# 
# plsr_fit = plsreg1(X, scale(Y), comps=comps, crosval=T)
# 
# library(ropls)
# opls(X, Y, plotL = FALSE)


# plsrdata <- data.frame(Y = Y)
# plsrdata$X <- Xnew
# comps <- 10
# plsr_fit = plsreg1(X, Y, comps=comps, crosval=T)
# coef <- plsr_fit$std.coefs
# vip <- VIP_plsreg(plsr_fit); vip <- vip[nrow(vip), ]
# 
# if (crossVal){
#   Q2 <- plsr_fit$Q2
#   list(coef = coef, vip = vip, Q2 = Q2)
# }else {
#   data.frame(coef = coef, vip = vip, x = seq_along(vip), 
#              xid = rep(1:(length(vip)/7), 7), 
#              mete = rep(MeteNames[-4], rep(length(vip)/7, 7)))#quickly return
# }
# 
# ## 测试plsr与plsreg1计算结果是否一致
# data <- data.frame(doy); data$x <- as.matrix(X)
# data$x <- scale(as.matrix(X))
# pls_fit0 <- plsr(doy~x, 10, data = data)
