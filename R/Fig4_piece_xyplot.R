rm(list = ls())
source('R/mainfunc/mainfunctions.R', encoding = 'UTF-8')
source('R/mainfunc/Trends_func.R', encoding = 'UTF-8')
load("data/AVHRR_EOS_TP.rda")
# load("data/Fig1_agrdata2.rda")
load("data/00basemap.rda")
load("data/Fig1_infoValid.rda")

pheno <- phenoEOS[c(1, 3, 7, 9)]

# mean <- mean(pheno)
pheno_mean <- array(NA, dim=c(22204, 32, 4))
for (i in 1:4){
  pheno_mean[,, i] <- as.matrix(pheno[[i]])
}

pheno$Mean <- apply(pheno_mean, c(1,2), mean, na.rm = T) %>% as.data.frame()
pheno %<>% llply(function(x) set_colnames(as.matrix(x), NULL))
# x <- llply(pheno, apply, 1, mean, na.rm = T, .id = NULL, .progress = "text") %>% 
#   do.call(cbind.data.frame, .)

clusters <- c(levels(Id_cluster)[-7], "TP")
info <- llply(pheno, function(x) {
  df <- aggregate(x, list(cluster = Id_cluster), mean, na.rm = T)[-7, -1] %>% 
    {set_colnames(data.frame(year=1982:2013, t(.)), c("year", clusters[-9]))} %>% 
    set_rownames(NULL)
  df$TP <- apply(x, 2, mean, na.rm =T)
  return(df)#return
})

x <- melt(info, "year") %>% set_colnames(c("year", "cluster", "doys", "method"))
# x <- melt_list(info, "method")
ggplot(x, aes(year, doys, color = cluster)) + 
  geom_point() + geom_line() + 
  facet_grid(method~cluster, scales = "free_y") + 
  stat_smooth()

## 计算分区突变情况

info_doys <- lapply(info, function(x) x[, -1])
## 末期物候指标突变均不显著
psi <- llply(info_doys, function(x) ldply(x, piecewise, .id = NULL))

save(info_doys, file = "Fig4_infoDoys.rda")

library(foreach)
cluster_trends <- foreach(doys = info_doys, brks = psi) %:% 
  foreach(doy = doys, brk = brks$psi, .combine = rbind, .final =.%>%set_rownames(clusters)) %do% {
    MankKendall.Trend(doy, brk, trim = 2)
    # listk(doy, brk)
  } %>% set_names(names(info_doys))
#fix the problems of foreach package' dimnames, in future

x <- as.matrix(MeanDoys[, -(1:2)]); header <- MeanDoys[, 1:2]
## Piecewise changepoint detect information, pvalue < 0.1 and satisfy line regress x length
info.piece <- adply(x, 1, piecewise, .id = NULL)
psi.piece <- info.piece$psi; conf <- info.piece$conf
TrendInfo_piece <- MankKendall.Trend_apply(x, psi.piece, conf, 0.1)

## cpt.mean changepoint
info.cpt <- adply(x, 1, NDVI3g_cpt, .id = NULL)
psi.cpt <- info.cpt$mean; conf <- 1 - info.cpt$meanConf#small better
TrendInfo_cpt <- MankKendall.Trend_apply(x, psi.cpt, conf, 0.1)
## ------准备输入的绘图数据------------------------
## senslope in cluster regions
itemNames <- c("(a) 突变前趋势", "(b) 突变后趋势", 
               "(c) 突变指标全年趋势", "(d) 未突变指标全年趋势", 
               "(e) 全部指标全年趋势")
TrendInfo <- TrendInfo_piece
senslope <- TrendInfo[, c(5, 10, 16, 18, 15)] %>% cbind(header,.)
colnames(senslope) <- c("cluster", "doys", itemNames)
## 根据Z值判断显著性水平
breaks <- qnorm(c(0, 0.025,0.05, 0.5, 0.95, 0.975, 1))
sign <- TrendInfo[, c(4, 9, 20, 22, 14) -3] %>% apply(., 2, function(x) findInterval(x, breaks))
TrendSign <- cbind(header, sign) %>% set_colnames(c("cluster", "doys", itemNames))
dataPlot <- cbind(melt(senslope, id.vars= c("cluster", "doys"), value.name = "slp"), 
                  sign = melt(TrendSign, id.vars= c("cluster", "doys"))[, -(1:3)])
pch <- c(25, 6, NA, NA, 2, 24)
dataPlot$pch <- pch[dataPlot$sign]

at <- c(-Inf, -15, -5, -2)
cols_des <- RColorBrewer::brewer.pal(9, "YlOrRd")[c(3, 5, 6, 7)] #%>% c(., "red")
cols_inc <- RColorBrewer::brewer.pal(9, "YlGnBu")[c(3, 4, 6, 7)]
cols.sd <- c(rev(cols_inc), cols_des)
at <- c(at, 0, -rev(at))

ncol <- length(at) - 1
cols <- colorRampPalette(c("firebrick1" ,"orange3", "darkgoldenrod2", "grey90", 
                           brewer.pal(9, "YlGnBu")[c(4, 6, 7)]))(ncol)
cols <- cols[-5]; cols[ncol] <- "green4"

nat <- length(at)
dataPlot$slplevel <- cut(dataPlot$slp*10, at)
## piece表格输出
dlist <- split(dataPlot[, -3], dataPlot$variable)
writelist_ToXlsx(dlist, fname = "Fig5-1.xlsx")
metricsCH <- c("返青期", "成熟期", "衰老期", "休眠期", 
               "上升期", "稳定期", "下降期", "衰退期", 
               "DES生长开始时间", "DES生长季结束时间", "DES生长季长度", "峰值时间", 
               "TRS2生长开始时间", "TRS2生长季结束时间", "TRS2生长季长度",
               "TRS5生长开始时间", "TRS5生长季结束时间", "TRS5生长季长度")