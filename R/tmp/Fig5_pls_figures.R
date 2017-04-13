library(stringr)

time_mode <- seq(ymd(sprintf("2010%02d01", beginMonth)), ymd("2011-12-31"), by = "month") 
time_fmt2 <- time_mode %>% format("%b"); time_fmt2[1:6] %<>% paste0("_0")
time_fmt <- month(time_mode)
## 抽取绘图需要的数据
#  foreach 大大简化了之前的代码
Q2 <- llply(plsr_fits, llply, `[[`, 'Q2')
Q2Info <- llply(Q2, ldply, function(x) summary(x[, "Q2"]))

coefs <- llply(plsr_fits, llply, `[[`, 'std.coefs')
VIP <- llply(plsr_fits, llply, function(fit) VIP_plsreg(fit) %>% {.[nrow(.), ]}, .progress = "text")

pls_list <- foreach(vipl = VIP, coefl = coefs, .final = .%>% set_names(clusters)) %:% 
  foreach(vip = vipl, coef = coefl, .final = .%>% set_names(methods)) %do% {
    names <- names(coef)
    vars <- substr(names, 1, 4) %>% factor(c("Tavg", "Tmax", "Tmin", "Prec", "Srad"))
    month <- str_extract(names, "(_|\\+)?\\d{1,4}") %>% {as.numeric(gsub("_", "-", .))}
    # month[month < 0] %<>% {abs(.) - 12}
    x <- data.frame(var = vars, xid = rep(1:(length(coef)/5), 5), 
                    month, coef, vip, row.names = NULL)
    x
  }
varnames <- c("var","xid", "month","coef","vip")
pls_df <- melt(pls_list, id.vars = varnames) %>%  set_names(c(varnames, "method", "cluster"))

clusterEN <- c("EBF", "Alpine Meadows", "Alpine Steppe",
               "Alpine Shrub Meadows", "Desert", "DBF",
               "Temperate Meadows", "Temperate Steppe", "TP")
# clusterEN <- c("Evergreen Broadleaf Forests", "Alpine Meadows", "Alpine Steppe", 
#                "Alpine Shrub Meadows", "Desert", "Deciduous Broadleaf Forest", 
#                "Temperate Meadows", "Temperate Steppe", "TP")
## 绘图主函数，先以分区为循环
set_fnames <- function(names){
  fun <- . %>% set_names(names)
  return(fun)
}

# for (i in seq_along(clusters)){
#   cluster <- clusters[i]
#   for (j in seq_along(methods)){
#     x <- pls_list[[i]][[j]]
#     pp <- vipCoef_plot(x, i, clusterEN[i], method)
#   }
# }
pp <- foreach(cluster = clusters, i = icount(), .final = set_fnames(clusters)) %:%
  foreach(method = methods, j = icount(), .final = set_fnames(methods)) %do% {
    x <- pls_list[[i]][[j]]
    vipCoef_plot(x, i, clusterEN[i], method)
  }
# p_need <- llply(p, `[[`, 1)

# library(gridExtra)
# # marrangeGrob(p, nrow = 3, ncol = 3)
# handle <- arrangeGrob(grobs = p, nrow = 2, ncol = 5)
# grid::grid.draw(handle)

# clusterEn <- c()
library(Cairo)

for (method in methods){
  cat(sprintf("[%s]", method), sep = "\n")
  p <- llply(pp, `[[`, method)
  handle <- gridExtra::arrangeGrob(grobs = p, nrow = 2, ncol = 5)

  CairoPDF(file = sprintf("Fig7_pls%s.pdf", method), width = 13.5, height = 8)
  grid::grid.draw(handle)
  dev.off()
}

# handle <- arrangeGrob(grobs = p, nrow = 2, ncol = 3);grid::grid.draw(handle)
