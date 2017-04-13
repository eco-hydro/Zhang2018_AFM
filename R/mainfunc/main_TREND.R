library(segmented)
Rcpp::sourceCpp('src/mktrend.cpp')
# 分段回归PieceWise Regression ------------------------------------------------
piecewise <- function(y){#速度势必会很慢，同时并行运算数据在内存中传递同是浪费时间
  psi <- NA; conf <- NA
  if (length(unique(y)) > 3){#如果全部为0,则把该点数据剔除
    x <- seq_along(y)
    lm_fit <- lm(y ~x)
    #test beta2 significant
    dtest <- davies.test(lm_fit, ~x, k = 10)
    conf <- dtest$p.value#pvalue越小越好
    #如果dtest返回的有psi，则将此psi作为segmented的起始值
    psi_init <- tryCatch(floor(dtest$statistic[[1]]), error = function(e) 16)
    piece_fit <- segmented(lm_fit, seg.Z = ~x, psi = psi_init, 
                           control = seg.control(stop.if.error = FALSE, n.boot = 0, it.max = 100))
    if (("psi" %in% names(piece_fit)) && (length(piece_fit$psi) > 1)){#判断是否存在突变点
      psi <- floor(piece_fit$psi[2])
      #psi <- ifelse(is.na(psi), 16, psi)#断点不存在时处理
      ## 为保证分段回归数据长度：如果psi小于5，或者大于28则进行调整
    }else{
      psi <- tryCatch(floor(dtest$statistic[[1]]), error = function(e) NA)
    }
    # dati <- data.frame(x, y)
    # slope1 <- summary(lm(y ~x, subset(dati, x <= psi)))$coefficients[2, c(1, 4)]
    # slope2 <- summary(lm(y ~x, subset(dati, x > psi)))$coefficients[2, c(1, 4)]
  }
  ## 突变前后采用什么方法计算斜率另做打算
  data.frame(psi = psi, conf = conf)#, slope1, slope2, summary(lm_fit)$coefficients[2, c(1, 4)])#返回突变点，突变前后斜率，及对应的概率
}
## 转折点检验
# 为保证分段回归数据的长度，对突变点位置进行调整
cpt_adjust <- function(x) {
  if (x[1] < 5) {
    x[1] <- 5
  } else if (x[1] > 28) {
    x[1] <- 28
  }
  x  #quickly return
}

# fume package
# Z		The original (non corrected) Mann-Kendall test Z statistic
# p.value The original (non corrected) Mann-Kendall test p-value
# Zc		The new Z statistic after applying the correction
# Corrected p.value 	Corrected p-value after accounting for serial autocorrelation
# tau		Mann-Kendall's tau statistic
# N/n*s	Value of the correction factor, representing the quotient of the number of samples N divided by the effective sample size (n*s)
# Sen slope 			The slope of the (linear) trend according to Sen test

mkTrend.rcpp <- function(x, ci = 0.95) {
  x = x
  z0 = z = NULL
  pval0 = pval = NULL
  S = 0
  Tau = NULL
  essf = NULL
  ci = ci
  if (is.vector(x) == FALSE) stop("Input data must be a vector")
  if (any(is.infinite(x))) {
    x <- x[-which(is.infinite(x))]
    warning("The input vector contains non-finite numbers. An attempt was made to remove them")
  }
  n <- length(x)
  #20161211 modified, avoid x length less then 5, return rep(NA,5) c(z0, pval0, z, pval, slp)
  if (n < 5) return(rep(NA, 5))
  S = Sf(x)

  ro <- acf(rank(lm(x ~ I(1:n))$resid), lag.max = (n - 1), plot = FALSE)$acf[-1]
  sig <- qnorm((1 + ci)/2)/sqrt(n)
  rof <- ifelse(abs(ro) > sig, ro, 0)#modified by dongdong Kong, 2017-04-03
  
  cte <- 2/(n * (n - 1) * (n - 2))
  ess = 0
  for (i in 1:(n - 1)) {
    ess = ess + (n - i) * (n - i - 1) * (n - i - 2) * rof[i]
  }
  essf = 1 + ess * cte
  var.S = n * (n - 1) * (2 * n + 5)/18
  if (length(unique(x)) < n) {
    aux <- unique(x)
    for (i in 1:length(aux)) {
      tie <- length(which(x == aux[i]))
      if (tie > 1) {
        var.S = var.S - tie*(tie - 1)*(2*tie + 5)/18
      }
    }
  }
  VS = var.S * essf
  if (S == 0) {
    z = 0
    z0 = 0
  }
  if (S > 0) {
    z = (S - 1)/sqrt(VS)
    z0 = (S - 1)/sqrt(var.S)
  } else {
    z = (S + 1)/sqrt(VS)
    z0 = (S + 1)/sqrt(var.S)
  }
  pval = 2 * pnorm(-abs(z))
  pval0 = 2 * pnorm(-abs(z0))
  Tau = S/(0.5 * n * (n - 1))

  slp <- senslope(x)
  #return(list(Z = z0, p.value = pval0, Zc = z, `Corrected p.value` = pval, tau = Tau, `N/N*s` = essf, `Sen's Slope` = slp))
  #c(z = z0, p = pval0, zc = z, padj = pval, slp = slp)
  c(z0, pval0, z, pval, slp)#quickly return
} 

#20161211 modified, avoid x length less then 5
mkTrend <- function(x, ci = 0.95) {
    x = x
    z = NULL
    z0 = NULL
    pval = NULL
    pval0 = NULL
    S = 0
    Tau = NULL
    essf = NULL
    ci = ci
    if (is.vector(x) == FALSE) {
        stop("Input data must be a vector")
    }
    if (any(is.infinite(x))) {
        x <- x[-which(is.infinite(x))]
        warning("The input vector contains non-finite numbers. An attempt was made to remove them")
    }
    n <- length(x)
    #20161211 modified, avoid x length less then 5, return rep(NA,5) c(z0, pval0, z, pval, slp)
    if (n < 5) return(rep(NA, 5))
    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            S = S + sign(x[j] - x[i])
        }
    }
    ro <- acf(rank(lm(x ~ I(1:n))$resid), lag.max = (n - 1), plot = FALSE)$acf[-1]
    sig <- qnorm((1 + ci)/2)/sqrt(n)
    rof <- ifelse(abs(ro) > sig, ro, 0)#modified by dongdong Kong, 2017-04-03
    
    rof <- rep(NA, length(ro))
    for (i in 1:(length(ro))) {
        if (ro[i] > sig || ro[i] < -sig) {
            rof[i] <- ro[i]
        } else {
            rof[i] = 0
        }
    }
    cte <- 2/(n * (n - 1) * (n - 2))
    ess = 0
    for (i in 1:(n - 1)) {
        ess = ess + (n - i) * (n - i - 1) * (n - i - 2) * rof[i]
    }
    essf = 1 + ess * cte
    var.S = n * (n - 1) * (2 * n + 5) * (1/18)
    if (length(unique(x)) < n) {
        aux <- unique(x)
        for (i in 1:length(aux)) {
            tie <- length(which(x == aux[i]))
            if (tie > 1) {
                var.S = var.S - tie * (tie - 1) * (2 * tie + 5) * (1/18)
            }
        }
    }
    VS = var.S * essf
    if (S == 0) {
        z = 0
        z0 = 0
    }
    if (S > 0) {
        z = (S - 1)/sqrt(VS)
        z0 = (S - 1)/sqrt(var.S)
    } else {
        z = (S + 1)/sqrt(VS)
        z0 = (S + 1)/sqrt(var.S)
    }
    pval = 2 * pnorm(-abs(z))
    pval0 = 2 * pnorm(-abs(z0))
    Tau = S/(0.5 * n * (n - 1))
    V <- rep(NA, times = (n^2 - n)/2)
    k = 0
    for (i in 2:n) {
      for (j in 1:(i-1)){
        # for (j in 1:(n - 1)) {
            k = k + 1
            V[k] = (x[i] - x[j])/(i - j)
        }
    }
    slp <- median(na.omit(V))
    #return(list(Z = z0, p.value = pval0, Zc = z, `Corrected p.value` = pval, tau = Tau, `N/N*s` = essf, `Sen's Slope` = slp))
    #c(z = z0, p = pval0, zc = z, padj = pval, slp = slp)
    c(z0, pval0, z, pval, slp)#quickly return
} 

MankKendall.Trend <- function(x, brk, addName = T, trim = 1){
  if (trim == 1){
    brk[brk < 5 | brk > 27] <- NA
  }else if (trim == 2){
    brk[brk < 5] <- 5
    brk[brk > 27] <- 27
  }
    #cat(sprintf("Running %dth\n", i))
  info <- numeric(15)*NA
  if (!is.na(brk)){
      x_before <- x[1:brk]
      x_after <- x[(brk+1):length(x)]
      ## 还是需要对空值进行处理
      info <- c(mkTrend(x_before), mkTrend(x_after), mkTrend(x))
  }
  if (addName) 
      names(info) <- c("before.Z", "before.p", "before.Zc", "before.p.adj", "before.SenSlope", 
              "after.Z", "after.p", "after.Zc", "after.p.adj", "after.SenSlope", 
              "Z", "p", "Zc", "padj", "slp")
  info#quickly return
}
MankKendall.Trend_apply <- function(x, psi, conf, sign = NULL, trim = 1){
  n <- nrow(x)
  if (!is.null(sign)) psi[conf > sign] <- NA
  n <- length(psi)
  ## 做线性回归时序列长度保证至少为5
  
  TrendNames <- c("before.Z", "before.p", "before.Zc", "before.p.adj", "before.SenSlope", 
                  "after.Z", "after.p", "after.Zc", "after.p.adj", "after.SenSlope", 
                  "Z", "p", "Zc", "padj", "slp")
  TrendInfo <- data.frame(matrix(NA, nrow = length(psi), ncol = 15, dimnames = 
                                   list(NULL, TrendNames)))
  for (i in 1:n){
    if (mod(i, 100) == 0) 
      cat(sprintf("[%4d]th ------------\n", i))
    brk <- psi[i]
    if (!is.na(brk)){
      x_before <- x[i, 1:brk]
      x_after <- x[i, (brk+1):ncol(x)]
      ## 还是需要对空值进行处理
      TrendInfo[i, 1:10] <- c(mkTrend(x_before), mkTrend(x_after))
    }
    TrendInfo[i, 11:15] <- mkTrend(x[i, ])
  }
  TrendInfo#quickly return

  Id_brk <- !is.na(psi)
  nbrk_slp <- nbrk_z <- brk_slp <- brk_z <- numeric(length(Id_brk))*NA
  nbrk_slp[!Id_brk] <- TrendInfo$slp[!Id_brk]
  nbrk_z[!Id_brk] <- TrendInfo$Z[!Id_brk]
  brk_slp[Id_brk] <- TrendInfo$slp[Id_brk]
  brk_z[Id_brk] <- TrendInfo$Z[Id_brk]
  cbind(TrendInfo, brk_slp, brk_z, nbrk_slp, nbrk_z)#quickly return
}

## Fig5 Trend before and after TP: 
# Trend_spplot.grid(gridClip2, Trend.before,New_names, at.before)
Trend_spplot.grid <- function(gridClip, x, zcol, at, cols = NULL){
  if (!is.data.frame(x)) x <- as.data.frame(x)
  gridClip@data <- x
  nat <- length(at);
  if (is.null(cols)) cols <- col_fun(length(at) - 1)#col_fun:globalenv variable
  New_names <- zcol
  spplot(gridClip2, New_names, as.table = T, col.regions = cols, 
         sp.layout = sp.cluster, strip = F, 
         xlim = bbox(cluster)[1, ], ylim =  bbox(cluster)[2, ], 
         layout = c(4, 5), 
         colorkey = list(space = "bottom", height  = 0.8, title = "doys/decades", 
                         labels = list(labels = at[2:(nat - 1)], 
                                       at = 2:(nat - 1) - 0.5,
                                       cex = 1.1, fontfamily = "Times", fontface = 1),
                         axis.line = list(col = 'black')),
         #legend=list(list(fun=grid::textGrob("days\n(decades)", x=1.09))), 
         panel = function (x, y, z, subscripts, ...,  sp.layout) 
         {
           sppanel(list(sp.layout), panel.number(), first = TRUE)
           ## 转折不显著的设置成灰色
           Id_na <- is.na(z); zmark <- z; zmark[Id_na] <- 0; zmark[!Id_na] <- NA
           panel.levelplot.raster(x, y, zmark, subscripts, col.regions = "grey80", interpolate = T)
           panel.levelplot.raster(x, y, z, subscripts, ..., interpolate = T)

           i <- panel.number()
           panel.text(76.6, 39.3, paste0("(",letters[i], ") ", New_names[i]), fontfamily = "Times", cex = 1.2, adj = 0)
           sppanel(list(sp.layout), panel.number(), first = FALSE)
           panel.addbarchart(z, subscripts, cols, showCluster = F)
         }, par.settings = list(axis.line = list(col = "transparent")))
  ## add colorkey legend title
  # trellis.focus("legend", side="bottom", clipp.off=TRUE, highlight=FALSE)
  # grid::grid.text(expression(paste("(", day/decade, ")")), 0.907, 0.5, hjust=0, vjust=0, 
  #                 gp=gpar(fontfamily = "Times", cex = 1.2))
  # trellis.unfocus()
}