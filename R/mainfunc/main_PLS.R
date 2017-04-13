library(pls)
# library(plsdepot)#already rewrite its PLSREG1 function
library(ropls)
source("E:/GitHub/statistics_Kong.R", encoding = "utf-8")

## ------------------------ GLOBAL VARIABLES -------------------------------
clusterEN <- c("EBF", "Alpine Meadow", "Alpine Steppe",
                 "Alpine Shrub Meadow", "Desert", "DBF",
                 "Temperate Meadow", "Temperate Steppe", "TP")

beginMonth = 7
time_mode <- seq(ymd(sprintf("2010%02d01", beginMonth)), ymd("2011-12-31"), by = "month") 
  time_fmt2 <- time_mode %>% format("%b"); time_fmt2[1:6] %<>% paste0("_0")
  time_fmt <- month(time_mode)

# vars_order <- vars_new[c(4:15, 1:3)]
# vars <- vars_new

## ------------------------------- GLOBAL FUNCTIONS --------------------------
set_fnames <- function(names) . %>% set_names(names)

## ------------------------ Partial least regression -------------------------
VIP_pls <- function(object) {
  # if (object$method != "oscorespls")
  #     stop("Only implemented for orthogonal scores algorithm.  Refit with 'method = \"oscorespls\"'")
  # if (nrow(object$Yloadings) > 1)
  #     stop("Only implemented for single-response models")
  SS <- c(object$Yloadings)^2 * colSums(object$scores^2)
  Wnorm2 <- colSums(object$loading.weights^2)
  SSW <- sweep(object$loading.weights^2, 2, SS / Wnorm2, "*")
  sqrt(nrow(SSW) * apply(SSW, 1, cumsum) / cumsum(SS))
}
## 对于plsdepot package, 返回的plsr model计算VIP
VIP_plsreg <- function(object){
  # Writed By Dongdong kong, 2016-05-30
  SS <- c(object$y.loads)^2 * colSums(object$x.scores^2)
  Wnorm2 <- colSums(object$raw.wgs^2)
  SSW <- sweep(object$raw.wgs^2, 2, SS / Wnorm2, "*")
  sqrt(nrow(SSW) * apply(SSW, 1, cumsum) / cumsum(SS))
}

## calculate each variables SSE, 2017-01-14
SSX_plsreg1 <- function(fit, X){
  xMN <- scale(X)
  ssxTotN <- sum(xMN^2, na.rm = TRUE)
  xpre <- with(fit, x.scores %*% t(x.loads))
  # R2X <- sum(xpre^2)/ssxTotN
  apply(xpre, 2, function(x) sum(x^2))
}

SSX_ropls <- function(fit, X){
  xMN <- scale(X)
  ssxTotN <- sum(xMN^2, na.rm = TRUE)
  xpre <- fit@scoreMN %*% t(fit@loadingMN)
  # R2X <- sum(xpre^2)/ssxTotN
  #contribution of each variables to SSX
  apply(xpre, 2, function(x) sum(x^2))
}

## plsr cross validation and fit
plsr_crossVal <- function(dl, crossVal = T){
  comps <- ifelse(crossVal, 10, dl$comps)
  #if ('comps' %in% names(dl)) comps = dl$comps else comps = 10

  Xnew <- dl$Xnew; Y <- dl$Y
  # plsrdata <- data.frame(Y = Y)
  # plsrdata$X <- Xnew
  plsr_fit = plsreg1(Xnew, Y, comps=comps, crosval=T)
  coef <- plsr_fit$std.coefs
  vip <- VIP_plsreg(plsr_fit); vip <- vip[nrow(vip), ]

  if (crossVal){
    Q2 <- plsr_fit$Q2
    list(coef = coef, vip = vip, Q2 = Q2)
  }else {
    data.frame(coef = coef, vip = vip, x = seq_along(vip), 
      xid = rep(1:(length(vip)/7), 7), 
      mete = rep(MeteNames[-4], rep(length(vip)/7, 7)))#quickly return
  }
}

## ----------------------- PLS select variables -------------------------
#recalculate RMSEE to promote the precision
get_RMSEE <- function(x){
  # .errorF <- function(x, y) sqrt(mean(drop((x - y)^2), na.rm = TRUE))
  # RMSEE <- sqrt(.errorF(x@suppLs$yMCN, x@suppLs$yPreMN)^2 * nrow(X)/(nrow(X) - x@summaryDF$pre - 1))
  RMSEE <- sum((x@suppLs$yMCN - x@suppLs$yPreMN)^2, na.rm = TRUE)/(nrow(X) - x@summaryDF$pre - 1)
  return(sqrt(RMSEE))
  
  # PRESS <- sum((x@suppLs$yMCN - x@suppLs$yPreMN)^2, na.rm = TRUE)
  # Yorigin <- x@suppLs$yMCN
  # SS <- sum((Yorigin - mean(Yorigin))^2)
  # 1 - PRESS/SS
  
  # #this two result is same
  # Ypre <- x@suppLs$xModelMN %*% tcrossprod(x@weightStarMN, x@cMN)
  # Ypre <- x@scoreMN %*% t(x@cMN)
  # #adjust mean and std
  # YpreR <- scale(scale(Ypre, FALSE, 1/x@ySdVn), -x@yMeanVn, FALSE)
}

valid_opls <- function(X, Y){
  ## 选择一种方法去剔除冗余变量
  errorId <- list()
  orthoI <- 0
  printL <- plotL <- FALSE
  
  nvar <- rev(2:ncol(X))
  fits <- foreach(i = nvar, k = icount()) %do%{
    #acording to PRESS and Q2 select comps
    opls_fit <- tryCatch(opls(X, Y, orthoI = orthoI, plotL = plotL, printL = printL), 
                         error = function(e) {
                           message(sprintf("[%d] %s", i, e))
                           errorId[[k]] <<- k
                           tryCatch(opls(X, Y, predI = 2, plotL = plotL, printL = printL),
                                    error = function(e) message(sprintf("[%d] %s", i, e)))
                         })
    
    SSXi <- SSX_ropls(opls_fit, X)
    Id_del <- which.min(SSXi)
    X <- X[, -Id_del, drop = F]
    opls_fit
  }
  fits %<>% set_names(nvar)
  Id_del <- unlist(errorId)
  if (length(Id_del)) fits <- fits[-Id_del]
  
  Q2 <- ldply(fits, function(x) x@summaryDF, .id = "nvar")
  Q2$RMSEE <- laply(fits, get_RMSEE)
  
  return(Q2)
} 

valid_pls_tradition <- function(X, Y, comps = 2){
	# method <- toupper(method[1])
	plsr_fit <- plsreg1(X, Y, comps = comps)
	Info <- data.frame(nvar = ncol(X), predI = 1:2, plsr_fit$Q2) %>% 
			mutate(PRESS_adj = PRESS/(nrow(X) - predI - 1))
	vip <- VIP_plsreg(plsr_fit) %>% {.[nrow(.), ]}
	coef <- plsr_fit$std.coefs

	coef = data.frame(vip, coef)#RETURN
	return(list(coef = coef, validInfo = Info))
}

## 继续强化，查看opls是否表现出众
#  检查valid_pls每次计算的结果为什么会产生波动
valid_pls <- function(X, Y, minVar = 10, method = c("SSX", "VIP"), plsModel = FALSE){
	method <- toupper(method[1])

	nvar <- rev(minVar:ncol(X)); vars <- colnames(X)
	errorId <- list()

	fits <- foreach(i = nvar, k = icount()) %do%{
	#acording to PRESS and Q2 select comps
	plsr_fit <- tryCatch(
	  plsreg1(X, Y, comps = NULL),
	  error = function(e) {
	    errorId[[k]] <<- k
	    # message(sprintf("[%d]error: ", i))
	    message(sprintf("[%d] %s", i, e))
	    plsreg1(X, Y, comps = 2)
	  }
	)
	#delete unimportant variables based on "SSX" or "pls"
	if (method == "SSX"){
	  SSXi <- SSX_plsreg1(plsr_fit, X)
	  Id_del <- which.min(SSXi)
	}else if(method == "VIP"){
	  Id_del <- VIP_plsreg(plsr_fit) %>% {which.min(.[nrow(.), ])}
	}else{
	  stop("method should be 'SSX' or 'VIP'!\n")
	}

	X <- X[, -Id_del, drop = F]
	plsr_fit
	}
	fits %<>% set_names(nvar)
	Id_del <- unlist(errorId)
	# if (length(Id_del)) fits <- fits[-Id_del]

	Q2 <- llply(fits, `[[`, "Q2") %>% ldply(function(x) cbind.data.frame(predI = 1:nrow(x), x), .id = "nvar")
	if (nrow(Q2) > 0) Q2$nvar %<>% {as.numeric(as.character(.))}

	#according to Q2 select significant variables
	get_sign <- function(x){
	  limq2 <- 0.05
	  Id <- which(x$Q2 > limq2)
	  x[Id, ]
	}
	info <- ddply(Q2, .(nvar), get_sign)
	# if (nrow(info) > 0) info %<>% mutate(PRESS_adj = PRESS/(nrow(X) - predI - 1))
	#delete nvars which not all components Q2 is significant
	if (nrow(info) > 0){
		info %<>% mutate(PRESS_adj = PRESS/(nrow(X) - predI - 1))
		Id_nvar <- split(info$predI, info$nvar) %>% {which(laply(., function(x) 1 %in% x))}
		if (length(Id_nvar) == 0) {
		  warning("all fits the first components is not reliable!")
		  Id_nvar <- 1:nrow(info)
		}
		info <- subset(info, nvar %in% unique(info$nvar)[Id_nvar]) %>% set_rownames(NULL)
		# info <- info[Id_var, ]
	}
	# if (plsModel) {
    Id.best <- which(nvar == info$nvar[which.min(info$PRESS_adj)])
    plsr_fit <- fits[[Id.best]]
    vip <- VIP_plsreg(plsr_fit) %>% {.[nrow(.), ]}
    coef <- plsr_fit$std.coefs
    COEF <- VIP <- {numeric(length(vars))*NA} %>% set_names(vars)
    Id.var <- match(names(vip), vars)
    VIP[Id.var] <- vip
    COEF[Id.var] <- coef
    coef = data.frame(vip = VIP, coef = COEF)#RETURN
	# }else{
	# info#RETURN
	# }
	return(list(coef = coef, validInfo = info))
}

VIP_tidy <- function(x, hSOS = FALSE){
  if (hSOS) x <- x[-1, ]
  names <- rownames(x)
  vars <- substr(names, 1, 4) %>% factor(c("Tavg", "Tmax", "Tmin", "Prec", "Srad"))
  month <- stringr::str_extract(names, "(_|\\+)?\\d{1,4}") %>% {as.numeric(gsub("_", "-", .))}
  # month[month < 0] %<>% {abs(.) - 12}
  data.frame(var = vars, xid = rep(1:(nrow(x)/5), 5), 
                  month, x, row.names = NULL)#return
}

## get significant PLS variables id group, namely vip > 1
get_GroupId <- function(x, critical = 1){
  # sign <- (vip >= 1) %>% as.numeric
  vip <- x$vip; vip[is.na(vip)] <- 0
  runs <- rle(vip >= critical)
  runs$values[runs$values] <- cumsum(runs$values[runs$values])#continue serise tag
  group <- inverse.rle(runs)
  group[group == 0] <- NA

  # Id <- which(vip >= 1); group <- rleid(vip >= 1)
  # runs$lengths[runs$values] < 5
  # group[vip < 1] <- NA 
  # group#quickly return
  ## return continue serise begin position and end position
  tmp <- aggregate(seq_along(group), list(group = group), function(x) c(xmin = min(x) - 0.5, xmax = max(x) + 0.5))
  if (nrow(tmp) > 0){
    data.frame(xid = 1:nrow(tmp), group =tmp$group, tmp$x)
  }else{
    data.frame(xid = 1, group = 1, xmin=1, xmax = 1)[-1, ]
  }
}

vipCoef_plot <- function(x, i=1, cluster = "", method="", critical = 1){
  # PLS variables significant rectangle data
  sign <- ddply(x, .(var), get_GroupId, critical)
  yrange <- range(x$coef, na.rm = T)*100
  ylims <- c(floor(yrange[1]), ceiling(yrange[2]))/100
  sign %<>% data.frame(ymax = ylims[2], ymin = ylims[1])
  sign$var %<>% factor(c("Tavg", "Tmax", "Tmin", "Prec", "Srad"))
  
  nmonth <- nrow(x)/5
  breaks <- seq(1, nmonth, 2)
  breaks_labels <- time_fmt[seq(1, nmonth, 2)]
  delta <- if (length(breaks) == 10) c(-0.1, 0.1) else c(0, 0)
  
  	## adjust colors
	cols <- matrix(c(244, 117, 85, 
	          164, 174, 39, 
	          0, 168, 108, 
	          0, 157, 227, 
	          157, 64, 131), byrow = T, ncol = 3)/255 
	cols <- rgb(cols[, 1], cols[, 2], cols[, 3])

  p <- ggplot(x) + 
    geom_bar(aes(xid, coef, fill= var), stat = "identity", position = "identity", width = 0.85) + 
    # geom_area(aes(xid, coef, fill = mete)) +
    #, scales = "free_y", modified 20170216
    facet_grid(var~.) + 
    geom_hline(yintercept = 0, color = "grey10") +
    geom_vline(xintercept = 6.5, color = "red", lty = 2) + 
    geom_rect(data = sign, aes(xmin = xmin, xmax=xmax, ymin = ymin, ymax= ymax), 
              fill = "lightblue",color="grey40", alpha=0.2, size = 0.4) + 
    scale_x_continuous(breaks = breaks, labels = breaks_labels) +
    scale_y_continuous(breaks = c(ylims[1], ylims[2]), labels = c(ylims[1], ylims[2])) + 
    scale_fill_manual(values = cols) + #20170316
    coord_cartesian(xlim = c(1, nmonth) + delta) +
    guides(fill = F) + 
    ggtitle(sprintf('(%s) %s', letters[i], cluster)) + 
    theme(strip.text = element_blank(), 
          axis.title = element_blank(), 
          plot.title = element_text(family = "Times", hjust = -0.2, size = 14),
          axis.text.x = element_text(angle = 0), 
          # plot.title = element_text(family  = "ST", size = 14, hjust = 0), 
          plot.margin = unit(c(1,1,1,1)*0.1, "cm"))
    # guides(fill = guide_legend(nrow = 1)) + ylab('Model coefficient') + 
    # theme(strip.text = element_blank()) + 
    
    # ggtitle(sprintf('(%s) %s-%s', letters[i], cluster, method))
  return(p)
}

PLS_main <- function(doys, dataIn.pls, Idetrend = TRUE, hSOS = TRUE){
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
    suppressMessages(valid_pls(X[Id_nona, ], Y[Id_nona], minVar = 10, method = "VIP", plsModel = F))
  }

  validInfo <- lapply(PLS_list, `[[`, "validInfo")
  plsCoef <- lapply(PLS_list, `[[`, "coef")
  phenoCoef <- lapply(plsCoef, `[`, i = 1:2, j=) %>% melt_list("cluster") %>% cbind(var = c("SOS", "EOS_1"), .)

  validInfo.opt <- ldply(validInfo, function(x) subset(x, nvar == x$nvar[which.min(x$PRESS_adj)]), .id = "cluster") %T>% print
  plsCoef %>% lapply(na.omit)

  plsCoef.tidy <- llply(plsCoef, VIP_tidy, hSOS)

  ## 同时把上年EOS和当年SOS放进自变量，查看植被自身物候的变化or气候变化对其的影响哪个较大
  suppressWarnings({
      pp <- foreach(cluster = clusters, vip = plsCoef.tidy, i = icount()) %do% 
      vipCoef_plot(vip, i, clusterEN[i], critical = 0.8)

      pp <- pp[-9]
      # CairoPDF(file = filename, width = 13.5, height = 8)
      CairoPDF(file = filename, width = 11, height = 7)
      handle <- gridExtra::arrangeGrob(grobs = pp, nrow = 2)
      grid::grid.draw(handle)
      dev.off()
    })
  return(list(coef = plsCoef, validInfo = validInfo.opt))
}

PLS_main.SOS <- function(doys, dataIn.pls, Idetrend = TRUE, hSOS = TRUE){
  methods <- "Mean"
  # Idetrend <- TRUE#FALSE
  # hSOS <- TRUE
  filename <-  c("Fig7_plsmean", if (hSOS) "SOS" else NULL, 
                 if (Idetrend) "detrend" else NULL) %>% 
    paste(collapse = "_", sep = "") %>% paste0(., ".pdf")

  progress <- create_progress_bar("text")
  progress$init(length(clusters) * length(methods));

  for (i in seq_along(doys)){
    doy_cluster = doys[[i]]
    metes = dataIn.pls[[i]]
    
    X <- metes$SOS
    Y = doy_cluster$SOS %>% matrix(ncol = 1)
    
    if (Idetrend){
      X %<>% {pracma::detrend(as.matrix(.))}
      Y %>% pracma::detrend()
    }
    Id_nona <- which(complete.cases(X) & complete.cases(Y))
    
    # progress$step()
    print(i)
    tmp <- suppressMessages(valid_pls(X[Id_nona, ], Y[Id_nona], minVar = 10, method = "VIP", plsModel = F))
  }
  PLS_list <- foreach(doy_cluster=doys, metes=dataIn.pls, .final = .%>% set_names(clusterEN), i = icount()) %do% {
    X <- metes$SOS
    Y = doy_cluster$SOS %>% matrix(ncol = 1)

    if (Idetrend){
      X %<>% {pracma::detrend(as.matrix(.))}
      Y %>% pracma::detrend()
    }
    Id_nona <- which(complete.cases(X) & complete.cases(Y))
    
    progress$step()
    print(i)
    suppressMessages(valid_pls(X[Id_nona, ], Y[Id_nona], minVar = 10, method = "VIP", plsModel = F))
  }

  validInfo <- lapply(PLS_list, `[[`, "validInfo")
  plsCoef <- lapply(PLS_list, `[[`, "coef")
  phenoCoef <- lapply(plsCoef, `[`, i = 1:2, j=) %>% melt_list("cluster") %>% cbind(var = c("SOS", "EOS_1"), .)

  validInfo.opt <- ldply(validInfo, function(x) subset(x, nvar == x$nvar[which.min(x$PRESS_adj)]), .id = "cluster") %T>% print
  plsCoef %>% lapply(na.omit)

  plsCoef.tidy <- llply(plsCoef, VIP_tidy, hSOS)

  ## 同时把上年EOS和当年SOS放进自变量，查看植被自身物候的变化or气候变化对其的影响哪个较大
  suppressWarnings({
      pp <- foreach(cluster = clusters, vip = plsCoef.tidy, i = icount()) %do% 
      vipCoef_plot(vip, i, clusterEN[i], critical = 0.8)

      pp <- pp[-9]
      # CairoPDF(file = filename, width = 13.5, height = 8)
      CairoPDF(file = filename, width = 11, height = 7)
      handle <- gridExtra::arrangeGrob(grobs = pp, nrow = 2)
      grid::grid.draw(handle)
      dev.off()
    })
  return(list(coef = plsCoef, validInfo = validInfo.opt))
}

## 对该程序进行改进，优化挑选变量的方法
#  Writed By Dongdong Kong, 2017-01-07
#  ------------------------------------
#  origin limq2 = 0.0975, 2017-01-15
plsreg1 <- function (predictors, response, comps = 2, crosval = TRUE, limq2 = 0.05) 
{
    X = as.matrix(predictors)
    n = nrow(X)
    p = ncol(X)
    if (p < 2) stop("predictors must contain more than one column")
    if (is.null(colnames(X))) colnames(X) = paste0("X", 1:p)
    if (is.null(rownames(X))) rownames(X) = 1:n
    Y = as.matrix(response)
    if (ncol(Y) != 1) stop("\nresponse must be a single variable")
    if (any(is.na(response))) stop("\nresponse must not contain missing values")
    if (nrow(X) != nrow(Y)) stop("\npredictors and response have different number of rows")
    if (is.null(colnames(Y))) colnames(Y) = "Y"
    if (is.null(rownames(Y))) rownames(Y) = 1:n
    na.miss <- ifelse(any(is.na(X)), TRUE, FALSE)

    if (!is.null(comps)) {
        nc = comps
        #nc >= 1, modified 2017-01-15
        if (mode(nc) != "numeric" || length(nc) != 1 || nc <= 1 || (nc%%1) != 0 || nc > min(n, p)) nc = min(n, p)
        if (nc == n) nc = n - 1
    } else {
        if (na.miss) {
            crosval = FALSE
            nc = 2
        } else {
            crosval <- ifelse(n>=10, TRUE, FALSE)
            nc = min(n, p)
        }
    }
    if (!is.logical(crosval)) crosval = FALSE
    Xx = scale(X)
    Yy = scale(Y)
    ## 进行计算时已经进行了标准化处理
    X.old = Xx
    Y.old = Yy
    Hot <- Th <- Uh <- matrix(NA, n, nc)
    Ph <- Wh <- matrix(NA, p, nc)
    ch = rep(NA, nc)
    hlim = rep(NA, nc)
    if (crosval) {
        RSS = c(n - 1, rep(NA, nc))
        PRESS = rep(NA, nc)
        Q2 = rep(NA, nc)
        sets_size = c(rep(n%/%10, 9), n - 9 * (n%/%10))
        #Q2, PRESS et.al, changes' reason
        obs = 1:n
        # obs = sample(1:n, size = n)#2017-01-15
        segments = vector("list", length = 10)
        ini = cumsum(sets_size) - sets_size + 1
        fin = cumsum(sets_size)
        for (k in 1:10) segments[[k]] = obs[ini[k]:fin[k]]
    }
    w.old = rep(1, p)
    t.new = rep(1, n)
    p.new = rep(NA, p)
    h = 1
    repeat {
        if (na.miss) {
            for (j in 1:p) {
                i.exist = which(complete.cases(X[, j]))
                w.old[j] = sum(X.old[i.exist, j] * Y.old[i.exist])
            }
            w.new = w.old/sqrt(sum(w.old^2))
            for (i in 1:n) {
                j.exist = which(complete.cases(X[i, ]))
                t.new[i] = sum(X.old[i, j.exist] * w.new[j.exist])
            }
            for (j in 1:p) {
                i.exist = intersect(which(complete.cases(X[, j])), which(complete.cases(t.new)))
                p.new[j] = sum(X.old[i.exist, j] * t.new[i.exist])/sum(t.new[i.exist]^2)
            }
            c.new = t(Y.old) %*% t.new/sum(t.new^2)
            u.new = Y.old/as.vector(c.new)
        }
        if (!na.miss) {
            w.old = t(X.old) %*% Y.old/sum(Y.old^2)
            w.new = w.old/sqrt(sum(w.old^2))
            t.new = X.old %*% w.new                	#X-scores,  t
            p.new = t(X.old) %*% t.new/sum(t.new^2)	 
            c.new = t(Y.old) %*% t.new/sum(t.new^2) #Y-weights, c
            u.new = Y.old/as.vector(c.new)        	#an updated set of Y-scores
            if (crosval) {
                RSS[h + 1] = sum((Y.old - t.new %*% c.new)^2)
                press = rep(0, 10)
                for (i in 1:10) {
                  aux = segments[[i]]
                  Xy.aux = t(X.old[-aux, ]) %*% Y.old[-aux]
                  wh.si = Xy.aux %*% sqrt(solve(t(Xy.aux) %*% Xy.aux))
                  th.si = X.old[-aux, ] %*% wh.si
                  ch.si = t(Y.old[-aux]) %*% th.si %*% solve(t(th.si) %*% th.si)
                  ch.si = as.vector(ch.si)
                  Yhat.si = ch.si * X.old[aux, ] %*% wh.si
                  press[i] = sum((Y.old[aux] - Yhat.si)^2)
                }
                PRESS[h] = sum(press)
                Q2[h] = 1 - PRESS[h]/RSS[h]
            }
        }
        Y.old = Y.old - (t.new %*% c.new)		#remove the present component from X and Y
        X.old = X.old - (t.new %*% t(p.new))
        Th[, h] = t.new
        Ph[, h] = p.new
        Wh[, h] = w.new
        Uh[, h] = u.new
        ch[h] = c.new
        Hot[, h] = (n/(n - 1)) * t.new^2/(sum(t.new^2)/(n - 1))
        hlim[h] = qf(0.95, h, n - h) * (h * (n^2 - 1))/(n * (n - h))
        if (is.null(comps) && crosval) {
            if (Q2[h] < limq2 || h == nc) break
        } else {
            if (h == nc) break
        }
        h = h + 1
    }
    if (crosval) {
        q2cum = rep(NA, h)
        for (k in 1:h) q2cum[k] = prod(PRESS[1:k])/prod(RSS[1:k])
        Q2cum = 1 - q2cum
        Q2cv = cbind(PRESS[1:h], RSS[1:h], Q2[1:h], rep(limq2, h), Q2cum)
        dimnames(Q2cv) = list(1:h, c("PRESS", "RSS", "Q2", "LimQ2", "Q2cum"))
        if (is.null(comps)) h = h - 1
    }
    if (!crosval) Q2cv = NULL
    Th = Th[, 1:h]
    Ph = Ph[, 1:h]
    Wh = Wh[, 1:h]
    Uh = Uh[, 1:h]
    ch = ch[1:h]
    Ws = Wh %*% solve(t(Ph) %*% Wh)
    Bs = as.vector(Ws %*% ch)
    if (!na.miss) {
        Br = Bs * (rep(apply(Y, 2, sd), p)/apply(X, 2, sd))
        cte = as.vector(colMeans(Y) - Br %*% apply(X, 2, mean))
        y.hat = as.vector(X %*% Br + cte)
        cor.xyt = cor(cbind(Xx, y = Yy), Th)
    } else {
        mu.x <- attributes(Xx)$"scaled:center"
        sd.x <- attributes(Xx)$"scaled:scale"
        X.hat = Th %*% t(Ph) %*% diag(sd.x, p, p) + matrix(rep(mu.x, 
            each = n), n, p)
        Br = Bs * (rep(apply(Y, 2, sd), p)/sd.x)
        cte = as.vector(colMeans(response) - Br %*% mu.x)
        y.hat = as.vector(X.hat %*% Br + cte)
        cor.xyt = matrix(NA, p + 1, h)
        for (j in 1:p) {
            i.exist <- which(complete.cases(X[, j]))
            cor.xyt[j, ] = cor(Xx[i.exist, j], Th[i.exist, ])
        }
        cor.xyt[p + 1, ] = cor(Yy, Th)
    }
    resid = as.vector(Y - y.hat)
    R2 = as.vector(cor(Th, Yy))^2
    R2Xy = t(apply(cor.xyt^2, 1, cumsum))
    T2hot = rbind(hlim[1:h], t(apply(Hot[, 1:h], 1, cumsum)))
    dimnames(Wh) = list(colnames(X), paste0("w", 1:h))
    dimnames(Ws) = list(colnames(X), paste0("w*", 1:h))
    dimnames(Th) = list(rownames(X), paste0("t", 1:h))
    dimnames(Ph) = list(colnames(X),paste0("p", 1:h))
    dimnames(Uh) = list(rownames(Y), paste0("u", 1:h))
    names(ch) = paste0("c", 1:h)
    dimnames(T2hot) = list(c("T2", rownames(X)), paste0("H", 1:h))
    names(Bs) = colnames(X)
    names(Br) = colnames(X)
    names(resid) = rownames(Y)
    names(y.hat) = rownames(Y)
    names(R2) = paste0("t", 1:h)
    colnames(R2Xy) = paste0("t", 1:h)
    dimnames(cor.xyt) = list(c(colnames(X), colnames(Y)), colnames(Th))
    res = list(
      x.scores = Th, x.loads = Ph,
      y.scores = Uh, y.loads = ch,
      cor.xyt = cor.xyt,
      raw.wgs = Wh, mod.wgs = Ws,
      std.coefs = Bs, reg.coefs = c(Intercept = cte, Br),
      R2 = R2,
      R2Xy = R2Xy,
      y.pred = y.hat,
      resid = resid,
      T2 = T2hot,
      Q2 = Q2cv,
      y = response
    )
    class(res) = "plsreg1"
    return(res)
}
## ----------------------------- END PLS --------------------------------
