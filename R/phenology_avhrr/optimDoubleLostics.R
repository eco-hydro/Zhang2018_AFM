## 获得一个生长周期的物候指数
get_Metric_TRS <- function(x, t = index(x), tout = t, trs=0.6){
  if (all(is.na(x))) return(NULL)
  fit.spline <- splinefit(x, t, tout)
  # trs2 <- PhenoTrs(fit.spline$predicted, approach = "White", trs = 0.5)
  # approach is Trs or White
  ans <- lapply(c(0.2, 0.5, 0.6), function(tr) PhenoTrs(fit.spline$predicted, approach = "White", tr)) %>%
    do.call(c,.)
  return(ans)
}

get_MetricsBatch <- function(x, t = index(x), tout = t){
  if (all(is.na(x))) return(NULL)
  
  fit.spline <- splinefit(x, t, tout)
  fit.beck <- FitDoubleLogBeck(x, t, tout)
  fit.elmore <- FitDoubleLogElmore(x, t, tout)
  # fit.klos <- FitDoubleLogKlLight(x, t, tout)
  
  trs2.spline <- PhenoTrs(fit.spline$predicted, approach = "White", trs = 0.2)
  trs5.spline <- PhenoTrs(fit.spline$predicted, approach = "White", trs = 0.5)
  der.spline <- PhenoDeriv(fit.spline$predicted)
  
  # kl2 <- PhenoKl(fit.klos$params, fit.klos)
  # beck2 <- PhenoKl(fit.beck$params, fit.beck)
  # elmore2 <- PhenoKl(fit.elmore$params, fit.elmore)
  # fit_list <- list(spl = fit.spline, kl = fit.klos, beck = fit.beck, elmr = fit.elmore) %>% lapply(`[[`, "predicted")
  fit_list <- list(spl = fit.spline, 
                   beck = fit.beck, 
                   elmr = fit.elmore) %>% lapply(`[[`, "predicted")
  ## ------------------------------------------------
  ## 1. script model 
  ## Zhang XiaoYang, 2003 RSE
  # zhang.spline <- Phenokl.kong(fit.spline$predicted)
  # zhang.kl <- Phenokl.kong(fit.klos$predicted)
  # zhang.beck <- Phenokl.kong(fit.beck$predicted)
  # zhang.elmore <- Phenokl.kong(fit.elmore$predicted)
  # 
  # gu.spline <- PhenoGu.kong(fit.spline$predicted)[1:4]
  # gu.kl <- PhenoGu.kong(fit.klos$predicted)[1:4]
  # gu.beck <- PhenoGu.kong(fit.beck$predicted)[1:4]
  # gu.elmore <- PhenoGu.kong(fit.elmore$predicted)[1:4]
  ## ------------------------------------------------
  ## 2. batch model
  zhang <- lapply(fit_list, Phenokl.kong)
  gu <- lapply(fit_list, function(xpred) PhenoGu.kong(xpred)[1:4])
  return(list(trs2=trs2.spline, trs5=trs5.spline, der=der.spline, zhang = zhang, gu = gu))
  
}
#' @param outdir save calculated phenology in oudir txts
getYears_phenology <- function(xlist, dates_doy, outdir, PhenoFUN = get_MetricsBatch, ...){
  metrics.year <- list()
  nyear <- 32
  xx <- xlist$xx
  for (j in 1:nyear) {
    # fprintf("\tj=%dth------------\n", j)
    t <- dates_avhrr.doy[j,]
    tout <- t[1]:t[length(t)]
    xs <- zoo(xx[j,], t)
    # xl <- spline(doy, xx[i, ], xout = tout)$y %>% zoo(., tout)
    # a <- FitDoubleLogElmore(xs, t,tout)
    # b <- FitDoubleLogBeck(xs, t, tout)
    # be caution that trycatch warning can also lead to function stoped
    metrics.year[[j]] <- tryCatch(
      PhenoFUN(xs, t, tout, ...),
      warning = function(w) {
        cat(sprintf("[%d]WARNING %s", xlist$i, w))
        suppressWarnings(PhenoFUN(xs, t, tout, ...))
      },
      error = function(e) sprintf("[%d]ERROR %s", xlist$i, e)
    )
  }
  #20161122 names(metrics.year) <- 1982:2013
  names(metrics.year) <- seq(1982, by=1, length.out = length(metrics.year))#in order the last is NULL
  Id_NULL <- which(sapply(metrics.year, length) == 0)
  if (length(Id_NULL)) metrics.year %<>% .[-Id_NULL]
  if (length(metrics.year)){
    df <- ldply(metrics.year, function(x) round(unlist(x)), .id = "year")
    fname <- sprintf("%spheno_avhrr_%05d.txt", outdir, xlist$i)
    write.table(df,fname,col.names = T,row.names = F,sep = "\t")
  }else{
    warning(sprintf("[%d] Phenological metrics is null!\n", i))
  }
  # return(df)
}

## optim for double logistic functions
splinefit <- function(x, t = index(x), tout = t, plot = FALSE, ...){
	xpred.out <- spline(t, x, xout = tout)$y %>% zoo(., tout)
	return(list(predicted = xpred.out, params = NULL, formula = NULL))
}

optim_pheno <- function(x, t, tout, prior, FUN, formula, parnames, method){
    # FUN <- doubleLog.beck
    # pars.lst <- alply(prior, 1, optimx,
    #                   fn=.error, fun = FUN, x = x, t = t,
    #                   control = list(maxit = 1000, all.methods = TRUE, dowarn = FALSE))
    # pars.sort <- lapply(pars.lst, function(pars) pars[with(pars, order(convcode, value, xtimes)), ])
    # print(pars.sort)
    
    opt.df <- adply(prior, 1, optimx, .id=NULL, 
                    fn=.error, fun = FUN, x = x, t = t, method = method, 
                    control = list(maxit = 1000, dowarn = FALSE))
    best <- which.min(opt.df$value)
    
    # test for convergence if maximum iterations where reached - restart from best with more iterations
    if (opt.df$convcode[best] == 1) {
      # repeat with more maximum iterations if it did not converge
      par <- opt.df[best, parnames] %>% unlist #best par generally
      opt <- optimx(par, .error, fun = FUN, x = x, t = t, method = method, 
                    control = list(maxit = 1000))
    } else if (opt.df$convcode[best] == 0) {
      opt <- opt.df[best, ]
    }
    # return NA in case of no convergence
    if (opt$convcode != 0) {
      par <- numeric(length(parnames))*NA
      xpred <- rep(NA, length(tout))
    } else {
      par <- opt[, parnames] %>% as.numeric()
      xpred <- FUN(par, tout)
    }
    xpred.out <- zoo(xpred, order.by = tout)
    names(par) <- parnames
    return(list(predicted = xpred.out, params = par, formula = formula))
}

FitDoubleLogBeck <- function(x, t = index(x), tout = t, plot = FALSE, ...) {
    if (any(is.na(x))) 
        stop("NA in the time series are not allowed: fill them with e.g. na.approx()")
    
    n <- length(x)
    avg <- mean(x, na.rm = TRUE)
    mx <- max(x, na.rm = TRUE)
    mn <- min(x, na.rm = TRUE)
    ampl <- mx - mn

    doy <- quantile(t, c(0.25, 0.75), na.rm = TRUE)
    formula <- expression(mn + (mx - mn) * (1/(1 + exp(-rsp * (t - sos))) + 1/(1 + exp(rau * (t - eos)))))
    parnames <- c("mn", "mx", "sos", "rsp", "eos", "rau")
    prior <- rbind(
        c(mn, mx, doy[1], 0.5, doy[2], 0.5), 
        c(mn, mx, doy[2], 0.5, doy[1], 0.5), 
        c(mn - ampl/2, mx + ampl/2, doy[1], 0.5, doy[2], 0.5), 
        c(mn - ampl/2, mx + ampl/2, doy[2], 0.5, doy[1], 0.5)) %>% set_colnames(parnames)
    
    FUN <- doubleLog.beck
    RESULT <- optim_pheno(x, t, tout, prior, FUN, formula, parnames, method = "nlminb")
    return(RESULT)#list(predicted = xpred.out, params = par, formula = formula)
}

FitDoubleLogElmore <- function(x, t = index(x), tout = t, return.par = FALSE, plot = FALSE, ...) {
  if (any(is.na(x))) 
    stop("NA in the time series are not allowed: fill them with e.g. na.approx()")

  n <- length(na.omit(x))
  if (n < 7) {
    if (return.par)  return(rep(NA, 7))
    if (!return.par) return(rep(NA, length(tout)))
  }
  n <- length(x)
  avg <- mean(x, na.rm = TRUE)
  mx <- max(x, na.rm = TRUE)
  mn <- min(x, na.rm = TRUE)
  ampl <- mx - mn
  doy <- quantile(t, c(0.25, 0.75), na.rm = TRUE)
  
  formula <- expression(m1 + (m2 - m7 * t) * ((1/(1 + exp(((m3/m4) - t)/(1/m4)))) - (1/(1 + exp(((m5/m6) - t)/(1/m6))))))
  parnames <- paste0("m", 1:7)
  prior <- rbind(
    c(mn, mx - mn, 200, 1.5, 300, 1.5, 0.002), 
    c(mn, mx - mn, 100, 0.5, 200, 0.9, 0.002), 
    c(mn, mx - mn, 50, 0.5, 300, 1.2, 0.05), 
    c(mn, mx - mn, 300, 2, 350, 2.5, 0.05)) %>% set_colnames(parnames)
  
  FUN <- doubleLog.elmore
  RESULT <- optim_pheno(x, t, tout, prior, FUN, formula, parnames, method = "nlminb")
  return(RESULT)#list(predicted = xpred.out, params = par, formula = formula)
}

FitDoubleLogZhang <- function(x, t = index(x), tout = t, ..., method = "BFGS") {
  if (any(is.na(x))) 
    stop("NA in the time series are not allowed: fill them with e.g. na.approx()")
  n <- length(x)
  avg <- mean(x, na.rm = TRUE)
  mx <- max(x, na.rm = TRUE)
  mn <- min(x, na.rm = TRUE)
  ampl <- mx - mn
  doy.max <- which.max(x)
  doy <- quantile(t, c(0.25, 0.75), na.rm = TRUE)
  
  t1 <- doy[1] + 0.5 * (doy.max - doy[1])
  t2 <- doy.max + 0.5 * (doy[2] - doy.max)
  b1 <- 10  #ok
  b2 <- 10  #ok
  c1 <- 1  # ok 
  c2 <- 1
  
  # parnames <- c("y0", )
  prior <- rbind(
    c(y0, a1, a2, t01, t02, b1, b2, c1, c2), 
    c(y0, a1, a2, t01, t02, b1, b2, 1.2, c2), 
    c(y0, 0.05, 0.05, t01, t02, 0.5, b2, c1, c2), 
    c(y0, a1, a2, doy[1], t02, b1, b2, c1, c2), 
    c(y0, a1, a2, t01, doy[2], 5, 5, c1, c2))
  formula <- expression(y0 + (a1/(1 + exp(-(t - t1)/b1))^c1) - (a2/(1 + exp(-(t - t2)/b2))^c2))
}


FitDoubleLogKlLight <- function(x, t = index(x), tout = t, ..., method = "BFGS") {
    if (any(is.na(x))) 
        stop("NA in the time series are not allowed: fill them with e.g. na.approx()")

    n <- length(x)
    avg <- mean(x, na.rm = TRUE)
    mx <- max(x, na.rm = TRUE)
    mn <- min(x, na.rm = TRUE)
    ampl <- mx - mn
    
    doy <- quantile(t, c(0.25, 0.75), na.rm = TRUE)
    a1 <- 0  #ok
    a2 <- 0  #ok
    b1 <- mn  #ok
    b2 <- 0  #ok
    c <- 0.2 * max(x)  # ok 
    ## very slightly smoothed spline to get reliable maximum
    # tmp <- smooth.spline(x, df = 0.5 * length(x))#, find error: 20161104, fix tomorrow
    # doy.max <- which.max(tmp$y)
    doy.max <- which.max(x)#wait for test, kongdd
    
    B1 <- 4/(doy.max - doy[1])
    B2 <- 3.2/(doy[2] - doy.max)
    m1 <- doy[1] + 0.5 * (doy.max - doy[1])
    m2 <- doy.max + 0.5 * (doy[2] - doy.max)
    m1.bis <- doy[1]
    m2.bis <- doy[2]
    q1 <- 0.5  #ok
    q2 <- 0.5  #ok
    v1 <- 2  # ok
    v2 <- 2  # ok
    
    prior <- rbind(
        c(a1, a2, b1, b2, c, B1, B2, m1, m2, q1, q2, v1, v2), 
        c(a1, a2, b1, 0.01, 0, B1, B2, m1, m2.bis, q1, 1, v1, 4), 
        c(a1, a2, b1, b2, c, B1, B2, m1.bis, m2, q1, q2, v1, v2), 
        c(a1, a2, b1, b2, c, B1, B2, m1, m2.bis, q1, q2, v1, v2), 
        c(a1, a2, b1, b2, c, B1, B2, m1.bis, m2, q1, q2, v1, v2))
    opt.l <- apply(prior, 1, optim, .error, fun = doubleLog.klos, 
                   x=x, t=t, 
                   method = method, control = list(maxit = 1000), hessian = hessian)  # fit from different prior values
    opt.df <- cbind(cost = unlist(llply(opt.l, function(opt) opt$value)), 
                    convergence = unlist(llply(opt.l, function(opt) opt$convergence)), 
                    ldply(opt.l, function(opt) opt$par))
    pos <- which(opt.df$V6 < 0 | opt.df$V7 < 0)
    if (length(pos) != 0) {
        opt.df <- opt.df[-pos, ]
        opt.l <- opt.l[-pos]
    }
    best <- which.min(opt.df$cost)
    if (opt.df$convergence[best] == 1) {
        # if maximum iterations where reached - restart from best with more iterations
        opt <- opt.l[[best]]
        opt <- optim(opt.l[[best]]$par, .error, fun = doubleLog.klos, x = x, t = t, method = method, control = list(maxit = 700), hessian = hessian)
        prior <- rbind(prior, opt$par)
        xpred <- doubleLog.klos(opt$par, t)
    } else if (opt.df$convergence[best] == 0) {
        opt <- opt.l[[best]]
        prior <- rbind(prior, opt$par)
        xpred <- doubleLog.klos(opt$par, t)  ## perche questo restituisce nan?
    }
    if (opt$convergence != 0) {
        opt$par[] <- NA
        xpred <- rep(NA, length(tout))
    } else {
        xpred <- doubleLog.klos(opt$par, tout)
    }
    xpred.out <- zoo(xpred, order.by = tout)
    names(opt$par) <- c("a1", "a2", "b1", "b2", "c", "B1", "B2", "m1", "m2", "q1", "q2", "v1", "v2")
    
    if (hessian) {
        opt.new <- optim(opt$par, .error, fun = doubleLog.klos, x = x, t = t, method = method, hessian = TRUE)
        vc <- .qr.solve(opt$hessian)
        npar <- nrow(vc)
        s2 <- opt.df$cost[best]^2/(n - npar)
        std.errors <- sqrt(diag(vc) * s2)  # standard errors
    }
    
    fit.formula <- expression((a1 * t + b1) + (a2 * t^2 + b2 * t + c) * (1/(1 + q1 * exp(-B1 * (t - m1)))^v1 - 1/(1 + 
        q2 * exp(-B2 * (t - m2)))^v2))
    output <- list(predicted = xpred.out, params = opt$par, formula = fit.formula)
    if (hessian) 
        output <- list(predicted = xpred.out, params = opt$par, formula = fit.formula, stdError = std.errors)
    return(output)
}


doubleLog.gu <- function(par, t) {
    y0 = par[1]
    a1 <- par[2]
    a2 <- par[3]
    t01 <- par[4]
    t02 <- par[5]
    b1 <- par[6]
    b2 <- par[7]
    c1 <- par[8]
    c2 <- par[9]
    xpred <- y0 + (a1/(1 + exp(-(t - t01)/b1))^c1) - (a2/(1 + exp(-(t - t02)/b2))^c2)
    # xpred <- (a1*t + b1) + (a2*t^2 + b2*t + c)*(1/(1+q1*exp(-B1*(t-m1)))^v1 - 1/(1+q2*exp(-B2*(t-m2)))^v2)
    return(xpred)
}

doubleLog.klos <- function(par, t) {
    a1 <- par[1]
    a2 <- par[2]
    b1 <- par[3]
    b2 <- par[4]
    c <- par[5]
    B1 <- par[6]
    B2 <- par[7]
    m1 <- par[8]
    m2 <- par[9]
    q1 <- par[10]
    q2 <- par[11]
    v1 <- par[12]
    v2 <- par[13]
    xpred <- (a1 * t + b1) + (a2 * t^2 + b2 * t + c) * (1/(1 + q1 * exp(-B1 * (t - m1)))^v1 - 1/(1 + q2 * exp(-B2 * 
        (t - m2)))^v2)
    return(xpred)
}

doubleLog.beck <- function(par, t) {
    mn <- par[1]
    mx <- par[2]
    sos <- par[3]
    rsp <- par[4]
    eos <- par[5]
    rau <- par[6]
    
    # try(if (eos < sos) return(rep(99, length(t))), silent = T)
    xpred <- mn + (mx - mn) * (1/(1 + exp(-rsp * (t - sos))) + 1/(1 + exp(rau * (t - eos))))
    return(xpred)
}

doubleLog.elmore <- function(par, t) {
    m1 <- par[1]
    m2 <- par[2]
    m3 <- par[3]
    m4 <- par[4]
    m5 <- par[5]
    m6 <- par[6]
    m7 <- par[7]
    m3l <- m3/m4
    m4l <- 1/m4
    m5l <- m5/m6
    m6l <- 1/m6
    # xpred <- m1 + (m2 - m7 * t) * ( (1/(1 + exp(((m3/m4) - t)/(1/m4)))) - (1/(1 + exp(((m5/m6) - t)/(1/m6)))) )
    xpred <- m1 + (m2 - m7 * t) * ((1/(1 + exp((m3l - t)/m4l))) - (1/(1 + exp((m5l - t)/m6l))))
    return(xpred)
}
# only fit part of growing season NDVI, before or after peak NDVI
doubleLog.zhang2 <- function(par, t){
    c <- par[1]
    a <- par[2]
    b <- par[3]
    d <- par[4]
    xpred <- c/(1 + exp(a + b * t)) + d
    return(xpred)
}
# fit whole growing season NDVI
doubleLog.zhang <- function(par, t){
    t0 <- par[1]
    c1 <- par[2]
    a1 <- par[3]
    b1 <- par[4]
    d1 <- par[5]
    c2 <- par[6]
    a2 <- par[7]
    b2 <- par[8]
    d2 <- par[9]
    xpred1 <- c1/(1 + exp(a1 + b1 * t[t <= t0])) + d1
    xpred2 <- c2/(1 + exp(a2 + b2 * t[t > t0])) + d2
    xpred <- c(xpred1, xpred2)
    return(xpred)
}

.error <- function(par, x, t, 
    fun = c(doubleLog.elmore, doubleLog.beck, doubleLog.klos, doubleLog.gu), 
    weights = NULL){
    # FUN <- match.fun(fun)
    if (any(is.infinite(par))) return(99999)
    xpred <- fun(par, t = t)
    ifelse(is.null(weights), 
        sse <- sum((xpred - x)^2, na.rm = TRUE),
        sse <- sum((xpred - x)^2 * weights, na.rm=TRUE))
    return(sse)
}

.qr.solve <- function(a, b, tol = 1e-07, LAPACK = TRUE) {
    if (!is.qr(a)) a <- qr(a, tol = tol, LAPACK = LAPACK)
    nc <- ncol(a$qr)
    nr <- nrow(a$qr)
    if (a$rank != min(nc, nr)) stop("singular matrix 'a' in solve")
    if (missing(b)) {
        if (nc != nr) stop("only square matrices can be inverted")
        b <- diag(1, nc)
    }
    res <- qr.coef(a, b)
    res[is.na(res)] <- 0
    res
}
