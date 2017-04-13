PhenoTrs <- function(x, approach = c("White", "Trs"), trs = 0.5, min.mean = 0.1, ...) {
    if (all(is.na(x))) return(c(sos = NA, eos = NA, los = NA))
    
    # get statistical values
    n <- index(x)[length(x)]
    avg <- mean(x, na.rm = TRUE)
    x2 <- na.omit(x)
    avg2 <- mean(x2[x2 > min.mean], na.rm = TRUE)
    peak <- max(x, na.rm = TRUE)
    mn <- min(x, na.rm = TRUE)
    ampl <- peak - mn
    
    # get peak of season position
    pop <- median(index(x)[which(x == max(x, na.rm = TRUE))])
    # select (or scale) values and thresholds for different methods
    approach <- approach[1]
    if (approach == "White") {
        # scale annual time series to 0-1
        ratio <- (x - mn)/ampl
        # trs <- 0.5
        trs.low <- trs - 0.1
        trs.up <- trs + 0.1
    }
    if (approach == "Trs") {
        ratio <- x
        a <- diff(range(ratio, na.rm = TRUE)) * 0.1
        trs.low <- trs - a
        trs.up <- trs + a
    }
    # identify greenup or dormancy period
    .Greenup <- function(x, ...) {
        ratio.deriv <- c(NA, diff(x))
        greenup <- rep(NA, length(x))
        greenup[ratio.deriv > 0] <- TRUE
        greenup[ratio.deriv < 0] <- FALSE
        return(greenup)
    }
    greenup <- .Greenup(ratio)
    
    ## 添加限制因素，区分上班年下半年
    # select time where SOS and EOS are located (around trs value)
    bool <- ratio >= trs.low & ratio <= trs.up
    
    # get SOS, EOS, LOS
    soseos <- index(x)
    half <- 210#上半年与下半年的区分
    
    # fixed 2017-01-04, according to TP phenology property
    # sos <- round(median(soseos[greenup & bool], na.rm = TRUE))
    # eos <- round(median(soseos[!greenup & bool], na.rm = TRUE))

    sos <- round(median(soseos[greenup & bool & soseos < half], na.rm = TRUE))
    eos <- round(median(soseos[!greenup & bool & soseos > half], na.rm = TRUE))
    los <- eos - sos#los < 0 indicate that error
    # los[los < 0] <- n + (eos[los < 0] - sos[los < 0])
    
    # get MGS, MSP, MAU
    # mgs <- mean(x[ratio > trs], na.rm = TRUE)
    # msp <- mau <- NA
    # if (!is.na(sos)) {
    #     id <- (sos - 10):(sos + 10)
    #     id <- id[(id > 0) & (id < n)]
    #     msp <- mean(x[which(index(x) %in% id == TRUE)], na.rm = TRUE)
    # }
    # if (!is.na(eos)) {
    #     id <- (eos - 10):(eos + 10)
    #     id <- id[(id > 0) & (id < n)]
    #     mau <- mean(x[which(index(x) %in% id == TRUE)], na.rm = TRUE)
    # }
    # metrics <- c(sos = sos, eos = eos, los = los, pop = pop, mgs = mgs, rsp = NA, rau = NA, peak = peak, msp = msp, mau = mau)
    metrics <- c(sos = sos, eos = eos, los = los)
    # if (plot) { if (approach == 'White') PlotPhenCycle(x, metrics=metrics, trs=trs, ...)  if (approach == 'Trs')
    # PlotPhenCycle(ratio, metrics=metrics, trs=trs, ...)  }
    return(metrics)
    ### The function returns a vector with SOS, EOS, LOS, POP, MGS, rsp, rau, PEAK, MSP and MAU. }
}

PhenoDeriv <- function(x, ...) {
    if (all(is.na(x))) 
        return(c(sos = NA, eos = NA, los = NA, pop = NA, mgs = NA, rsp = NA, rau = NA, peak = NA, msp = NA, mau = NA))
    n <- index(x)[length(x)]
    avg <- mean(x, na.rm = TRUE)
    x2 <- na.omit(x)
    # avg2 <- mean(x2[x2 > min.mean], na.rm=TRUE)
    peak <- max(x, na.rm = TRUE)
    mn <- min(x, na.rm = TRUE)
    ampl <- peak - mn
    
    # get peak of season position
    pop <- median(index(x)[which(x == max(x, na.rm = TRUE))])
    # return NA if amplitude is too low or time series has too many NA values if (!calc.pheno) { if (avg < min.mean) { #
    # return for all metrics NA if mean is too low return(c(sos=NA, eos=NA, los=NA, pop=NA, mgs=NA, rsp=NA, rau=NA,
    # peak=NA, msp=NA, mau=NA)) } else { # return at least annual average if annual mean > min.mean return(c(sos=NA,
    # eos=NA, los=NA, pop=pop, mgs=avg2, rsp=NA, rau=NA, peak=peak, msp=NA, mau=NA)) } }
    
    # calculate derivative
    xd <- c(NA, diff(x))
    
    # get SOS and EOS
    soseos <- index(x)
    rsp <- max(xd, na.rm = TRUE)
    rau <- min(xd, na.rm = TRUE)
    sos <- median(soseos[xd == rsp], na.rm = TRUE)
    eos <- median(soseos[xd == rau], na.rm = TRUE)
    los <- eos - sos
    los[los < 0] <- n + (eos[los < 0] - sos[los < 0])
    
    # get MGS
    if (sos < eos) {
        mgs <- mean(x[index(x) %in% sos:eos], na.rm = TRUE)
    } else {
        cut1 <- as.vector(window(x, end = eos))
        cut2 <- as.vector(window(x, start = sos))
        mgs <- mean(c(cut1, cut2), na.rm = TRUE)
    }
    # get MSP, MAU
    msp <- mau <- NA
    if (!is.na(sos)) {
        id <- (sos - 10):(sos + 10)
        id <- id[(id > 0) & (id < n)]
        msp <- mean(x[which(index(x) %in% id == TRUE)], na.rm = TRUE)
    }
    if (!is.na(eos)) {
        id <- (eos - 10):(eos + 10)
        id <- id[(id > 0) & (id < n)]
        mau <- mean(x[which(index(x) %in% id == TRUE)], na.rm = TRUE)
    }
    # metrics <- c(sos = sos, eos = eos, los = los, pop = pop, mgs = mgs, rsp = rsp, rau = rau, peak = peak, msp = msp, 
    #     mau = mau)
    metrics <- c(sos = sos, eos = eos, los = los, pop = pop)
    # if (plot) { PlotPhenCycle(x, metrics=metrics, ...)  }
    return(metrics)
    ### The function returns a vector with SOS, EOS, LOS, POP, MGS, RSP, RAU, PEAK, MSP and MAU. }
}

PhenoGu <- function(x, fit, ...) {
    if (is.null(x)) {
        x <- fit$predicted
        spline.eq <- smooth.spline(x, df = length(x))
        der1 <- predict(spline.eq, d = 1)$y
        t <- index(x)
        days <- index(x)
        values <- as.vector(x)
    } else {
        names(x) <- names(fit$params)
        retrieved.formula <- fit$formula
        days <- index(fit$predicted)
        t <- index(fit$predicted)
        values <- as.vector(fit$predicted)
        D1 <- D(retrieved.formula, "t")
        # D2 <- D(D1, "t")
        ## e1 <- parent.frame()
        der1 <- eval(D1, envir = as.list(x))
    }
    # if (length(which(is.na(der1)==TRUE))!=0 | length(which(is.infinite(der1)==TRUE))!=0) { der1[is.na(der1)] <- 0
    # der1[is.infinite(der1)] <- 0 warning('Check your fitting because the first derivative contains NA or infinite
    # values \n They were set at 0!') }
    if (any(is.na(der1))) der1 %<>% na.approx(rule = 2)
    if (any(is.infinite(der1))){
        metrics <- rep(NA, 9)
    } else {
        ## extract parameters,parameters <- fit$params
        # get statistical values
        prr <- max(der1, na.rm = T)
        ## peak recovery date
        prd <- days[which.max(der1)]
        ## peak senescence rate
        psr <- min(der1, na.rm = T)
        ## peak senescence date
        psd <- days[which.min(der1)]
        ## gcc @ prd
        y.prd <- values[which(days == prd)]
        ## time peak recovery rate tPRD <- parameters['t01'] + parameters['b1']*log(parameters['c1']) tPSD <-
        ## parameters['t02'] + parameters['b2']*log(parameters['c2']) ## recovery line
        rl.y <- prr * days + y.prd - prr * prd
        rl.eq <- lm(rl.y ~ days)
        ## gcc @ psd
        y.psd <- values[which(days == psd)]
        ## senenscence line
        sl.y <- psr * days + y.psd - psr * psd
        sl.eq <- lm(sl.y ~ days)
        baseline <- min(values, na.rm = T)
        maxline <- max(values, na.rm = T)
        ## upturn day is the intersection between rl and x axis
        UD <- (baseline - rl.eq$coefficients[1])/rl.eq$coefficients[2]
        ## recession day is the intersection between sl and x axis
        RD <- (baseline - sl.eq$coefficients[1])/sl.eq$coefficients[2]
        ## stabilization day, intersection between maxline and rl
        SD <- (maxline - rl.eq$coefficients[1])/rl.eq$coefficients[2]
        ## downturn day, intersection between maxline and sl
        DD <- (maxline - sl.eq$coefficients[1])/sl.eq$coefficients[2]
        ## subset data between SD and DD
        sub.time <- days[which(days >= SD & days <= DD)]
        sub.gcc <- values[which(days >= SD & days <= DD)]
        if (length(sub.time) > 3) {
            ## compute a linear fit
            plateau.lm <- lm(sub.gcc ~ sub.time)
            M <- matrix(c(coef(plateau.lm)[2], coef(sl.eq)[2], -1, -1), nrow = 2, ncol = 2)
            intercepts <- as.matrix(c(coef(plateau.lm)[1], coef(sl.eq)[1]))
            interception <- -solve(M) %*% intercepts
            DD <- interception[1, 1]
        }
        ## calculate area under the curve cut.x <- days[which(days>=UD & days<=RD)] cut.y <- offset.y[which(days>=UD &
        ## days<=RD)] the.fun <- function(t) {eval(retrieved.formula, envir=as.list(params))}
        plateau.slope <- try(plateau.lm$coefficients[2], silent = TRUE)
        plateau.intercept <- try(plateau.lm$coefficients[1], silent = TRUE)
        if (class(plateau.slope) == "try-error") {
            plateau.slope <- NA
            plateau.intercept <- NA
        }
        metrics <- c(UD, SD, DD, RD, maxline, baseline, prr, psr, plateau.slope)
        # if (length(which(diff(metrics[1:4])<0)!=0)) { metrics <- rep(NA, 9) warning('Threshold do not respect expected
        # timing: set to NA') }
    }
    names(metrics) <- c("UD", "SD", "DD", "RD", "maxline", "baseline", "prr", "psr", "plateau.slope")
    return(metrics)
}

PhenoGu.kong <- function(xpred, ...) {
  if (all(is.na(xpred))){
    metrics <- rep(NA, 9)
  } else {
    days <- t <- index(xpred)
    values <- as.vector(xpred)
    der1 <- c(NA, diff(values)); der1 %<>% na.approx(rule = 2)
    
    ## extract parameters, parameters <- fit$params
    # get statistical values
    prr <- max(der1, na.rm = T)
    ## peak recovery date
    prd <- days[which.max(der1)]
    ## peak senescence rate
    psr <- min(der1, na.rm = T)
    ## peak senescence date
    psd <- days[which.min(der1)]
    ## gcc @ prd
    y.prd <- values[which(days == prd)]
    ## time peak recovery rate tPRD <- parameters['t01'] + parameters['b1']*log(parameters['c1']) tPSD <-
    ## parameters['t02'] + parameters['b2']*log(parameters['c2']) ## recovery line
    rl.y <- prr * days + y.prd - prr * prd
    rl.eq <- lm(rl.y ~ days)
    ## gcc @ psd
    y.psd <- values[which(days == psd)]
    ## senenscence line
    sl.y <- psr * days + y.psd - psr * psd
    sl.eq <- lm(sl.y ~ days)
    baseline <- min(values, na.rm = T)
    maxline <- max(values, na.rm = T)
    ## upturn day is the intersection between rl and x axis
    UD <- (baseline - rl.eq$coefficients[1])/rl.eq$coefficients[2]
    ## recession day is the intersection between sl and x axis
    RD <- (baseline - sl.eq$coefficients[1])/sl.eq$coefficients[2]
    ## stabilization day, intersection between maxline and rl
    SD <- (maxline - rl.eq$coefficients[1])/rl.eq$coefficients[2]
    ## downturn day, intersection between maxline and sl
    DD <- (maxline - sl.eq$coefficients[1])/sl.eq$coefficients[2]
    ## subset data between SD and DD
    sub.time <- days[which(days >= SD & days <= DD)]
    sub.gcc <- values[which(days >= SD & days <= DD)]
    if (length(sub.time) > 3) {
      ## compute a linear fit
      plateau.lm <- lm(sub.gcc ~ sub.time)
      M <- matrix(c(coef(plateau.lm)[2], coef(sl.eq)[2], -1, -1), nrow = 2, ncol = 2)
      intercepts <- as.matrix(c(coef(plateau.lm)[1], coef(sl.eq)[1]))
      interception <- -solve(M) %*% intercepts
      DD <- interception[1, 1]
    }
    ## calculate area under the curve cut.x <- days[which(days>=UD & days<=RD)] cut.y <- offset.y[which(days>=UD &
    ## days<=RD)] the.fun <- function(t) {eval(retrieved.formula, envir=as.list(params))}
    plateau.slope <- try(plateau.lm$coefficients[2], silent = TRUE)
    plateau.intercept <- try(plateau.lm$coefficients[1], silent = TRUE)
    if (class(plateau.slope) == "try-error") {
      plateau.slope <- NA
      plateau.intercept <- NA
    }
    metrics <- c(UD, SD, DD, RD, maxline, baseline, prr, psr, plateau.slope)
    # if (length(which(diff(metrics[1:4])<0)!=0)) { metrics <- rep(NA, 9) warning('Threshold do not respect expected
    # timing: set to NA') }
  }
  names(metrics) <- c("UD", "SD", "DD", "RD", "maxline", "baseline", "prr", "psr", "plateau.slope")
  return(metrics)
}

PhenoKl <- function(x, fit, ...) {
    # x <- ifelse(uncert==TRUE, x$uncertainty, x$fit) x are the parameters condition to understand if we enter the
    # function with a Spline fitting or else.
    if (is.null(x)) {
        x <- fit$predicted
        spline.eq <- smooth.spline(x, df = length(x))
        der1 <- predict(spline.eq, d = 1)$y
        der2 <- predict(spline.eq, d = 2)$y
        t <- index(x)
        values <- as.vector(x)
    } else {
        names(x) <- names(fit$params)
        retrieved.formula <- fit$formula
        t <- index(fit$predicted)
        values <- as.vector(fit$predicted)
        D1 <- D(retrieved.formula, "t")
        D2 <- D(D1, "t")
        ## e1 <- parent.frame()
        der1 <- eval(D1, envir = as.list(x))
        der2 <- eval(D2, envir = as.list(x))
        ## in case for NA values
        if (any(is.na(der1))) der1 %<>% na.approx(rule = 2)
        if (any(is.na(der2))) der2 %<>% na.approx(rule = 2)
    }
    if (all(is.na(x))) 
        metrics <- rep(NA, 4) 
    else {
        k <- der2/(1 + der1^2)^(3/2)
        # if (length(which(is.na(k) == TRUE)) != 0 | length(which(is.infinite(k) == TRUE)) != 0) {
        if (any(is.na(k)) | any(is.infinite(k))) {
            metrics <- rep(NA, 4)
        } else {
            spline.k <- smooth.spline(k, df = 0.1 * length(k))
            der.k <- predict(spline.k, d = 1)$y
            ## find maxima of derivative of k ## split season
            half.season <- which.max(values) + 20
            increasing.k <- try(der.k[1:half.season])
            increasing.k.d <- try(t[1:half.season])
            decreasing.k <- try(der.k[half.season:length(k)])
            decreasing.k.d <- try(t[half.season:length(k)])
            check.list <- list(increasing.k, decreasing.k)
            classes <- lapply(check.list, class)
            if (length(which(classes == "try-error")) != 0) 
                metrics <- rep(NA, 4) else {
                ## subset before first min
                subset1 <- increasing.k[1:which.min(increasing.k)]
                subset1.d <- increasing.k.d[1:which.min(increasing.k)]
                p1 <- subset1.d[which.max(subset1)]
                ## subset between first min and mid season
                subset2 <- increasing.k[which.min(increasing.k):length(increasing.k)]
                subset2.d <- increasing.k.d[which.min(increasing.k):length(increasing.k)]
                p2 <- subset2.d[which.max(subset2)]
                ## subset between mid season and max
                subset3 <- decreasing.k[1:which.max(decreasing.k)]
                subset3.d <- decreasing.k.d[1:which.max(decreasing.k)]
                p3 <- subset3.d[which.min(subset3)]
                ## rest of the season
                subset4 <- decreasing.k[which.max(decreasing.k):length(decreasing.k)]
                subset4.d <- decreasing.k.d[which.max(decreasing.k):length(decreasing.k)]
                p4 <- subset4.d[which.min(subset4)]
                metrics <- c(p1, p2, p3, p4)
            }
        }
    }
    ## extract parameters
    names(metrics) <- c("Greenup", "Maturity", "Senescence", "Dormancy")
    return(metrics)
}

#' @description Extract phenology using curvature as Zhang Xioayang, 2003
#' @param xpred an zoo object returned by double logistic function fitting
Phenokl.kong <- function(xpred, ...) {
  if (all(is.na(xpred))) 
      metrics <- rep(NA, 4) 
  else {
      t <- index(xpred)
      values <- as.vector(xpred)
      der1 <- c(NA, diff(values)); der1 %<>% na.approx(rule = 2)
      der2 <- c(NA, diff(der1)); der2 %<>% na.approx(rule = 2)
      k <- der2/(1 + der1^2)^(3/2)

      if (any(is.na(k)) | any(is.infinite(k))) {
        metrics <- rep(NA, 4)
      } else {
        spline.k <- smooth.spline(k, df = 0.1 * length(k))
        der.k <- predict(spline.k, d = 1)$y
        ## find maxima of derivative of k ## split season
        half.season <- which.max(values) + 20
        increasing.k <- try(der.k[1:half.season])
        increasing.k.d <- try(t[1:half.season])
        decreasing.k <- try(der.k[half.season:length(k)])
        decreasing.k.d <- try(t[half.season:length(k)])
        check.list <- list(increasing.k, decreasing.k)
        classes <- lapply(check.list, class)
        if (length(which(classes == "try-error")) != 0) 
          metrics <- rep(NA, 4) else {
            ## subset before first min
            subset1 <- increasing.k[1:which.min(increasing.k)]
            subset1.d <- increasing.k.d[1:which.min(increasing.k)]
            p1 <- subset1.d[which.max(subset1)]
            ## subset between first min and mid season
            subset2 <- increasing.k[which.min(increasing.k):length(increasing.k)]
            subset2.d <- increasing.k.d[which.min(increasing.k):length(increasing.k)]
            p2 <- subset2.d[which.max(subset2)]
            ## subset between mid season and max
            subset3 <- decreasing.k[1:which.max(decreasing.k)]
            subset3.d <- decreasing.k.d[1:which.max(decreasing.k)]
            p3 <- subset3.d[which.min(subset3)]
            ## rest of the season
            subset4 <- decreasing.k[which.max(decreasing.k):length(decreasing.k)]
            subset4.d <- decreasing.k.d[which.max(decreasing.k):length(decreasing.k)]
            p4 <- subset4.d[which.min(subset4)]
            metrics <- c(p1, p2, p3, p4)
          }
      }
    }
  ## extract parameters
  names(metrics) <- c("Greenup", "Maturity", "Senescence", "Dormancy")
  return(metrics)
}