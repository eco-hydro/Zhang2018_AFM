Id_error <- read.table("Id_cal.txt", header = T)$x

NDVI_error <- NDVI_list[Id_error]
indir <- "E:/Phenology/AVHRR/"

result <- list()
for (i in seq_along(NDVI_error)) {
  cat(sprintf("[%d]--------------------------\n", i))
  xlist <- NDVI_error[[i]]
  
  metrics.year <- list()
  nyear <- 32
  xx <- xlist$xx
  plot(as.numeric(t(xx))[1:(24 * 6)], type = "b")
  ## warning lead to script stoped
  for (j in 1:nyear) {
    fprintf("\tj=%dth------------\n", j)
    t <- dates_avhrr.doy[j,]
    tout <- t[1]:t[length(t)]
    xs <- zoo(xx[j,], t)
    # xl <- spline(doy, xx[i, ], xout = tout)$y %>% zoo(., tout)
    # a <- FitDoubleLogElmore(xs, t,tout)
    # b <- FitDoubleLogBeck(xs, t, tout)
    
    # be caution that trycatch warning can also lead to function stoped
    metrics.year[[j]] <- tryCatch(get_MetricsBatch(xs, t, tout),
                                  warning = function(w) cat(sprintf("WARNING[%d] %s", xlist$i, w)), 
                                  error = function(e) cat(sprintf("ERROR[%d] %s", xlist$i, e)))
  }
  names(metrics.year) <- seq(1982, by=1, length.out = length(metrics.year))#in order the last is NULL
  Id_NULL <- which(sapply(metrics.year, length) == 0)
  if (length(Id_NULL)) metrics.year %<>% .[-Id_NULL]
  df <- ldply(metrics.year, function(x) round(unlist(x)), .id = "year")
  # df <- ldply(metrics.year, function(x) round(unlist(x)))
  fname <- sprintf("%spheno_avhrr_%d.txt", outdir, xlist$i)
  write.table(df,fname,col.names = T,row.names = F,sep = "\t")
}
# result[[i]] <- getYears_phenology(xlist, dates_avhrr.doy, indir)

tryCatch(get_MetricsBatch(xs, t, tout),
         # warning = function(w){ w; return(val)}, 
         warning = function(w) {
           cat(sprintf("WARNING[%d] %s", xlist$i, w))
           invokeRestart("muffleWarning")
           },
         error = function(e) cat(sprintf("ERROR[%d] %s", xlist$i, e)))

withCallingHandlers(e,
         # warning = function(w){ w; return(val)}, 
         # warning = function(w) cat(sprintf("WARNING[%d] %s", xlist$i, w)),
         error = function(e) cat(sprintf("ERROR[%d] %s", xlist$i, e)))