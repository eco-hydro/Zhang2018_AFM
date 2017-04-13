rm(list=ls())
source('R/mainfunc/main.R', encoding = 'UTF-8')
source('R/mainfunc/main_VSI.R', encoding = 'UTF-8')

load("data/Vegeteation_Sensitivity.rda")
load("../data/00basemap.rda")

VSI <- VSI[-Id_del]

## try to visual PCALM coef directly
varCoef <- lapply(VSI, `[[`, "varCoef")
varCoef.list   <- llply(1:12, function(i) ldply(varCoef,   function(x) x[i, -1], .id = NULL), .progress = "text") %>% set_names(1:12)
pVal <- ldply(VSI,  function(x) as.numeric(as.matrix(x$pVal)), .id = NULL)

x <- do.call(cbind.data.frame, varCoef.list[4:10]) %>% set_colnames(paste0("M", colnames(.)))

opts <- llply(varCoef.list, function(x) {
  meteVars <- c("Tmax", "Tavg", "Tmin", "Prec", "Srad")
  opt <- apply(x, 1, function(x) which.max(abs(x))) 
  opt[sapply(opt, length) == 0] <- NA
  opt <- factor(names(unlist(opt)), meteVars)
  return(opt)
}, .progress = "text") %>% do.call(cbind.data.frame, .)

gridClip@data <- opts

CairoPNG(filename = "test01.png", 8.50, 11, dpi = 250, pointsize = 12, units = "in")
# spplot(gridClip, as.table = T, layout = c(5, 12))
print(p)
dev.off()

writeGDAL(gridClip, "dominate")

x <- varCoef.list$`1`[, -1]

## test the relation of sdDetVec and meanVec