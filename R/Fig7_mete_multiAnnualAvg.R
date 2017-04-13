rm(list = ls())
source('R/mainfunc/main.R', encoding = 'UTF-8')
source('R/mainfunc/main_PLS.R', encoding = 'UTF-8')

# load("data/phenology_list20170105.rda")
# load("../TP_phenoAgr/agrmet25_TP.rda")
load("../data/00basemap.rda")

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

x <- dataIn.trim$Prec
y <- array(x, dim = c(nrow(x), 12, ncol(x)/12))

a <- aaply(y, c(1, 3), sum, na.rm = T, .progress = "text")

library(data.table)
vars <- fread("data/TP_mete_multiAnnual.txt")[-Id_del, ] %>% as.data.frame()

gridClip@data <- vars[, 1, drop = FALSE]

source('R/mainfunc/main_PLOT.R', encoding = 'UTF-8')



rain = colorRampPalette(colors = RColorBrewer::brewer.pal())

brks_list <- list(
  precp = c(seq(0, 700, 100), 2000),
  srad = c(100, seq(2500, 3000, 100), 3200),
  Tmin = c(-20, -10, -7, -5, -2, 0, 5),
  Tmax = c(-10, 0, 4, 8, 10, 16),
  Tavg = c(-15, -5, 0, 5, 16)
)
# xx <- mapply(cut, x = vars, breaks = brks_list, include.lowest = F, ordered_result = T) %>% data.frame()
xx <- mapply(findInterval, x = vars, vec = brks_list) %>% data.frame()
myplot(2)


Names <- c("Precp", "Srad", "Tmin", "Tmax", "Tavg")
myplot <- function(i){
  gridClip@data <- xx[, i, drop = F]
  breaks <- brks_list[[i]]
  
  cols <- rainbow(length(breaks) - 1, start = 0.7, end = 0.1)
  if (i == 1) cols %<>% rev
  
  colors <- get_color(breaks)# %>% list2env(environment())
  # cols <- colors$cols
  colorkey <- colors$colorkey
  if (i > 2) cols <- rev(colors$cols)
  
  sp_layout = if (i < 3) list(sp_poly) else list(sp_cluster)
  spplot(gridClip, 
         at = seq(length(breaks)) - 0.5,
         col.regions = cols, 
         colorkey = colorkey,
         sp.layout = sp_layout,
         xlim = xlim, ylim = ylim, 
         par.settings = list(layout.heights = list(top.padding = 0, bottom.padding = 0), 
                             layout.widths = list(left.padding = 0, right.padding = 0, axis.left = 0, axis.right = 0)), 
         panel = function (x, y, z, subscripts, ..., sp.layout) 
         {
           sppanel(list(sp.layout), panel.number(), first = TRUE)
           
           if (i > 2) {
             panel.levelplot.raster(x, y, z, subscripts, ...)
             panel.contourplot(x, y, z, subscripts, contour = TRUE, labels = F, col = "red", 
                               at = c(-Inf, 0, Inf), col.regions = "transparent")
           }else{
             panel.contourplot(x, y, z, subscripts, contour = TRUE, labels = F, col = "black", ...)
           }
           
           sppanel(list(sp.layout), panel.number(), first = FALSE)
           #not significant subplot
           
           # panel.text(76, 40.4, paste0(zcol[i]), #english name: New_names[i])
           # fontfamily = "Times", cex = 1.3, font = 2, adj = 0)#"(",letters[i], ") ", 
           panel.text(76.6, 40.2, Names[i], fontfamily = "Times", cex = 1.3, font = 2, adj = 0)
           panel.addbarchart(z, subscripts, cols, showCluster = F, A = 9, by = 0.9, ntick = 2, origin.x = 76,...)
         })
}


library(gridExtra)

ps <- llply(seq_along(xx), myplot)
p <- do.call(grid.arrange, c(ps, nrow = 2, ncol = 3))
p

## 计算逐月平均

lapply(c("layout.heights", "layout.widths"), trellis.par.get) %>% str()

