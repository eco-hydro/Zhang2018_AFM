library(lubridate)
library(sp)
library(maptools)
library(floodmap)
library(magrittr)
# library(rlist)
# library(zoo)
library(plyr)
library(data.table)
library(lattice)
library(latticeExtra)
# library(grid)
library(Cairo)
library(ggplot2)

library(foreach)
library(iterators)
source("E:/GitHub/statistics_Kong.R", encoding = "utf-8")

## ------------------------ GLOBAL VARIABLES -------------------------------
Id_del <- c(16035, 18443, 20925, 21407, 21618, 21732, 21827)

vars_origin <- c("TRS2", "TRS5", "DER", "Zhang.spl.Senescence", "Zhang.spl.Dormancy", 
          "Zhang.beck.Senescence", "Zhang.beck.Dormancy",
          "Zhang.elmr.Senescence", "Zhang.elmr.Dormancy",
          "Gu.spl.DD", "Gu.spl.RD", 
          "Gu.beck.DD", "Gu.beck.RD",
          "Gu.elmr.DD", "Gu.elmr.RD")
vars <- c("TRS5", "TRS6", "DER", "Zhang.spl.Senescence", "Zhang.spl.Dormancy", 
          "Zhang.beck.Senescence", "Zhang.beck.Dormancy",
          "Zhang.elmr.Senescence", "Zhang.elmr.Dormancy",
          "Gu.spl.DD", "Gu.spl.RD", 
          "Gu.beck.DD", "Gu.beck.RD",
          "Gu.elmr.DD", "Gu.elmr.RD")
# vars_order <- vars_new[c(4:15, 1:3)]
# vars <- vars_new
oyg.colors <- colorRampPalette(c("green4", "yellow",  "red"))
# ggplot for plsr model coef theme initial
mytheme <- theme_grey(base_size = 14, base_family = "Times") +
  theme(
    # legend.position = c(0.02, 0.02), legend.justification = c(0, 0), 
    # legend.position="bottom", 
    text = element_text(colour = "black", size = 14), 
    # axis.title.x = element_blank(),
    # plot.title = element_text(hjust = 0),
    axis.title = element_text(face = "bold", size = 14), 
    axis.text = element_text(colour = "black", size = 14), 
    strip.text = element_text(colour = "black", size = 14),
    # axis.text.x = element_text(angle = 0, size = 14), 
    # legend.margin = unit(0.2, "cm"), 
    legend.text = element_text(family = "Times", size = 14, face = "bold"), 
    legend.title = element_blank(), 
    # panel.grid.minor.x = element_blank(), 
    panel.grid.major = element_line(colour = "white"), 
    panel.grid.minor = element_line(size = 0.2))
theme_set(mytheme)
Sys.setlocale("LC_TIME", "english")#

## ------------------------------- GLOBAL FUNCTIONS --------------------------
## 对物候指标进行进一步处理, 首先可以忽略los, pop
#  计算的物候指标中，如果末期物候指标小于初期物候指标，以及小于零的物候指标设置成空值
#  x: data.table class
filter_rational.dt <- function(x){
  x[x <= 0 | x > 366] <- NA
  x_list <- lapply(method.group, function(i) x[, i, with = FALSE])
  
  ## can adapt data.table class
  res <- lapply(x_list, function(x_df){
    Id_na <- apply(x_df, 1, order) %>% apply(2, function(x) {
      x_unique <- unique(diff(x))
      (length(x_unique) != 1) || (x_unique != 1)#should be attention that `|` different from `||`
    }) %>% which
    x_df[Id_na, ] <- NA#if Id_na is NULL have no any negtive effects
    x_df#需要对trs2结果进行改进
  })#quickly return, x_list.trim <-
  do.call(cbind, setNames(res, NULL))#return data.frame
}
filter_rational <- function(x){
  x[x <= 0 | x > 366] <- NA
  x_list <- lapply(method.group, function(i) x[, i])
  
  ## can adapt data.table class
  res <- lapply(x_list, function(x_df){
    Id_na <- apply(x_df, 1, order) %>% apply(2, function(x) {
      x_unique <- unique(diff(x))
      (length(x_unique) != 1) || (x_unique != 1)#should be attention that `|` different from `||`
    }) %>% which
    x_df[Id_na, ] <- NA#if Id_na is NULL have no any negtive effects
    x_df#需要对trs2结果进行改进
  })#quickly return, x_list.trim <-
  do.call(cbind, setNames(res, NULL))#return data.frame
}

## ----------------------- PLOT FUNCTIONS --------------------------
# change factor levels for plot
changeLevels <- function(x, at){
  result <- data.frame(x)
  for (i in 1:ncol(x))
    result[, i] <- cut(x[, i], breaks = at, include.lowest = T)
  result#quickly return
}

pos <- data.frame(x = c(99, 92.5, 86.5, 98,   81,   102.6,  100,  86), 
                  y = c(29, 33,   31,   33.2, 34.1, 32.8, 36.7, 27))
pos$name <- letters[1:nrow(pos)]

# plot(sp_cluster[[2]], axes = T)
# points(pos$x, pos$y, pch = 3, col = "red"); text(pos$x, pos$y, pos$name)

# plot(seq(75, 105, length.out = 100), seq(20, 40, length.out = 100), type = "n"); text(pos$x, pos$y, pos$name)

panel.addaxis <- function(prec, origin.x = 76, origin.y = 26.5, tck = 0.2, 
                          dx = 0.4, A = 15, by = 0.6, col, box.width = by - 0.1, ntick = 3, 
                          text.cex = 1, ...){
  # dx <- 0.4;#为了使y axis和barchart分开
  # tck <- 0.2#tick length
  ypos = prec*A + origin.y 
  xpos <- seq(origin.x + dx, by = by, length.out = length(prec))
  panel.barchart(x = xpos, y = ypos, horizontal = F, origin = origin.y, 
                 reference = F, col = col, box.width = box.width,...)
  
  ymax <- ceiling(max(prec)*10)/10
  # ymax <- round(max(prec), 1)
  if (ymax <= 0.2) ntick <- 1 #else ntick <- 3
  tick <- pretty(c(0, ymax), ntick);
  if (ymax >= 0.5 && ymax <= 0.6) tick <- c(0, 0.3, 0.6)
  if (length(tick) > 3) {
    print(str(listk(tick, ymax, ntick)))
  }

  tick_pos <- tick*A + origin.y
  ## draw axis
  panel.lines(rep(origin.x , 2), c(origin.y, tick_pos[length(tick_pos)]), col = "black")
  panel.lines(c(origin.x, xpos[length(xpos)] + by/2), rep(origin.y, 2), col = "black")
  ## draw axis line ticks
  for (i in seq_along(tick_pos)){
    panel.lines(c(0, -tck) + origin.x, rep(tick_pos[i], 2), col = "black")
    panel.text(origin.x -tck, tick_pos[i], tick[i], fontfamily = "Times", cex = text.cex, adj = 1)
  }
}

panel.addbarchart <- function(z, subscripts, cols.avg, showCluster = T, ...){
  z <- z[subscripts];
  Id_remain <- which(!is.na(z))
  # origin.x <- 76; origin.y <- 26.5; tck = 0.2; dx = 0.4; by = 0.6; A = 15; box.width = 0.5
  if (showCluster){
    ## for every single cluster region
    info <- aggregate(factor(z[Id_remain]), list(group = Id_cluster[Id_remain]), summary)[-7, ]
    perc <- apply(info$x, 1, function(x) x/sum(x))
    for(j in 1:8){
      yposi <- perc[, j]*6 + pos[j, 2]
      xposi <- seq(pos[j, 1], by = 0.5, length.out = length(yposi))
      panel.addaxis(perc[, j], origin.x = pos[j, 1], origin.y = pos[j, 2], dx = 0.3, 
                    by = 0.5, box.width = 0.4, 
                    col = cols.avg, ntick = 2, border = "grey40", A = 5,
                    lwd = 0.5, text.cex = 1)
    }
  }
  z <- z[Id_remain]; 
  #20161123 modify data_in <- summary(factor(z)); prec <- data_in/sum(data_in);
  data_in <- sapply(1:max(z), function(x) length(which(z == x)))
  prec <- data_in/sum(data_in)
  ## TP summary infomation
  panel.addaxis(prec, col = cols.avg, border = "transparent", ...)
}

cols.fun <- function(n){
  ncol <- n - 1
  cols <- colorRampPalette(c("firebrick1","orange3", "darkgoldenrod2", "grey90", 
                             RColorBrewer::brewer.pal(9, "YlGnBu")[c(4, 6, 7)]))(ncol)#,colors()[504]
  cols <- cols[-5]
  cols[ncol] <- "green4"; #cols <- rev(cols)
  return(cols)
}

# panel.EOS_spatial <- function (x, y, z, subscripts, interpolate, ...,  sp.layout) 
# {
#    # print(str(listk(x, y, z, subscripts, ...)))#debug code
#    sppanel(list(sp.layout), panel.number(), first = TRUE)
#    panel.levelplot.raster(x, y, z, subscripts, ..., interpolate = interpolate)
#    # panel.contourplot(x, y, z, subscripts, ..., contour = TRUE, labels = F)
#    sppanel(list(sp.layout), panel.number(), first = FALSE)
   
#    i <- panel.number()
#    panel.text(76, 40.4, paste0("(",letters[i], ") ", zcol[i]), #english name: New_names[i])
#               fontfamily = "Times", cex = 1.3, font = 2, adj = 0)
#    panel.addbarchart(z, subscripts, cols, showCluster = showCluster, A = A, by = 0.9, ntick = ntick, ...)
# }
EOS_TP.spatial <- function(x, zcol = colnames(x), 
  brks, cols = cols.fun(length(brks)), layout = c(4, 4), 
  legendName = "DOY", saveFIG = TRUE, file, 
  showCluster = F, interpolate = T, 
  # panel = panel.EOS_spatial, 
  A = 8, by=0.9, ntick = 2, ...){
  gridClip@data <- llply(x, cut, breaks = brks, include.lowest = T) %>% do.call(cbind.data.frame,.)
  p <- spplot(gridClip, zcol = zcol, col.regions = cols, 
         # scales = list(draw=T), 
         # checkEmptyRC = T, #addNAemptyRowsCols, has no effects on the results, since it add NA values at tail.
         as.table = T, layout = layout,
         sp.layout = list(sp_cluster, sp_points, sp_text),
         xlim = xlim, ylim = ylim, strip = F, 
         # colorkey = list(labels = list(labels = levels(cut(1, brks)), at = seq_along(brks),
         colorkey = list(labels = list(labels = brks[-c(1, length(brks))], at = seq_along(brks[-(1:2)]) + 0.5,
                                       cex = 1.3, fontfamily = "Times", fontface = 2), 
                         axis.line = list(col = 'black'), 
                         height = 0.9), 
         panel = function (x, y, z, subscripts, ...,  sp.layout) 
        {
           # print(str(listk(x, y, z, subscripts, ...)))#debug code
           sppanel(list(sp.layout), panel.number(), first = TRUE)
           panel.levelplot.raster(x, y, z, subscripts, ..., interpolate = interpolate)
           # panel.contourplot(x, y, z, subscripts, ..., contour = TRUE, labels = F)
           sppanel(list(sp.layout), panel.number(), first = FALSE)
           
           i <- panel.number()
           panel.text(76, 40.4, paste0("(",letters[i], ") ", zcol[i]), #english name: New_names[i])
                      fontfamily = "Times", cex = 1.3, font = 2, adj = 0)
           panel.addbarchart(z, subscripts, cols, showCluster = showCluster, A = A, by = 0.9, ntick = ntick, ...)
        },
         par.settings = list(axis.line = list(col = "transparent")))

  width = 900; height = 520; ratio = 3.5
  if (saveFIG) CairoPNG(filename = file, width*ratio, height*ratio, dpi = 250, pointsize = 12)
  ifelse(saveFIG, print(p), return(p))
  # grid::grid.text(legendName, x=unit(0.945, "npc"), y=unit(0.936, "npc"), 
  #           gp=grid::gpar(fontfamily = "Times", fontsize=14, fontface = "bold"))
  if (saveFIG) dev.off()
}
# grid::grid.text(legendName, x=unit(0.945, "npc"), y=unit(0.936, "npc"), 
#                 gp=grid::gpar(fontfamily = "Times", fontsize=14, fontface = "bold"))