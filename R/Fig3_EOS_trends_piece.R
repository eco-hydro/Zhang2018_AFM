rm(list = ls())
source('R/mainfunc/mainfunctions.R', encoding = 'UTF-8')
source('R/mainfunc/Trends_func.R', encoding = 'UTF-8')
load("data/AVHRR_EOS_TP.rda")
# load("data/Fig1_agrdata2.rda")
load("data/00basemap.rda")
load("data/Fig1_infoValid.rda")

pheno <- phenoEOS[c(1, 3, 7, 9)]

mean <- mean(pheno)
pheno_mean <- array(NA, dim=c(22204, 32, 4))
for (i in 1:4){
  pheno_mean[,, i] <- as.matrix(pheno[[i]])
}

pheno$Mean <- apply(pheno_mean, c(1,2), mean, na.rm = T) %>% as.data.frame()
pheno %<>% llply(function(x) set_colnames(as.matrix(x), NULL))
x <- llply(pheno, apply, 1, mean, na.rm = T, .id = NULL, .progress = "text") %>% 
  do.call(cbind.data.frame, .)

gridClip@data <- x
spplot(gridClip)


EOS_TP.spatial <- function(x, zcol = vars, brks, legendName = "DOY", saveFIG = TRUE, file, 
                           A = 8, by=0.9, ntick = 2, ...){
  ncol <- length(brks) - 1
  cols <- colorRampPalette(c("firebrick1","orange3", "darkgoldenrod2", "grey90", 
                             brewer.pal(9, "YlGnBu")[c(4, 6, 7)]))(ncol)#,colors()[504]
  cols <- cols[-5]
  cols[ncol] <- "green4"; #cols <- rev(cols)
  
  gridClip@data <- llply(x, cut, breaks = brks, include.lowest = T) %>% do.call(cbind.data.frame,.)
  p <- spplot(gridClip, zcol = zcol, col.regions = cols, 
              # scales = list(draw=T), 
              as.table = T,
              # sp.layout = list(sp_cluster, sp_points, sp_text),
              sp.layout = list(sp_cluster),
              xlim = xlim, ylim = ylim, strip = F, 
              # colorkey = list(labels = list(labels = levels(cut(1, brks)), at = seq_along(brks),
              colorkey = list(labels = list(labels = brks[-c(1, length(brks))], at = seq_along(brks[-(1:2)]) + 0.5,
                                            cex = 1.3, fontfamily = "Times", fontface = 2), 
                              axis.line = list(col = 'black'), 
                              height = 0.9), 
              panel = function (x, y, z, subscripts, ...,  sp.layout) 
              {
                sppanel(list(sp.layout), panel.number(), first = TRUE)
                panel.levelplot.raster(x, y, z, subscripts, ..., interpolate = T)
                # panel.contourplot(x, y, z, subscripts, ..., contour = TRUE, labels = F)
                sppanel(list(sp.layout), panel.number(), first = FALSE)
                
                i <- panel.number()
                panel.text(76, 40.4, paste0("(",letters[i], ") ", zcol[i]), #english name: New_names[i])
                           fontfamily = "Times", cex = 1.3, font = 2, adj = 0)
                # panel.text(76.6, 40.8 +heights[i], metricsEN[i], fontfamily = "Times", cex = 1.2, font = 2, adj = 0)
                panel.addbarchart(z, subscripts, cols, showCluster = F, A = A, by = 0.9, ntick = ntick, ...)
              },
              par.settings = list(axis.line = list(col = "transparent")))
  
  width = 900; height = 520; ratio = 3.5
  if (saveFIG) CairoPNG(filename = file, width*ratio, height*ratio, dpi = 250, pointsize = 12)
  print(p)
  grid::grid.text(legendName, x=unit(0.945, "npc"), y=unit(0.936, "npc"), 
            gp=grid::gpar(fontfamily = "Times", fontsize=14, fontface = "bold"))
  if (saveFIG) dev.off()
}

# brks <- c(10, seq(240, 300, 10), 366)
brks <- c(0.01, 0.05, 0.1, 0.5, 10) %>% {c(-rev(.), 0, .)}
EOS_TP.spatial(x, zcol = colnames(x), brks, saveFIG = F)



EOS_TP.spatial(x, zcol = 5, brks, saveFIG = F)

## 进行并行计算
options(warn = -1)
tmp <-clusterEvalQ(cl, {
  # library(magrittr)
  library(segmented)
  # library(plyr)
  # library(optimx)
  # library(zoo)
  #close all warnings, due to optimx check initial start point, can ignore it
  options(warn =-1)
})
clusterExport(cl, "piecewise", envir = globalenv())

# for (i in seq_along(NDVI_list)) NDVI_list[[i]] <- list(xx=NDVI_list[[i]], i=i)

psiInfo <- list()
for (i in 1:5){
  cat(sprintf("[%d]th --------\n", i))
  x_list <- alply(pheno[[i]], 1, function(x) x) 
  psiInfo[[i]] <- parLapplyLB(cl,x_list, piecewise) %>% do.call(rbind, .)
}
save(psiInfo, file = "data/Fig3_FivePresent_psiInfo.rda")

load("data/Fig3_FivePresent_psiInfo.rda")
# psiInfo <- list()
info <- list()
for (i in 2:5){
  cat(sprintf("[%d]th --------\n", i))
  x_list <- alply(pheno[[i]], 1, function(x) x) 
  # psiInfo[[i]] <- parLapplyLB(cl,x_list, piecewise) %>% do.call(rbind, .)
  info[[i]] <- MankKendall.Trend_apply(pheno[[i]], psiInfo[[i]]$psi, psiInfo[[i]]$conf, trim=2)
}
info %<>% set_names(names(pheno))

save(info, file = "data/Fig3_FivePresent_TrendInfo.rda")
load("data/Fig3_FivePresent_TrendInfo.rda")

# final plot EOS TRENDS ---------------------------------------------------
rm(list = ls())
source('R/mainfunc/mainfunctions.R', encoding = 'UTF-8')
load("data/00basemap.rda")
load("data/Fig3_FivePresent_TrendInfo.rda")


# brks <- c(2, 7, 15, 50) %>% c(-rev(.), .)
# p <- EOS_TP.spatial(trend, brks = brks, saveFIG = F, layout =c(3, 2), 
#                     showCluster = T, interpolate = T)

## Mann Kendall Z values
brks <- c(1.65, 1.96, 5) %>% c(-rev(.), 0, .)
trend_all <- llply(info, `[[`, "Z") %>% do.call(cbind.data.frame,.)
trend_before <- llply(info, `[[`, "before.Z") %>% do.call(cbind.data.frame,.)
trend_after <- llply(info, `[[`, "after.Z") %>% do.call(cbind.data.frame,.)

EOS_TP.spatial(trend_all, brks = brks, saveFIG = F, layout =c(3, 2), 
               showCluster = F, interpolate = T, cols = cols)

EOS_TP.spatial(trend_before, brks = brks, saveFIG = F, layout =c(3, 2), 
               showCluster = F, interpolate = T, cols = cols)

EOS_TP.spatial(trend_after, brks = brks, saveFIG = F, layout =c(3, 2), 
               showCluster = F, interpolate = T, cols = cols)