# cut <- 

## ----------------------- main plot functions


stations.sp@data <- info[[2]]
id <- 4:6
spplot(stations.sp[, -seq(4, 14, 2)], main = "MAE", 
       col.regions = cols[id], edge.col = "black",
       as.table = TRUE, strip = F, layout = c(3, 3), 
       xlim = xlim, ylim = ylim, 
       cuts = brks, pch = pch, cex = cex, 
       key = list(points=list(pch=pch, cex=cex, fill=cols),
                text = list(levels(cut(1, brks))),
                font=2, cex=1.1, fontfamily='Times',columns=2),
       key.space = list(x = 0.72, y = 0.1, corner = c(0, 0)),
       sp.layout = list(sp_cluster), 
       panel = function(sp.layout, x, y, subscripts, col, ...){
         panel.pointsplot(sp.layout, x, y, subscripts, col, ...)
         i <- panel.number()
         panel.text(76, 40.4, paste0("(",letters[i], ") ", 
                                     vars[-seq(4, 14, 2)][i]), #english name: New_names[i])
                    fontfamily = "Times", cex = 1.3, font = 2, adj = 0)
       }, 
       par.settings = list(axis.line = list(col = "transparent")))

stations.sp@data <- info[[3]]
spplot(stations.sp[, -seq(4, 14, 2)], main = "RMSE", 
       col.regions = cols[4:6],  #col = "black", 
       edge.col = "black",
       # fill = cols.fun(brks), 
       as.table = TRUE, strip = F, layout = c(3, 3), 
       xlim = xlim, ylim = ylim, 
       # at = c(-20, -10, 0, 10, 20), 
       cuts = brks, col = "black",
       pch = pch, cex = cex, 
       key=list(points=list(pch=pch, cex=cex, col="black", fill=cols),
                text = list(levels(cut(1, brks))),
                font=2, cex=1.1,fontfamily='Times', columns=2),
       key.space = list(x = 0.72, y = 0.1, corner = c(0, 0)),
       sp.layout = list(sp_cluster), 
       panel = function(sp.layout, x, y, subscripts, col, ...){
         panel.pointsplot(sp.layout, x, y, subscripts, col, ...)
         i <- panel.number()
         panel.text(76, 40.4, paste0("(",letters[i], ") ", 
                                     vars[-seq(4, 14, 2)][i]), #english name: New_names[i])
                    fontfamily = "Times", cex = 1.3, font = 2, adj = 0)
       }, 
       par.settings = list(axis.line = list(col = "transparent")))