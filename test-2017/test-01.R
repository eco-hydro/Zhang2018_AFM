alpha <- c(0.1, 0.05)
pc_cor <- get_critical(alpha, 32 - 2)
pc_acf <- get_critical(alpha, 31 - 2)
# brks_acf <- qnorm(1 - alpha/2)/sqrt(32)


brks_acf <- c(pc_cor, 1)
brks_z <- c(qnorm(1 - alpha/2), 10)
brks_slp <- c(0.1, 0.5, 1, 10)

VSI_show <- function(x, brks){
	brks %<>% round(3) %>% {c(-rev(.), 0, .)}
	zcol <- colnames(x)
	nmonths <- 7
	nvar <- ncol(x)/nmonths
	A = 8; by=0.9; ntick = 2
	# brks <- c(Rc, 1) %>% {c(-rev(.), 0, .)}

	# ncol <- length(brks) + 1
	# cols <- colorRampPalette(c("firebrick1","orange3", "darkgoldenrod2", "grey90", 
	#                            brewer.pal(9, "YlGnBu")[c(4, 6, 7)]))(ncol)#,colors()[504]
	# # cols <- cols[-ceiling(ncol/2)]
	# cols[ncol+1] <- "green4"; #cols <- rev(cols)
	# cols <- cols[-c(3, 5, 7)]
	# cols[3:4] <- "grey90"

	ncol <- length(brks) - 1
	# "grey90", YlGnBu:4, 6, 7
	cols <- colorRampPalette(c("firebrick1","orange3", "darkgoldenrod2", "grey90",
	                         brewer.pal(9, "YlGnBu")[c(4, 6, 7)], "green4"))(ncol)#,colors()[504]
	# cols <- cols[-ceiling(ncol/2)]
	# cols[ncol+1] <- "green4"; #cols <- rev(cols)
	# show_col(cols)

	gridClip@data <- llply(x, cut, breaks = brks, include.lowest = T) %>% do.call(cbind.data.frame,.)
	p <- spplot(gridClip, zcol = zcol, col.regions = cols, 
	          # scales = list(draw=T), 
	          as.table = T, 
	          # layout = c(nvar, nmonths),
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
	            panel.levelplot.raster(x, y, z, subscripts, ..., interpolate = F)
	            # panel.contourplot(x, y, z, subscripts, ..., contour = TRUE, labels = F)
	            sppanel(list(sp.layout), panel.number(), first = FALSE)
	            
	            i <- panel.number()
	            panel.text(76, 40.4, paste0(zcol[i]), #english name: New_names[i])
	                       fontfamily = "Times", cex = 1.3, font = 2, adj = 0)#"(",letters[i], ") ", 
	            # panel.text(76.6, 40.8 +heights[i], metricsEN[i], fontfamily = "Times", cex = 1.2, font = 2, adj = 0)
	            panel.addbarchart(z, subscripts, cols, showCluster = T, A = A, by = 0.9, ntick = ntick, origin.x = 77, ...)
	          },
	          par.settings = list(axis.line = list(col = "transparent")))

	# CairoPNG(filename, 14, 11, dpi = 300, pointsize = 12, units = "in")
	# # spplot(gridClip, as.table = T, layout = c(5, 12))
	print(p)
	# dev.off()
	# width = 900; height = 520; ratio = 3.5
	# if (saveFIG) CairoPNG(filename = file, width*ratio, height*ratio, dpi = 250, pointsize = 12)
	# print(p)
	# grid::grid.text(legendName, x=unit(0.945, "npc"), y=unit(0.936, "npc"), 
	#                 gp=grid::gpar(fontfamily = "Times", fontsize=14, fontface = "bold"))
	# if (saveFIG) dev.off()
}