# -------------- global variables
library(grid)
load("x.rda")
load("../data/00basemap.rda")
bbox <- bbox(gridClip)

# functions

# grid.grabExpr <- function (expr, warn = 2, wrap = FALSE, ...) 
# {
#   pdf(file = NULL)
#   on.exit(dev.off())
#   eval(expr)
#   res <- grid:::grabDL(warn, wrap, ...)
#   return(res)
# }

panel1 <- function (x, y, z, subscripts, ...,  sp.layout) {
  sppanel(list(sp.layout), panel.number(), first = TRUE)
  panel.levelplot.raster(x, y, z, subscripts, ...)
  # panel.contourplot(x, y, z, subscripts, ..., contour = TRUE, labels = F)
  sppanel(list(sp.layout), panel.number(), first = FALSE)
}

panel2 <- function (x, y, z, subscripts, ...,  sp.layout) {
  dev.new()
  panel1(x, y, z, subscripts, ...,  sp.layout = sp.layout)
  # a2 <- grid.grab()
  dev.off()
  
  p <- panel1(x, y, z, subscripts, ...,  sp.layout = sp.layout)
  
  i <- panel.number()
  # vp <- 
  zc <- z; zc[abs(zc) < brks[1]] <- NA
  
  print(current.viewport())
  print(current.vpPath())
  # current.vpTree()
  
  bbox <- bbox(gridClip)


  str(list(...))
  debug(grid.grabExpr)
  # undebug(grid.grabExpr)
  
  a <- grid.grab()
  
  
  
  vp <- vpStack(
    viewport(
      name = "ozlay",
      layout = grid.layout(1, 1, widths = diff(bbox[1,]), heights = diff(bbox[2,]), respect = T)
    ),
    viewport(
      name = "ozvp", 
      layout.pos.row = 1, layout.pos.col = 1,
      xscale = bbox[1,], yscale = bbox[2,], 
      clip = TRUE
    )
  )
  
  pd <-  grid.grabExpr()
  all.equal(a, pd)
  
  pd <-  grid.grabExpr(panel1(x, y, z, subscripts, ...,  sp.layout = sp.layout))
  
  pd <-  grid.grabExpr(grid.rect(gp = gpar(col = "red",fill = "black")))
  
  downViewport("plot_01.panel.1.1.vp")
  gTree(
    name = "add",
    vp = viewport(x = -0.01, y = 1.1, just = c("left", "top"), width = 0.55, height = 0.55),
    childrenvp = vp,
    children = gList(a),
    cl = "ozGrob"
  ) %>% grid.draw()
  # grid.draw(pp)
  
  upViewport(0)
  # pd <- grid.grabExpr(print({
  #   sppanel(list(sp.layout), panel.number(), first = TRUE)
  #   panel.levelplot.raster(x, y, zc, subscripts, ...)
  #   # panel.contourplot(x, y, z, subscripts, ..., contour = TRUE, labels = F)
  #   sppanel(list(sp.layout), panel.number(), first = FALSE)
  # }))
}

spplot(obj = gridClip,
         zcol = zcol,
         col.regions = cols,
         # scales = list(draw=T),
         as.table = T,
         # layout = c(nvar, nmonths),
         # sp.layout = list(sp_cluster, sp_points, sp_text),
         sp.layout = list(sp_cluster),
         xlim = xlim, ylim = ylim, 
         strip = strip.custom(factor.levels = paste0("(", letters[1:ncol(x)], ") ", zcol)), 
         # colorkey = list(labels = list(labels = levels(cut(1, brks)), at = seq_along(brks),
         colorkey = colorkey, 
         interpolate = T,
         par.strip.text = list(fontfamily = "Times", cex = 1.1, font = 2), 
         panel = panel2)
# pushViewport(viewport(x = -0.01, y = 1.1, just = c("left", "top"),width = 0.55, height = 0.55))
# pushViewport(vp)
# panel1(x, y, zc, subscripts, ...,  sp.layout = sp.layout)
