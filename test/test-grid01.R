library(grid)
library(gridExtra)

source('R/mainfunc/main_PLOT.R', encoding = 'UTF-8')
load("../data/00basemap.rda")
bbox <- gridClip@bbox
# p <- VSI_show(x, brks_z)

## fig. 8 auto correlation and pearson correlation
#  kcor, acor

## long-term solution is to change spplot
x <- data.frame(acor, kcor)
brks <- list(c(Rc_acf, 1), 
             c(Rc_cor, 1))
VSI_show.subList(x, brks)


xInfo <- signPortion(x$acor, Rc_acf) %T>% print
info <- aggregate(x$acor, list(cluster = Id_cluster), signPortion, brks = c(Rc_acf, 1)) %$% 
  cbind.data.frame(cluster, x) %T>% print 

signPortion <- function(x, brks){
  Id_cluster[Id_cluster == "热带雨林"] <- NA
  brks <- c(brks, Inf) 
  breaks <- round(brks, 3) %>% {c(-rev(.), 0, .)}
  xx <- cut(x, breaks)
  
  info <- aggregate(xx, list(cluster =  Id_cluster), summary) %$% 
    cbind.data.frame(cluster = c(as.character(cluster), "TP"), 
                     apply(rbind(x, colSums(x)), 1, function(x) round(x/sum(x), 4)*100) %>% t)
  info
  # %>% count()
  # info$per <- round(info$freq/sum(info$freq), 4)*100
  # setNames(info$per, info$x)
  # setNames(info$per, info$freq)
}

ports <- llply(trend, signPortion, brks = Z)
writelist_ToXlsx(ports, fname = "Pheno_trend.xlsx")
signPortion(trend$EOS, Z)

# p1 <- VSI_show(x[, 1, drop = F], c(Rc_acf, 1), zcol = 1)
# p2 <- VSI_show(x[, 2, drop = F], c(Rc_cor, 1), zcol = 1)

# gridExtra::marrangeGrob(p1, p2, nrow = 1, ncol = 2)
# grid.arrange(list(p1, p2), nrow = 1)

## tests script
grid.rect(gp = gpar(fill = "transparent", col = "red"))

current.viewport()
current.vpPath()
current.vpTree()
downViewport("plot_09.panel.1.1.off.vp")
upViewport(0)

showViewport(newpage = T, leaves = T)

grid.ls(viewports = T, fullNames = T)
current.viewport()
current.vpTree()

# according to grep search assigned name viewport
vps <- grid.ls(grob = F, viewports = T, print = F)
vps_panel <- str_subset(unique(vp$name), ".+panel\\.1\\.1\\.vp") %>% sort(decreasing = T)
