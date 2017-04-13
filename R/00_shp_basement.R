library(magrittr)
library(maptools)
library(data.table)
source("E:/GitHub/statistics_Kong.R", encoding = "utf-8")

Id_del <- c(16035, 18443, 20925, 21407, 21618, 21732, 21827)#need to cautious about this

olddir <- getwd()
setwd("E:/GitHub/Vegetation/data/")
## 青藏高原范围:25-40, 73-105, exactly as matlab data input
get_grid <- function(range, cellsize){
  lat_range <- range[1:2]; long_range <- range[3:4]
  # lat_range <- c(25, 40); long_range <- c(73, 105)
  dims <- c(diff(long_range), diff(lat_range))/cellsize
  grid <- GridTopology(cellcentre.offset = c(long_range[1], lat_range[1]) + cellsize/2, 
                       cellsize = c(1, 1) * cellsize, cells.dim = dims)
  grid <- SpatialPixelsDataFrame(grid, data = data.frame(id = 1:prod(dims)), proj4string = prj)
  return(grid)
  # grid@data <- data.frame(as.numeric(df))
  # spplot(grid)
}

cellsize <- 1/12
range <- c(25, 40, 73, 105)
grid <- get_grid(range, cellsize)

Id_avhrr <- fread("Id_avhrr.txt")$V1#[-Id_del]
gridClip <- grid[Id_avhrr, ][-Id_del, ]
gridClip$id <- 1:length(gridClip)

## --------------------- 青藏高原面图层 ---------------------------
TP_poly <- readShapePoly("../data/shp/TP_poly.shp", proj4string = prj)
sp_poly <- list("sp.polygons", TP_poly, fill = "transparent", col = "black", first = F)
xlim = bbox(TP_poly)[1, ] + c(0, 0.3); ylim =  bbox(TP_poly)[2, ] + c(0, 1.8)
# shp <- fortify(TP_poly)
## 植被分区
cluster <- readShapePoly("shp/TP_vegZoneSolve84.shp", proj4string = CRS("+proj=longlat +ellps=WGS84"))
Id_cluster <- over(gridClip, cluster)$Name
# Id_cluster.num <- as.numeric(Id_cluster)
x <- factor(as.character(Id_cluster), levels = levels(Id_cluster)[-7]) %>% as.numeric(); x[is.na(x)] <- -1
fwrite(data.frame(x), file="cluster_avhrr.txt", col.names = F)

Id_cluster.list <- lapply(levels(Id_cluster), function(x) which(Id_cluster == x)) %>% set_names(levels(Id_cluster))
sp_cluster <- list("sp.polygons", cluster, fill = "transparent", col = "black", first = F)

## --------------------- 农业物候站点 ---------------------------
#  56312 station not in research regions
stations <- read.table("stations.txt", header = T, sep = "\t", stringsAsFactors = F)[-25, ]
stations$VId <- sprintf("v%02d_%d", 1:nrow(stations), stations$stationId)

stations.sp <- stations
coordinates(stations.sp) <- ~long + lat
proj4string(stations.sp) <- prj
# Id_clip <- over(stations, gridClip)$id

loc <- coordinates(stations.sp)
loc[, 2] %<>% + 0.3
sp_text <- list("sp.text", loc, col = "black", cex = 0.6,
                txt = paste0("", 1:nrow(stations)))
sp_points <- list("sp.points", stations.sp, col = "red", cex = 0.5)
# spplot(gridClip, sp.layout = list(sp_text, sp_points, sp_poly))

Id_clip <- over(stations.sp, gridClip)$id
## 把这些数据均放在一个env里面
# .options <- listk(prj, gridClip, stations, stations.sp, sp_text, 
#                   sp_points, sp_cluster)
save(prj, gridClip, stations, stations.sp, Id_clip, sp_text, sp_points, sp_cluster, xlim, ylim, 
     Id_cluster, Id_cluster.list,
     file = "00basemap.rda")

setwd("../TP_phenology/")