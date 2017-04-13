kong <- function(){
  spplot(gridClip, sp.layout = sp_points,  
         panel =  function (x, y, z, subscripts, ..., sp.layout) 
         {
           sppanel(list(sp.layout), panel.number(), first = TRUE)
           panel.levelplot(x, y, z, subscripts, ...)
           
           loc <- coordinates(stations)
           panel.text(x = loc[, 1], y = loc[, 2] + 0.5, 
                      label = stations$stationId, col = "red")
           sppanel(list(sp.layout), panel.number(), first = FALSE)
         })
}
