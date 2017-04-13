spplot.points = function(obj, zcol = names(obj), ..., names.attr, 
                         scales = list(draw = FALSE), xlab = NULL, ylab = NULL, 
                         aspect = mapasp(obj,xlim,ylim), panel = panel.pointsplot,
                         sp.layout = NULL, identify = FALSE, formula,
                         xlim = bbexpand(bbox(obj)[1,], 0.04), 
                         ylim = bbexpand(bbox(obj)[2,], 0.04),
                         edge.col = "transparent", colorkey = FALSE,
                         col.regions = get_col_regions()) 
{
  
  if (is.null(zcol)) stop("no names method for object")
  dots = list(...)
  sdf = obj
  if (!is.character(zcol)) 
    zcol = names(sdf)[zcol]
  # create formula:
  if (missing(formula)) {
    if (length(zcol) > 1) {
      formula = as.formula(paste(paste(dimnames(coordinates(sdf))[[2]][2:1], 
                                       collapse = "~"), "|name"))
      sdf = spmap.to.lev(sdf, zcol = zcol, names.attr = names.attr)
    } else {
      if (!is.character(zcol)) 
        zcol = names(sdf)[zcol]
      ccn = dimnames(coordinates(sdf))[[2]]
      formula = as.formula(paste(ccn[2], "~", ccn[1]))
    }
  }
  scales = longlat.scales(obj, scales, xlim, ylim)
  args.xyplot = append(list(formula, data = as(sdf, "data.frame"), 
                            panel = panel, aspect = aspect, scales = scales, 
                            xlab = xlab, ylab = ylab, sp.layout = sp.layout,
                            xlim = xlim, ylim = ylim, edge.col = edge.col,
                            col.regions = col.regions), dots)
  z = create.z(as(obj, "data.frame"), zcol)
  args.xyplot = fill.call.groups(args.xyplot, z = z, edge.col = edge.col, 
                                 colorkey = colorkey, ...)
  # debug:
  #print(args.xyplot)
  plt = do.call(xyplot, args.xyplot)
  if (!(is.logical(identify) && identify==FALSE) && interactive()) {
    print(plt)
    if (!(is.numeric(identify) && length(identify) == 2))
      identify = c(1,1)
    trellis.focus("panel", identify[1], identify[2])
    labels = row.names(as(sdf, "data.frame"))
    cat("left-mouse to identify points; right-mouse to end\n")
    cc = coordinates(obj)
    ret = panel.identify(cc[,1], cc[,2], labels)
    trellis.unfocus()
    return(ret)
  } else
    plt
}

create.z = function(df, zcol) {
  if (is.logical(df[[zcol[1]]])) {
    z = stack(df[zcol])[[1]]
    z = as.factor(z)
  } else if (is.numeric(df[[zcol[1]]]))
    z = stack(df[zcol])[[1]]
  else if (is.factor(df[[zcol[1]]])) {
    lev = levels(df[[zcol[1]]])
    z = factor(as.vector(sapply(df[zcol], as.character)), levels = lev)
  } else
    stop("no support for variable of this type")
  z
}

fill.call.groups <-
  function (lst, z, ..., cuts = ifelse(identical(FALSE, colorkey), 5, 100), 
            #col.regions = trellis.par.get("regions")$col, 
            legendEntries = "", pch, cex = 1, do.fill = TRUE, do.log = FALSE, 
            key.space = ifelse(identical(FALSE, colorkey), "bottom", "right"), 
            cex.key, edge.col, colorkey) 
  {
    dots = list(...)
    col.regions = lst$col.regions
    if (is.numeric(z)) {
      if (length(cuts) > 1) 
        ncuts = length(cuts) - 1
      else ncuts = cuts
      if (ncuts != length(col.regions)) {
        cols = round(1 + (length(col.regions) - 1) * (0:(ncuts - 
                                                           1))/(ncuts - 1))
        fill = col.regions[cols]
      } else 
        fill = col.regions
      valid = !is.na(z)
      if (length(cuts) == 1) {
        if (do.log) {
          lz = log(z)
          cuts = c(min(z[valid]), exp(seq(min(lz[valid]), 
                                          max(lz[valid]), length = cuts + 1))[2:(cuts)], 
                   max(z[valid]))
        }
        else cuts = seq(min(z[valid]), max(z[valid]), length = cuts + 
                          1)
      }
      groups = cut(as.matrix(z), cuts, dig.lab = 4, include.lowest = TRUE)
    } else if (is.factor(z)) {
      if (length(col.regions) == 1) 
        col.regions = rep(col.regions, nlevels(z))
      if (length(col.regions) < nlevels(z)) 
        stop("number of colors smaller than number of factor levels")
      if (length(col.regions) > nlevels(z)) {
        ncuts = nlevels(z)
        cols = round(1 + (length(col.regions) - 1) * (0:(ncuts - 
                                                           1))/(ncuts - 1))
        col.regions = col.regions[cols]
      }
      if (!missing(cuts)) 
        stop("ncuts cannot be set for factor variable")
      groups = z
      fill = col.regions
    } else stop("dependent of not-supported class")
    n = nlevels(groups)
    
    # deal with col:
    lst$groups = fill[groups]
    #print(lst$col)
    
    # deal with pch:
    if (edge.col != "transparent") { # WITH border: use fill
      if (missing(pch)) 
        pch = rep(ifelse(do.fill, 21, 1), n)
      lst$col = rep(edge.col, length.out = length(groups))
    } else { # no border: use col instead of fill
      if (missing(pch)) 
        pch = rep(ifelse(do.fill, 16, 1), n)
      lst$col = lst$groups
    }
    
    if (length(pch) == 1)
      pch = rep(pch, n)
    lst$pch = pch[groups]
    
    # deal with cex:
    if (missing(cex))
      cex = rep(1, n)
    if (length(cex) == 1)
      cex = rep(cex, n)
    if (length(cex) == n) {
      cex.key = cex
      lst$cex = cex[groups]
    } else if (missing(cex.key))
      cex.key = mean(cex, na.rm = TRUE)
    
    # do key:
    if (is.list(colorkey))
      lst$legend = colorkey
    else if (isTRUE(colorkey)) {
      lst$legend = list(
        right = list(
          fun = draw.colorkey,
          args = list(
            key = list(
              col = col.regions, 
              at = cuts
            ), 
            draw = FALSE
          )
        )
      )
      if (is.character(key.space)) 
        names(lst$legend) = key.space
    } else {
      if (!identical(dots$auto.key, FALSE)) { # xxx
        if (missing(legendEntries)) 
          legendEntries = levels(groups)
        if (!is.null(dots$key)) 
          lst$key = dots$key
        else { 
          if(is.list(dots$auto.key))
            lst$key = dots$auto.key
          else
            lst$key = list()
          if (edge.col != "transparent") {
            lst$key = append(lst$key,
                             list(points = list(
                               pch = rep(pch, length.out = n), 
                               col = rep(edge.col, length.out = n), 
                               fill = fill, 
                               cex = rep(cex.key, length.out = n)
                             ), 
                             text = list(legendEntries)
                             ))
          } else {
            lst$key = append(lst$key,
                             list(points = list(
                               pch = rep(pch, length.out = n), 
                               col = rep(fill, length.out = n), 
                               cex = rep(cex.key, length.out = n)
                             ), 
                             text = list(legendEntries)
                             ))
          }
        }
        if (is.character(key.space)) 
          lst$key$space = key.space
        else if (is.list(key.space)) 
          lst$key = append(lst$key, key.space)
        else warning("key.space argument ignored (not list or character)")
        # print(lst$key)
      }
      if (!is.null(dots$auto.key)) 
        lst$auto.key <- dots$auto.key
    }
    return(lst)
  }
# environment(spplot.points) <- environment(sp:::spplot.points)
# setMethod("spplot", signature("SpatialPointsDataFrame"), spplot.points)