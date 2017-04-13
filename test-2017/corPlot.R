function (r, numbers = FALSE, colors = TRUE, n = 51, main = NULL, 
          zlim = c(-1, 1), show.legend = TRUE, labels = NULL, n.legend = 10, 
          keep.par = TRUE, select = NULL, pval = NULL, cuts = c(0.001, 
                                                                0.01), cex, MAR, upper = TRUE, diag = TRUE, ...) 
{
  if (keep.par) 
    op <- par(no.readonly = TRUE)
  if (missing(MAR)) 
    MAR <- 5
  if (is.null(main)) {
    main <- "Correlation plot"
  }
  if (!is.matrix(r) & (!is.data.frame(r))) {
    if ((length(class(r)) > 1) & (class(r)[1] == "psych")) {
      switch(class(r)[2], omega = {
        r <- r$schmid$sl
        nff <- ncol(r)
        r <- r[, 1:(nff - 3)]
      }, cor.ci = {
        pval <- 2 * (1 - r$ptci)
        r <- r$rho
      }, fa = {
        r <- r$loadings
      }, pc = {
        r <- r$loadings
      }, principal = {
        r <- r$loadings
      })
    }
  }
  r <- as.matrix(r)
  if (min(dim(r)) < 2) {
    stop("You need at least two dimensions to make a meaningful plot")
  }
  minx <- min(r, na.rm = TRUE)
  maxx <- max(r, na.rm = TRUE)
  if ((minx < -1) | (maxx > 1)) 
    r <- cor(r, use = "pairwise")
  if (is.null(n)) {
    n <- dim(r)[2]
  }
  nf <- dim(r)[2]
  nvar <- dim(r)[1]
  if (!upper) 
    r[col(r) > row(r)] <- NA
  if (!diag) 
    r[col(r) == row(r)] <- NA
  if (nf == nvar) 
    r <- t(r)
  if (missing(pval) | is.null(pval)) {
    pval <- matrix(rep(1, nvar * nf), nvar)
  }
  else {
    if (length(pval) != nvar * nf) {
      pr = matrix(0, nvar, nf)
      pr[row(pr) > col(pr)] <- pval
      pr <- pr + t(pr)
      diag(pr) <- 0
      pval <- pr
    }
    pval <- con2cat(pval, cuts = cuts)
    pval <- (length(cuts) + 1 - pval)/length(cuts)
    pval <- t(pval)
  }
  if (is.null(labels)) {
    if (is.null(rownames(r))) 
      rownames(r) <- paste("V", 1:nvar)
    if (is.null(colnames(r))) 
      colnames(r) <- paste("V", 1:nf)
  }
  else {
    rownames(r) <- colnames(r) <- labels
  }
  max.len <- max(nchar(rownames(r)))/6
  if (is.null(zlim)) {
    zlim <- range(r)
  }
  if (colors) {
    gr <- colorRampPalette(c("red", "white", "blue"))
    colramp <- gr(n)
  }
  else {
    colramp <- grey((n:0)/n)
  }
  if (nvar != nf) {
    r <- t(r)
  }
  if (!is.null(select)) {
    r <- r[select, select]
    pval <- pval[select, select]
    nvar <- length(select)
  }
  ord1 <- seq(nvar, 1, -1)
  if (nf == nvar) {
    r <- r[, ord1]
    pval <- pval[, ord1]
  }
  else {
    r <- r[, ord1]
    pval <- t(pval[ord1, ])
  }
  par(mar = c(MAR + max.len, MAR + max.len, 4, 0.5))
  if (show.legend) {
    layout(matrix(c(1, 2), nrow = 1), widths = c(0.9, 0.1), 
           heights = c(1, 1))
  }
  image(r, col = colramp, axes = FALSE, main = main, zlim = zlim)
  box()
  if (nf < nvar) {
    at1 <- (0:(nf - 1))/(nf - 1)
  }
  else {
    at1 <- (0:(nvar - 1))/(nvar - 1)
  }
  at2 <- (0:(nvar - 1))/(nvar - 1)
  if (max.len > 0.5) {
    axis(2, at = at2, labels = colnames(r), las = 1, ...)
    axis(1, at = at1, labels = rownames(r), las = 2, ...)
  }
  else {
    axis(2, at = at2, labels = colnames(r), ...)
    axis(1, at = at1, labels = rownames(r), las = 1, ...)
  }
  if (numbers) {
    rx <- rep(at1, ncol(r))
    ry <- rep(at2, each = nrow(r))
    rv <- round(r, 2)
    if (missing(cex)) 
      cex = 10/max(nrow(r), ncol(r))
    text(rx, ry, rv, cex = pval * cex, ...)
  }
  if (show.legend) {
    leg <- matrix(seq(from = zlim[1], to = zlim[2], by = (zlim[2] - 
                                                            zlim[1])/n), nrow = 1)
    par(mar = c(MAR, 0, 4, 3))
    image(leg, col = colramp, axes = FALSE, zlim = zlim)
    at2 <- seq(0, 1, 1/n.legend)
    labels = seq(zlim[1], zlim[2], (zlim[2] - zlim[1])/(length(at2) - 
                                                          1))
    axis(4, at = at2, labels = labels, las = 2, ...)
  }
  if (keep.par) 
    par(op)
}