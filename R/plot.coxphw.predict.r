plot.coxphw.predict <- function(x, addci = TRUE,  xlab = NULL, ylab = NULL, ...)

{
  if(!inherits(x, "coxphw.predict")) stop("Must be an 'coxphw.predict' object")

   refline <- 0
   if(x$exp) {
     if(is.null(ylab)) { ylab <- "relative hazard" }
     refline <- 1
   }

  if(is.null(xlab)) { xlab <- x$x}
  if(is.null(ylab)) { ylab <- "log relative hazard" }
  if(addci) { ylimit <- range(x$estimates[,2:4]) } else { ylimit <- range(x$estimates[,2]) }
  plot(x$estimates[,1], x$estimates[,2], lty = 1, xlab = xlab, ylab = ylab, ylim = ylimit, type = "l", las = 1, ...)
  if(addci) {
    lines(x = x$estimates[,1], y = x$estimates[,3], lty = 2)
    lines(x = x$estimates[,1], y = x$estimates[,4], lty = 2)
  }
  abline(refline, 0, col = "gray", lty = 2)
}
