predict.coxphw <- function(object, type = c("shape", "slice.time", "slice.z", "slice.x", "lp", "risk"), x = NULL,
                           newx = NA, refx = NA, z = NULL, at = NULL,
                           exp = FALSE, se.fit = FALSE, pval = FALSE, digits = 4, verbose = FALSE, ...)
{

  ### Version 2015-12-20: new names for arguments.
  ### Version 2015-05-28: plotshape function was divied into predict.coxphw and plot.coxphw.predict
  ### version 2014-02-25:
  ### version 2010-04-29: NOT COMPATIBLE TO OLDER VERSIONS!!!
  ###
  ###
  ### shape: newx refx x
  ### slice.x: x z newx
  ### slice.z: x z at newx refx
  ### slice.time: x newx z
  ###
  ### object:   coxphf.beta object
  ### newx: range of x axis
  ### ref:    reference value (only for ref.type=="value" mode)
  ### treatment: variable which is in interaction with fp variable (use "")
  ### ref.type:  "value" for simple fp term (no interaction),
  ###            "interaction.time" for interaction for treatment with fp(time)
  ###            "interaction.treat" for interaction of treatment with fp(variable)
  ###            any value to specify the level of the treatment variable for which the fp of "variable" should be plotted
  ### variable:  name of fp variable (use "")

  if (missing(type)) stop ("type must be specified.")
  type <- match.arg(type)

  if (type == "shape") {
    if (any(is.na(newx), is.na(refx), is.null(x))) stop ("If type = shape => x, refx and newx must be specified.")
#    if (is.na(newx) || is.na(refx) || is.null(x)) stop ("If type = shape => x, refx and newx must be specified.")
    ref.type <- "value"
    z <- NULL
  } else
    if (type == "slice.time") {
      if (any(is.na(newx), is.null(z), is.null(x))) stop ("If type = slice.time => x, z and newx must be specified.")
#      if (is.na(newx) || is.null(z) || is.null(x)) stop ("If type = slice.time => x, z and newx must be specified.")
      ref.type <- "interaction.time"
    } else
      if (type == "slice.z") {
       if (any(is.na(newx), is.null(x), is.null(z))) stop ("If type = slice.z => x, z, at and newx must be specified.")
#       if (is.na(newx) || is.null(x) || is.null(z)) stop ("If type = slice.z => x, z, at and newx must be specified.")
        ref.type <- at
      } else
        if (type == "slice.x") {
          if (any(is.na(newx), is.null(x), is.null(z))) stop ("If type = slice.x => x, z and newx must be specified.")
#          if (is.na(newx) || is.null(x) || is.null(z)) stop ("If type = slice.x => x, z and newx must be specified.")
          ref.type <- "interaction.treat"
        }


  if(type == "lp")   {
    output <- object$linear.predictor
    if(verbose) print(output)
  }
  else if(type == "risk") {
    output <- exp(object$linear.predictor)
    if(verbose) print(output)
  } else {

    data<-object$dataline      # a dataline from original data set of analysis
    if (is.null(x)) stop("No x specified.\n")
    if (ref.type=="interaction.time" & is.null(z)) stop("No interaction variable for time specified.\n")
    if (is.na(refx)) refx <- min(newx)
    if (any(newx <= 0) & ref.type=="interaction.time") warning(paste("model.frame() might cancel lines with ", x, "<=0..."))
    n.x <- length(newx)
    data <- data[1,]
    data[,x]<-NA
    if(ref.type=="interaction.time" | ref.type=="interaction.treat") data[,z]<-1
    else if(ref.type != "value") data[,z]<-ref.type
    i.x <- which(is.na(data))
    plotData <- data[rep(1, n.x), , drop = FALSE]
    plotData[, i.x] <- newx
    rownames(plotData) <- 1:n.x
    objectP <- decomposeSurv(object$formula, plotData, sort = FALSE) #, object$offset, PTpreset = object$PTcoefs)
    kk <- ncol(objectP$mm1)
    objectP$mm1 <- objectP$mm1[, object$ind[1:kk], drop = FALSE]
#    objectP$covnames <- objectP$covnames[object$ind]
    objectP$timedata <- objectP$timedata[, object$ind[-(1:kk)], drop = FALSE]
#    plot.y <- (cbind(objectP$mm1, objectP$timedata)) %*% object$coefficients

    x1 <- cbind(objectP$mm1, objectP$timedata)

    if (ref.type != "interaction.time" & ref.type != "interaction.treat"){
        plotDataRef <- plotData
        plotDataRef[, i.x] <- refx
        objectPRef <- decomposeSurv(object$formula, plotDataRef, sort = FALSE) #, object$offset, PTpreset = object$PTcoefs)
        objectPRef$mm1 <- objectPRef$mm1[, object$ind[1:kk], drop = FALSE]
        objectPRef$timedata <- objectPRef$timedata[, object$ind[-(1:kk)], drop = FALSE]

        x0 <- cbind(objectPRef$mm1, objectPRef$timedata)
        diff <- x1-x0
        }
   if (ref.type == "interaction.treat"){
        plotDataRef <- plotData
        plotDataRef[, z] <- 0
        objectPRef <- decomposeSurv(object$formula, plotDataRef, sort = FALSE) #, object$offset, PTpreset = object$PTcoefs)
        objectPRef$mm1 <- objectPRef$mm1[, object$ind[1:kk], drop = FALSE]
        objectPRef$timedata <- objectPRef$timedata[, object$ind[-(1:kk)], drop = FALSE]

        x0 <- cbind(objectPRef$mm1, objectPRef$timedata)
        diff <- x1-x0
   }

   if (ref.type == "interaction.time") { diff <- x1 }

   plot.y <- diff %*% object$coefficients
   se2 <- sapply(1:nrow(diff), function(I) t(diff[I,]) %*% object$var %*% diff[I,])
   gammavt <- sqrt(qchisq(1-object$alpha, df=1) * se2)
   cilower <- plot.y - gammavt
   ciupper <- plot.y + gammavt

   if (pval) { p <- 1 - pchisq(plot.y^2 / se2, df = 1) } else { p <- NA }

   if (exp) {
     plot.y <- exp(plot.y)
     cilower <- exp(cilower)
     ciupper <- exp(ciupper)
   }

   estimates <- data.frame(newx, plot.y, cilower, ciupper)
   if (exp)  { dimnames(estimates)[[2]] <- c(x, "HR", paste(c("HR lower", "HR upper"), 1-object$alpha, sep=" ")) }
   else      { dimnames(estimates)[[2]] <- c(x, "coef", paste(c("coef lower", "coef upper"), 1-object$alpha, sep=" ")) }

   if (se.fit) { estimates <- cbind(estimates, se = sqrt(se2)) }
   if (pval) { estimates <- cbind(estimates, p) }
   if (verbose) print(round(estimates, digits = digits))

   output <- list(estimates = estimates, p = p, alpha = object$alpha, exp = exp, x = x)
   attr(output, "class") <- c("coxphw.predict")

  }
  return(invisible(output))
}
