print.coxphw.predict <- function(x, ...)
  
{
  if(!inherits(x, "coxphw.predict")) stop("Must be an 'coxphw.predict' object")
  
  print(x$estimates)
  
  invisible(x)
}