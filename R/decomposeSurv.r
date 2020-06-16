"decomposeSurv" <- function(
  formula,
  data,
  sort=FALSE
)
### Decomposes survival formula with time-factor interactions and offset in formula.
### MP, 2015/09-12.
{
  if (any(is.na(data)))
    stop("Data frame <data> must not contain any missing values.")

  effectsOriginal = colnames(model.matrix(formula, data = data))

  ######## Helper Functions Start #######
  ## Convenience formula function. 
  addTermsToFormula = function(terms) {
    as.formula(paste(as.character(formula)[2], "~",
                     as.character(formula)[3], "+",
                     paste(terms, sep="+")))
  }
  
  ## Convenience Split functions by special characters or simply ':'.
  f <- function(str) {
    for(chars in c("(", ")", ":", " ", ",", "*", "^"))
      str <- unlist(strsplit(str, split=chars, fixed=TRUE))
    str
  }
  f2 = function(str) {
    unlist(strsplit(str, split=":", fixed=TRUE))
  }
  ######## Helper Functions End #########
  
  ## Add plain effects of time-effect interactions to construct whole model matrix.
  repeat {
    ## For fractional polynomials add here fp() beside offset.
    terms <- terms(formula, "offset", data=data) 
    fac <- attr(terms, "factors")
  
    needed <- rownames(fac)[!rownames(fac) %in% colnames(fac) & rowSums(fac) > 0]
    if(length(needed) == 0) break
    formula <- addTermsToFormula(needed)
  }
  
  ## Extract offset term from formula. E.g. ~ a1 + offset(a2) -> offsetFac=='a2'.
  offset = attr(terms, "specials")$offset
  if (!is.null(offset)) {
    offsetFac = gsub("offset(", "", rownames(fac)[offset], fixed=TRUE)
    offsetFac = substr(offsetFac, 1, nchar(offsetFac) - 1) # Remove last bracket ')'.
    
    needed <- offsetFac[!offsetFac %in% colnames(fac)]
    formula.withOffset <- addTermsToFormula(offsetFac)
    mm.withOffset <- model.matrix(formula.withOffset, data = data)
    
    offset.values = mm.withOffset[, offsetFac, drop=FALSE]
    
  } else 
    offset.values <- 0
  
  ## Offset.
  ##if(length(offset) != 0)
  ##  offset.values <- offset
  ##else
  ##  offset.values <- 0
  
  ## Construct 3-col response:
  resp <- model.extract(model.frame(formula, data=data), "response")
  if(ncol(resp) == 2)
    resp <- cbind(start=rep(0, nrow(resp)), resp)
  
  ## Sort by STOPtime and Event.
  if(sort) {
    sort <- order(resp[, 2],  -resp[, 3])
    data <- data[sort, , drop=FALSE]
    resp <- resp[sort, ]
  }
  
  ## Model-Matrix.
  mm <- model.matrix(formula, data = data)
  ## Exclude intercept.
  mm1 <- mm[, -1, drop=FALSE]	
  cn = colnames(mm1)
  stopName <- sapply(rownames(fac)[1], f)[2] # Name of stoptime.
  
  ## Extract factors. E.g. for 'x1:log(time)' -> x1, log(time).
  facs = sapply(cn, f2)
  nFacs = unlist(lapply(facs, length))
  
  ## Extract tokens. E.g. for 'x1:log(time)' -> x1, log, time.
  tokens = lapply(cn, f)
  nTime = unlist(lapply(tokens, function(z) sum(z == stopName)))

  ## TRUE for columns holding interaction between time and non-time.
  isTimeInteraction = (nFacs > 1) & (nTime > 0)
  NTDE <- sum(isTimeInteraction)

  ## TRUE for columns holding simply a time-effect (e.g. log(time)).
  isPureTime = (nFacs == 1) & (nTime == 1)

  timedata <- matrix(0, nrow(data), 0)
  timeind <- c()

  for (i in which(isTimeInteraction)) {
    ## For a:b:log(time):d return [F,F,T,F]:
    hasStopName = unlist(lapply(sapply(facs[[i]], f), function(z) any(z == stopName)))
    if (sum(hasStopName) > 1) {
      browser()
      stop("For each time-by-covariate interaction term, please specify only one function of time.")
    }

    ## For a:b:log(time):d return 'log(time)':
    timeEffectName = facs[[i]][hasStopName]

    ## Return TRUE for column in mm1 holding plain time effect:
    indTimeData = (nFacs == 1) & unlist(lapply(facs, function(z) z[1] == timeEffectName))

    ## Add time-specific data:
    timedata = cbind(timedata, mm1[, indTimeData, drop=FALSE])

    if (sum(!hasStopName) > 1) {
      stop("Each time-by-covariate interaction term should be like: f(time) : g(x).")
    }

    ## For a:b:log(time):d return 'log(time)':
    nontimeEffectName = facs[[i]][!hasStopName]

    ## Return TRUE for column in mm1 holding plain non-time effect:
    indNonTimeData = (nFacs == 1) & unlist(lapply(facs, function(z) z[1] == nontimeEffectName))

    timeind = c(timeind, which(indNonTimeData))

    if (!nontimeEffectName %in% effectsOriginal)
      stop("Main covariate effect of time-by-covariate interaction not included in formula.",
           "Time-by-covariate interactions of rank 1 are currently not allowed")
  }

  ## FALSE for columns holding time-effects or not being in original formula.
  ind <- c(nTime == 0 & colnames(mm1) %in% effectsOriginal,
           rep(TRUE, NTDE))
  if (!any(ind)) {
    stop("Nothing to fit. Maybe there is a time-effect without interaction?")
  }

  covnames <- c(colnames(mm1),
                paste(colnames(timedata), colnames(mm1)[timeind], sep=":")) # Time-interactions.

  ## return object
  list(
    resp=resp,                  # N x 3 - response matrix
    mm1=mm1,                    # model matrix without time effects
    
    NTDE=NTDE,                  # number time dep. effects
    timedata=timedata, 	        # matrix with time functions as columns
    timeind=timeind, 	       	  # indicator of time-dependent effects
    
    ## FP are not handled - these 4 parameters exist only for compatibility.
    NFP=0,                      # number of frac.polys
    fpnames=c(),                # names of the variables used as fractional poly's
    fpind=matrix(NA, 0, 0),     # matrix with frac.polyn. in each row, numbered by 1 to 5
    PTcoefs=c(),                # coefficients of pretransformations
    
    covnames=covnames,          # names of covariates
    ind=ind,			              # vector indicating which terms are really part of the formula
    offset.values=offset.values # offset values
  )
}
