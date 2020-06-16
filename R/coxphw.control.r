coxphw.control <- function(

iter.max       = 200,                  # max. iterations
maxhs          = 5,                    # half steps
xconv          = 1e-4,                 # convergence criterion for standardized parameter estimates
gconv          = 1e-4,                 # convergence criterion for first derivatives of log likelihood
maxstep        = 1,
round.times.to = 0.00001,              # whether times should be rounded (can be safer with FORTRAN interface), set to 0 for no rounding
add.constant   = 0,                    # add 1*add.constant to start times, 2*add.constant to stop times
pc             = TRUE,                 # transforms data by principal components to speed up convergence when evaluating fp models
pc.time        = TRUE,                 # transforms time variables by principal components to speed up convergence in models
normalize      = TRUE
#fp.iter        = 10                    # maximum number of iterations of large <fp> loop
)

{
  # Gather all of the control parameters for coxphw into one spot

  if (iter.max < 0) stop("Invalid value for iter.max")
  if (maxhs < 0) stop("Invalid value for maxhs")

  if (maxstep < 0) stop("Invalid value for maxstep")                         # GG?

  if (xconv <= 0) stop ("Invalid convergence criterion for standardized parameter estimates (xconv)")
  if (gconv <= 0) stop ("Invalid convergence criterion for first derivatives of log likelihood (gconv)")

  if (round.times.to < 0) stop("Invalid round.times.to")

#  if (is.logcial(pc) || is.logical(pc.time) || is.logical(normalize)) stop("Arguments pc, pc.time and normalize must be logical expressions (either TRUE or FALSE)")

#  if (fp.iter < 0) stop("Invalid value for fp.iter")

  list(iter.max = as.integer(iter.max), maxhs = as.integer(maxhs), xconv = xconv, gconv = gconv,
       maxstep = maxstep, round.times.to = round.times.to, add.constant = add.constant, pc = pc,
       pc.time = pc.time, normalize = normalize, fp.iter = 10)
}
