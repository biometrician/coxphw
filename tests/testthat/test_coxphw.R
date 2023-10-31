library(coxphw)

context("Test if the 'coxphw' function gives correct results")

skip_on_cran()

### COMMENTS & TO-DO ###########################################################
# predict, other examples

### EXAMPLES ###################################################################
### * example from coxphw.Rd file ##############################################
utils::data("gastric")
fit1 <- coxphw(Surv(time, status) ~ radiation, data = gastric, template = "AHR")
#summary(fit1)

# one time deal (first check the fit1 object and then some utility functions)
#names(fit1); fit1$coefficients; fit1$dfbeta.resid; fit1$alpha; fit1$var; fit1$df; fit1$iter; fit1$method.ties; fit1$n
#fit1$y; fit1$formula; fit1$fpind; fit1$PTcoefs; fit1$ind; fit1$cov.j; fit1$cov.lw; fit1$cov.ls; fit1$cov.method; fit1$w.matrix
#fit1$caseweights; fit1$Wald; fit1$means; fit1$linear.predictors; fit1$ci.lower; fit1$ci.upper; fit1$prob; fit1$offset.values
#fit1$dataline; fit1$x; fit1$template; fit1$betafix

results_fit1_1 <- list(names = c("coefficients", "dfbeta.resid", "alpha", "var", "df", "iter", "method.ties", "n", "y", "formula",
                                 "fpind", "PTcoefs", "ind", "cov.j", "cov.lw", "cov.ls", "cov.method", "w.matrix", "caseweights",
                                 "Wald", "means", "linear.predictors", "ci.lower", "ci.upper", "prob", "offset.values", "dataline",
                                 "x", "template", "betafix", "call"),
                       coefficients = c(radiation = 0.4625051),
                       dfbeta_dim = c(90, 1),
                       dfbeta_resid = c(-0.05440737, 0.03276852, 0.03238335, 0.03199270, 0.03159641, 0.03119431, -0.04604401, 0.02942577,
                                        0.02901710,  0.02860210),
                       alpha = 0.05,
                       var = matrix(0.05699834, 1, 1, dimnames = list("radiation", "radiation")),
                       df = 1,
                       iter = 3,
                       method.ties = "breslow",
                       n = 90,
                       n_dim = c(90, 3),
                       n = matrix(c(rep(0, 10), 1, 17, 42, 44, 48, 60, 63, 72, 74, 95, rep(1, 10)), 10, 3,
                                  dimnames = list(1:10, c("start", "time", "status"))),
                       formula = "Surv(time, status) ~ radiation",
                       fpind = "structure(logical(0), dim = c(0L, 0L))",
                       PTcoefs = NULL,
                       ind = TRUE,
                       cov.j = NULL,
                       cov.lw = matrix(c(0.05699834), 1, 1, dimnames = list("radiation", "radiation")),
                       cov.ls = matrix(c(0.06061221), 1, 1, dimnames = list("radiation", "radiation")),
                       cov.method = "Lin-Wei",
                       wmatrix_dim = c(79, 4),
                       w.matrix = matrix(c(1, 17, 42, 44, 48, 60, 63, 72, 74, 95,  1, 0.9888889, 0.9777778, 0.9666667, 0.9555556, 0.9444444,
                                           0.9333333, 0.9222222, 0.9111111, 0.9000000, 1,  1,  1,  1,  1,  1,  1,  1,  1,  1, 1.7394781,
                                           1.7201506, 1.7008230, 1.6814955, 1.6621680, 1.6428404, 1.6235129, 1.6041854, 1.5848578,  1.5655303),
                                         10, 4, dimnames = list(1:10, c("time", "w.raw", "w.obskm", "w"))),
                       caseweights = rep(1, 90),
                       wald = 3.752934,
                       means = c(radiation = 0.5),
                       lp = c(-0.2312526, rep(0.2312526, 5), -0.2312526, rep(0.2312526, 3)),
                       cil = c(radiation = 0.9945916),
                       ciu = c(radiation = 2.535608),
                       prob = c(radiation = 0.05271492),
                       offset = 0,
                       dataline = matrix(c(1, 1, 0), 1, 3, dimnames = list(1, c("time", "status", "radiation"))),
                       x_dim = c(90, 1),
                       x = matrix(c(0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1,
                                    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1,
                                    0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1), 90, 1, dimnames = list(1:90, "radiation")),
                       template = "AHR",
                       betafix = NULL)


#coef(fit1);  concord(fit1, digits = 4); confint(fit1); wald(fit1$coefficients, fit1$var)
results_fit1_2 <- list(coef = c(radiation = 0.4625051),
                     concord = matrix(c(0.6136, 0.4986, 0.7172), 1, 3, dimnames = list("radiation", c("concordance prob.", "lower 0.95", "upper 0.95"))),
                     confint = matrix(c(-0.005423054, 0.9304333), 1, 2, dimnames = list("radiation", c("2.5 %",  "97.5 %"))),
                     wald = c(chi2 = 3.75293360, df = 1, p = 0.05271492))


# let's test
test_that("check the OBJECT fit1 <- coxphw(Surv(time, status) ~ radiation, data = gastric, template = AHR)", {
  expect_equal(names(fit1),                    results_fit1_1[[1]],  tolerance = 1e-07, scale = 1)
  expect_equal(fit1$coefficients,              results_fit1_1[[2]],  tolerance = 1e-07, scale = 1)
  expect_equal(dim(fit1$dfbeta.resid),         results_fit1_1[[3]],  tolerance = 1e-07, scale = 1)
  expect_equal(fit1$dfbeta.resid[1:10,],       results_fit1_1[[4]],  tolerance = 1e-07, scale = 1)
  expect_equal(fit1$alpha,                     results_fit1_1[[5]],  tolerance = 1e-07, scale = 1)
  expect_equal(fit1$var,                       results_fit1_1[[6]],  tolerance = 1e-07, scale = 1)
  expect_equal(fit1$df,                        results_fit1_1[[7]],  tolerance = 1e-07, scale = 1)
  expect_equal(fit1$iter,                      results_fit1_1[[8]],  tolerance = 1e-07, scale = 1)
  expect_equal(fit1$method.ties,               results_fit1_1[[9]],  tolerance = 1e-07, scale = 1)
  expect_equal(fit1$n,                         results_fit1_1[[10]], tolerance = 1e-07, scale = 1)
  expect_equal(dim(fit1$y),                    results_fit1_1[[11]], tolerance = 1e-07, scale = 1)
  expect_equal(fit1$y[1:10,],                  results_fit1_1[[12]], tolerance = 1e-07, scale = 1)
  expect_equal(deparse(formula(fit1$formula)), results_fit1_1[[13]], tolerance = 1e-07, scale = 1)
  expect_equal(deparse(fit1$fpind),            results_fit1_1[[14]], tolerance = 1e-07, scale = 1)
  expect_equal(fit1$PTcoefs,                   results_fit1_1[[15]], tolerance = 1e-07, scale = 1)
  expect_equal(fit1$ind,                       results_fit1_1[[16]], tolerance = 1e-07, scale = 1)
  expect_equal(fit1$cov.j,                     results_fit1_1[[17]], tolerance = 1e-07, scale = 1)
  expect_equal(fit1$cov.lw,                    results_fit1_1[[18]], tolerance = 1e-07, scale = 1)
  expect_equal(fit1$cov.ls,                    results_fit1_1[[19]], tolerance = 1e-07, scale = 1)
  expect_equal(fit1$cov.method,                results_fit1_1[[20]], tolerance = 1e-07, scale = 1)
  expect_equal(dim(fit1$w.matrix),             results_fit1_1[[21]], tolerance = 1e-07, scale = 1)
  expect_equal(fit1$w.matrix[1:10,],           results_fit1_1[[22]], tolerance = 1e-07, scale = 1)
  expect_equal(fit1$caseweights,               results_fit1_1[[23]], tolerance = 1e-07, scale = 1)
  expect_equal(fit1$Wald,                      results_fit1_1[[24]], tolerance = 1e-06, scale = 1)
  expect_equal(fit1$means,                     results_fit1_1[[25]], tolerance = 1e-07, scale = 1)
  expect_equal(fit1$linear.predictors[1:10],   results_fit1_1[[26]], tolerance = 1e-07, scale = 1)
  expect_equal(fit1$ci.lower,                  results_fit1_1[[27]], tolerance = 1e-07, scale = 1)
  expect_equal(fit1$ci.upper,                  results_fit1_1[[28]], tolerance = 1e-06, scale = 1)
  expect_equal(fit1$prob,                      results_fit1_1[[29]], tolerance = 1e-07, scale = 1)
  expect_equal(fit1$offset.values,             results_fit1_1[[30]], tolerance = 1e-07, scale = 1)
  expect_equal(as.matrix(fit1$dataline),       results_fit1_1[[31]], tolerance = 1e-07, scale = 1)
  expect_equal(dim(fit1$x),                    results_fit1_1[[32]], tolerance = 1e-07, scale = 1)
  expect_equal(fit1$x,                         results_fit1_1[[33]], tolerance = 1e-07, scale = 1)
  expect_equal(fit1$template,                  results_fit1_1[[34]], tolerance = 1e-07, scale = 1)
  expect_equal(fit1$betafix,                   results_fit1_1[[35]], tolerance = 1e-07, scale = 1)
})

test_that("check UTILITY functions fit1 <- coxphw(Surv(time, status) ~ radiation, data = gastric, template = AHR", {
  expect_equal(coef(fit1),                        results_fit1_2[[1]], tolerance = 1e-07, scale = 1)
  expect_equal(concord(fit1, digits = 4),         results_fit1_2[[2]], tolerance = 1e-07, scale = 1)
  expect_equal(confint(fit1),                     results_fit1_2[[3]], tolerance = 1e-07, scale = 1)
  expect_equal(wald(fit1$coefficients, fit1$var), results_fit1_2[[4]], tolerance = 1e-07, scale = 1)
})

rm(list = ls())

