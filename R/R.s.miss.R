#' Surrogate marker validation with missing values
#' This function calculates the proportion of treatment effect on the primary outcome explained by the treatment effect on a surrogate marker, correcting for missing values in the surrogate marker. This function is intended to be used for a fully observed continuous outcome. The user must specify what type of estimation they would like (parametric or nonparametric estimation of the proportion explained, denoted by R) and what estimator they would like (see below for details). Missing values can be handled through inverse probability weighting (IPW) or semiparametric sieve maximum likelihood estimation (SMLE).
#'
#' @param sone numeric vector; surrogate marker for treated observations, assumed to be continuous.
#' @param szero numeric vector; surrogate marker for control observations, assumed to be continuous.
#' @param yone numeric vector; primary outcome for treated observations, assumed to be continuous.
#' @param yzero numeric vector; primary outcome for control observations, assumed to be continuous.
#' @param wone  (only for IPW) numeric vector; probabilities of having non-missing surrogates for treated observations, if \code{NULL} (the default) then SMLE will be used.
#' @param wzero (only for IPW) numeric vector; probabilities of having non-missing surrogates for control observations, if \code{NULL} (the default) then SMLE will be used.
#' @param type string specifying the type of estimation; choices are \code{"robust"} (the default) or \code{"model"}.
#' @param max.it (only for SMLE) scalar; maximum number of iterations allowed in the EM algorithm, the default is \code{1E4}.
#' @param tol (only for SMLE) scalar; tolerance between iterations in the EM algorithm used to define convergence, the default is \code{1E-3}.
#' @param conf.int logical; if \code{TRUE}, variance estimates and 95% confidence intervals are included in output; if \code{FALSE} (the default), only point estimates are.
#' @param ipw.formula (only for IPW) formula; if \code{conf.int = TRUE} and using IPW, the logistic regression formula needed to recalculate the probabilities of non-missing surrogates are needed.
#' @return A list is returned:
#' \item{delta}{Estimated average treatment effect on the primary outcome.}
#' \item{delta.s}{Estimated average treatment effect on the surrogate marker.}
#' \item{R.s}{Estimated proportion of treatment effect on the primary outcome explained by the treatment effect on a surrogate marker.}
#' \item{delta.var}{Estimated variance for \code{delta} (if \code{conf.int = TRUE}).}
#' \item{delta.s.var}{Estimated variance for \code{delta.s} (if \code{conf.int = TRUE}).}
#' \item{R.s.var}{Estimated variance for \code{R.s} (if \code{conf.int = TRUE}).}
#' \item{conf.int.normal.delta}{95\% confidence interval for \code{delta} based on the normal approximation (if \code{conf.int = TRUE}).}
#' \item{conf.int.quantile.delta}{95\% confidence interval for \code{delta} based on quantiles (if \code{conf.int = TRUE}).}
#' \item{conf.int.normal.delta.s}{95\% confidence interval for \code{delta.s} based on the normal approximation (if \code{conf.int = TRUE}).}
#' \item{conf.int.quantile.delta.s}{95\% confidence interval for \code{delta.s} based on quantiles (if \code{conf.int = TRUE}).}
#' \item{conf.int.normal.R.s}{95\% confidence interval for \code{R.s} based on the normal approximation (if \code{conf.int = TRUE}).}
#' \item{conf.int.quantile.R.s}{95\% confidence interval for \code{R.s} based on quantiles (if \code{conf.int = TRUE}).}
#' @export

R.s.miss = function(sone, szero, yone, yzero, wone = NULL, wzero = NULL,
                    type = "robust", max.it = 1E4, tol = 1E-3,
                    conf.int = FALSE, ipw.formula = NULL) {
  # Estimate parameters
  est_res = R.s.miss.estimate(sone = sone, szero = szero, ## surrogates outcomes
                              yone = yone, yzero = yzero, ## primary outcomes
                              wone = wone, wzero = wzero, ## weights (optional)
                              type = type, max.it = max.it, tol = tol) ## other arguments

  # Estimate standard errors
  ## Bootstrap resampling
  if (conf.int) {
    se_res = R.s.miss.se.boot(num_boot = 500, conv_res = est_res,
                              sone = sone, szero = szero, yone = yone, yzero = yzero,
                              type = type, max.it = max.it, tol = tol, ipw.formula = ipw.formula)
  } else {
    ### Otherwise, just return point estimates (without SEs/CIs)
    return(est_res)
  }

  ### Build final return list (based on R.s.estimate syntax)
  res_list = list(
    delta = as.numeric(est_res$delta),
    delta.s = as.numeric(est_res$delta.s),
    R.s = as.numeric(est_res$R.s),
    delta.var = as.numeric(se_res$var_delta),
    delta.s.var = as.numeric(se_res$var_delta.s),
    R.s.var = as.numeric(se_res$var_R.s),
    conf.int.normal.delta = se_res$norm_ci_delta,
    conf.int.quantile.delta = se_res$quant_ci_delta,
    conf.int.normal.delta.s = se_res$norm_ci_delta.s,
    conf.int.quantile.delta.s = se_res$quant_ci_delta.s,
    conf.int.normal.R.s = se_res$norm_ci_R.s,
    conf.int.quantile.R.s = se_res$quant_ci_R.s
  )

  ## And return it
  return(res_list)
}
