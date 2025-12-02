# Perturbation resampling SEs
R.s.miss.se.pert = function(num_pert, conv_res, sone, szero, yone, yzero,
                            ipw.formula, type, max.it, tol) {
  # Create (n0 + n1) x num_pert matrix of perturbations
  weight_perturb = matrix(rexp(n = num_pert * (length(yone) + length(yzero)), rate = 1),
                          ncol = num_pert)

  # Apply the estimation functions with it
  pert_quant = do.call(what = rbind,
                       args = apply(X = weight_perturb,
                                    MARGIN = 2, ## apply across columns
                                    FUN = R.s.miss.estimate,
                                    sone = sone,
                                    szero = szero,
                                    yone = yone,
                                    yzero = yzero,
                                    ipw.formula = ipw.formula,
                                    type = type,
                                    max.it = max.it,
                                    tol = tol)
  )

  ## Create separate vectors for perturbed quantities
  pert_delta = unlist(pert_quant[, "delta"])
  pert_delta.s = unlist(pert_quant[, "delta.s"])
  pert_R.s = unlist(pert_quant[, "R.s"])

  # Calculate two types of 95% confidence intervals
  ## Normal approximation
  ### Variance estimates for each quantity
  var_delta = var(pert_delta)
  var_delta.s = var(pert_delta.s)
  var_R.s = var(pert_R.s)

  ### Used to compute Wald-type confidence intervals
  norm_ci_delta = conv_res$delta + c(-1.96, 1.96) * sqrt(var_delta)
  norm_ci_delta.s = conv_res$delta.s + c(-1.96, 1.96) * sqrt(var_delta.s)
  norm_ci_R.s = conv_res$R.s + c(-1.96, 1.96) * sqrt(var_R.s)

  ## Quantile-based
  quant_ci_delta = as.vector(quantile(x = pert_delta, probs = c(0.025, 0.975)))
  quant_ci_delta.s = as.vector(quantile(x = pert_delta.s, probs = c(0.025, 0.975)))
  quant_ci_R.s = as.vector(quantile(x = pert_R.s, probs = c(0.025, 0.975)))

  ## Return
  return(
    list(
      var_delta = var_delta,
      var_delta.s = var_delta.s,
      var_R.s = var_R.s,
      norm_ci_delta = norm_ci_delta,
      quant_ci_delta = quant_ci_delta,
      norm_ci_delta.s = norm_ci_delta.s,
      quant_ci_delta.s = quant_ci_delta.s,
      norm_ci_R.s = norm_ci_R.s,
      quant_ci_R.s = quant_ci_R.s
    )
  )
}
