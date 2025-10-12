R.s.miss.se.boot = function(num_boot, conv_res, sone, szero, yone, yzero,
                            type, ipw, smle, max_it, tol, ipw_formula) {
  # Build long dataset: (y, z, s, w, m)
  long_dat = data.frame(y = c(yone, yzero),
                        z = rep(x = c(1, 0),
                                times = c(length(yone), length(yzero))),
                        s = c(sone, szero),
                        m = c(as.numeric(!is.na(sone)), as.numeric(!is.na(szero))))

  # Initialize empty dataframe to hold bootstrapped quantities
  boot_quant = data.frame(delta = rep(NA, num_boot),
                          delta.s = rep(NA, num_boot),
                          R.s = rep(NA, num_boot))

  # Bootstrap resample from long dataset and fit estimator to it
  for (b in 1:num_boot) {
    ## Indices for which rows to resample, preserving the treatment/control split
    ind_b = c(sample(x = which(long_dat$z == 0),
                     size = length(which(long_dat$z == 0)),
                     replace = TRUE),
              sample(x = which(long_dat$z == 1),
                     size = length(which(long_dat$z == 1)),
                     replace = TRUE))

    boot_quant[b, ] = boot_R.s.miss(data = long_dat,
                                    indices = ind_b,
                                    type = type,
                                    ipw_formula = ipw_formula,
                                    conv_res = NULL,
                                    max_it = max_it,
                                    tol = tol)
  }

  # Calculate two types of 95% confidence intervals
  ## Normal approximation
  ### Variance estimates for each quantity
  var_all = apply(X = boot_quant,
                  MARGIN = 2,
                  FUN = var)
  var_delta = var_all[1]
  var_delta.s = var_all[2]
  var_R.s = var_all[3]

  ### Used to compute Wald-type confidence intervals
  norm_ci_delta = conv_res$delta + c(-1.96, 1.96) * sqrt(var_delta)
  norm_ci_delta.s = conv_res$delta.s + c(-1.96, 1.96) * sqrt(var_delta.s)
  norm_ci_R.s = conv_res$R.s + c(-1.96, 1.96) * sqrt(var_R.s)

  ## Quantile-based
  quant_ci_all = apply(X = boot_quant,
                       MARGIN = 2,
                       FUN = function(x) quantile(x = x, probs = c(0.025, 0.975)))
  quant_ci_delta = as.vector(quant_ci_all[, 1])
  quant_ci_delta.s = as.vector(quant_ci_all[, 2])
  quant_ci_R.s = as.vector(quant_ci_all[, 3])

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
