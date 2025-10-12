R.s.miss.estimate = function(weight_perturb = NULL, sone, szero, yone, yzero, wone = NULL, wzero = NULL, conv_res = NULL, type, max_it, tol, ipw_formula = NULL) {
  ## Define sample sizes
  none = length(yone)
  nzero = length(yzero)

  ## Define non-missingness indicators for the two treatment groups
  mone = as.numeric(!is.na(sone))
  mzero = as.numeric(!is.na(szero))

  ### If weight_perturb provided, multiply surrogates and outcomes by it
  if (!is.null(weight_perturb)) {
    szero = szero * weight_perturb[1:nzero]
    yzero = yzero * weight_perturb[1:nzero]
    sone = sone * weight_perturb[-c(1:nzero)]
    yone = yone * weight_perturb[-c(1:nzero)]
  }

  # Define TRUE/FALSE use IPW based on non-null weights supplied
  ipw = (!is.null(wone) & !is.null(wzero)) || ## either weights were supplied
    !is.null(ipw_formula) ## or the model formula was

  # Estimate parameters and SEs
  if (ipw) {
    ## If ipw_formula is not NULL, re-calculate weights
    if (!is.null(ipw_formula) & is.null(wone) & is.null(wzero)) {
      ### Define vectors of variables for model (all of them just in case)
      m = c(mone, mzero)
      z = rep(x = c(1, 0), times = c(length(sone), length(szero)))
      s = c(sone, szero)
      y = c(yone, yzero)

      ### Fit the IPW model
      ipw_fit = glm(formula = as.formula(ipw_formula),
                    data = data.frame(m, z, s, y),
                    family = "binomial")

      ### Get estimated weights (probabilities of being non-missing) for each observation
      w = 1 / predict(object = ipw_fit,
                      type = "response")

      ### Split weights into vectors for treatment/control
      wone = w[1:length(sone)]
      wzero = w[-c(1:length(sone))]
    }

    ## Using IPW to handle missing data
    est_res = R.s.miss_ipw(sone = sone, szero = szero, ### surrogates outcomes
                           yone = yone, yzero = yzero, ### primary outcomes
                           wone = wone, wzero = wzero, ### weights (required)
                           type = type) ### type of PTE estimator
  } else if (type == "model") { # Wang & Taylor's approach
    est_res = R.s.miss_model_smle(sone = sone,
                                  szero = szero,
                                  yone = yone,
                                  yzero = yzero,
                                  nonparam = TRUE,
                                  conv_res = conv_res,
                                  max_it = max_it,
                                  tol = tol)
  }

  ### Return point estimates
  return(est_res)
}
