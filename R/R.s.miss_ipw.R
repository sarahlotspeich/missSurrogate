# Inverse probability weighting estimators -- type = "model" or "robust"
R.s.miss_ipw = function(sone, szero, yone, yzero, wone, wzero, type) {
  ## Define non-missingness indicators for the two treatment groups
  mone = as.numeric(!is.na(sone))
  mzero = as.numeric(!is.na(szero))

  ## Define sample sizes
  none = length(yone)
  nzero = length(yzero)

  ## Calculate PTE
  if (type == "robust") {
    delta = mean(yone) - mean(yzero)
    delta_S = delta.s.single.ipw(sone = sone,
                             szero = szero,
                             yone = yone,
                             yzero = yzero,
                             weight.1 = wone,
                             weight.0 = wzero)
    R_S = 1 - delta_S / delta

    ## Return
    list(delta = delta,
         delta.s = delta_S,
         R.s = R_S)
  } else if (type == "model") {
    # Create long version of *observed* data (with missingness)
    long_dat = data.frame(Y = c(yzero, yone), ## primary outcome
                          Z = c(rep(c(0, 1), times = c(nzero, none))), ## treatment group
                          S = c(szero, sone), ## surrogate marker
                          R = as.numeric(!is.na(c(szero, sone))), ## (non-)missingness indicator
                          W = c(wzero, wone)) ## weight
    long_dat = long_dat[order(long_dat$S, decreasing = TRUE), ] ## order to put non-missing first

    ## Model Y ~ S * Z (saturated)
    modYgivSZ = lm(formula = Y ~ Z * S,
                   data = long_dat,
                   weights = W)

    ## Separate coefficients
    beta0 = as.numeric(modYgivSZ$coefficients["(Intercept)"])
    beta1 = as.numeric(modYgivSZ$coefficients["Z"])
    beta2 = as.numeric(modYgivSZ$coefficients["S"])
    beta3 = as.numeric(modYgivSZ$coefficients["Z:S"])

    ## Model S ~ Z
    modSgivZ = lm(formula = S ~ Z,
                  data = long_dat,
                  weights = W)

    ## Separate coefficients
    alpha0 = as.numeric(modSgivZ$coefficients["(Intercept)"])
    alpha1 = alpha0 + as.numeric(modSgivZ$coefficients["Z"])

    ## Construct percent treatment effect explained
    delta = beta1 + (beta2 + beta3) * alpha1 - beta2 * alpha0
    delta_S = beta1 + beta3 * alpha0
    R_S = 1 - delta_S / delta

    ## Return
    list(delta = as.numeric(delta),
         delta.s = as.numeric(delta_S),
         R.s = as.numeric(R_S))
  }
}
