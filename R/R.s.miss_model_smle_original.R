# helper for sigma
calc.sigma <- function(lmobj) {
  if (is.null(lmobj)) return(NA_real_)
  res <- lmobj$residuals
  rdf <- lmobj$df.residual
  if (is.null(rdf) || rdf <= 0) return(NA_real_)
  sqrt(sum(res^2, na.rm = TRUE) / rdf)
}

# Sieve maximum likelihood estimator -- Only type = "model"
R.s.miss_model_smle_original = function(sone, szero, yone, yzero, conv.res, max.it = 1E4, tol = 1E-3, full.output = FALSE) {
  # Save useful constants
  N0 = length(szero) ## number in control group
  N1 = length(sone) ## number in treatment group
  N = N0 + N1 ## total sample size

  # Create long version of *observed* data (with missingness)
  long.dat = data.frame(Y = c(yzero, yone), ## primary outcome
                        Z = c(rep(c(0, 1), times = c(N0, N1))), ## treatment group
                        S = c(szero, sone), ## surrogate marker
                        R = as.numeric(!is.na(c(szero, sone)))) ## (non-)missingness indicator
  long.dat = long.dat[order(long.dat$S, decreasing = TRUE), ] ## order to put non-missing first
  long.dat$ID = 1:N ## create an ID

  ## Split the observed data by treatment group (with missingness)
  long.dat.z0 = long.dat[long.dat$Z == 0, ] ## control group
  long.dat.z1 = long.dat[long.dat$Z == 1, ] ## treatment group

  # Create initial values for betas
  ## Linear regression of Y ~ Z + S + X x Z
  prev.beta = rep(0, 4)
  prev.sigma = 0.1

  # Create even longer version of *complete* data (without missingness)
  cd.nonmiss = long.dat[long.dat$R == 1, ] ## patients with non-missing surrogate markers, both treatment groups
  cd.nonmiss.z0 = cd.nonmiss[cd.nonmiss$Z == 0, ] ## separate control group patients
  cd.nonmiss.z1 = cd.nonmiss[cd.nonmiss$Z == 1, ] ## separate treatment group patients

  ## Count and save unique non-missing values of the surrogate
  ### Control group
  S.z0 = unique(cd.nonmiss.z0$S) ## unique values of non-missing surrogates (ordered descendingly)
  m.z0 = length(S.z0) ## number of unique non-missing surrogates
  n0 = nrow(cd.nonmiss.z0) ## number of patients with non-missing surrogates
  miss0 = N0 - n0 ## number of patients with missing surrogates

  ### Treatment group
  S.z1 = unique(cd.nonmiss.z1$S) ## unique values of non-missing surrogates (ordered descendingly)
  m.z1 = length(S.z1) ## number of unique non-missing surrogates
  n1 = nrow(cd.nonmiss.z1) ## number of patients with non-missing surrogates
  miss1 = N1 - n1 ## number of patients with missing surrogates

  ## Create even longer version of *complete* data (filling in missingness)
  ### Complete data for Z = 0
  cd.miss.z0 = long.dat.z0[rep(x = (n0 + 1):N0, each = m.z0), ] ### create m0 copies of each patient with missing surrogate
  cd.miss.z0$S = rep(x = S.z0, times = miss0) ## try out different surrogate values

  ### Complete data for Z = 1
  cd.miss.z1 = long.dat.z1[rep(x = (n1 + 1):N1, each = m.z1), ] ### create m1 copies of each patient with missing surrogate
  cd.miss.z1$S = rep(x = S.z1, times = miss1) ### try out different surrogate values

  ## Combined complete dataset
  cd.miss = rbind(cd.miss.z0, cd.miss.z1) ### only those with missing surrogates
  cd = rbind(cd.nonmiss, cd.miss) ### all patients

  # Conditional distribution of S given Z
  ## Empirical probabilities of S | Z = 0
  prev.p.z0 = p0.z0 = matrix(data = 1 / m.z0,
                             nrow = m.z0,
                             ncol = 1)

  ## Empirical probabilities of S | Z = 1
  prev.p.z1 = p0.z1 = matrix(data = 1 / m.z1,
                             nrow = m.z1,
                             ncol = 1)

  # If converged values supplied, use them as initials
  if (!is.null(conv.res)) {
    ## Linear regression of Y ~ Z + S + X x Z
    prev.beta = conv.res$betas
    prev.sigma = conv.res$sigma
    # Conditional distribution of S given Z
    ## Empirical probabilities of S | Z = 0
    prev.p.z0 = p0.z0 = matrix(data = conv.res$p0,
                               ncol = 1)

    ## Empirical probabilities of S | Z = 1
    prev.p.z1 = p0.z1 = matrix(data = conv.res$p1,
                               ncol = 1)
  }

  # EM Algorithm
  converged = FALSE ## initialize as unconverged
  it = 1 ## initialize iteration counter
  while (!converged & it <= max.it) {
    ## E step ------------------------------------------------------------------
    ### Update the phi.ki = P(S=sk|Zi) for patients w/ missing surrogate -------
    #### Outcome model: P(Y|S,Z) -----------------------------------------------
    ##### Calculate mu = beta0 + beta1X + beta2Z + from previous betas ...
    mu.beta = prev.beta[1] + prev.beta[2] * cd.miss$Z +
      prev.beta[3] * cd.miss$S + prev.beta[4] * cd.miss$S * cd.miss$Z
    ##### and use it with previous sigma to compute P(Y|S,Z) from normal PDF ---
    pYgivSZ = dnorm(x = cd.miss$Y,
                    mean = mu.beta,
                    sd = prev.sigma)
    ############################################################################
    ### Conditional distribution of surrogate given treatment: P(S|Z) ----------
    pSgivZ = c(prev.p.z0[rep(x = 1:m.z0, times = miss0)], #### take from p.k0 if Z = 0
               prev.p.z1[rep(x = 1:m.z1, times = miss1)]) #### take from p.k1 if Z = 1
    ############################################################################
    ## Estimate conditional expectations ---------------------------------------
    ### Update numerator -------------------------------------------------------
    #### P(Y|S,Z)P(S|Z) --------------------------------------------------------
    phi.num = pYgivSZ * pSgivZ
    ### Update denominator -----------------------------------------------------
    #### Sum over P(Y|S,Z)P(S|Z) per patient -----------------------------------
    ##### Control group --------------------------------------------------------
    phi.num.z0 = phi.num[cd.miss$Z == 0]
    phi.denom.z0 = rowsum(x = phi.num.z0,
                          group = cd.miss.z0$ID,
                          reorder = FALSE)
    ##### Treatment group ------------------------------------------------------
    phi.num.z1 = phi.num[cd.miss$Z == 1]
    phi.denom.z1 = rowsum(x = phi.num.z1,
                          group = cd.miss.z1$ID,
                          reorder = FALSE)
    #### Avoid NaN resulting from dividing by 0 --------------------------------
    phi.denom.z0[phi.denom.z0 == 0] = 1
    phi.denom.z1[phi.denom.z1 == 0] = 1
    #### Make long version of phi.denom.z0/1 (same dimension as phi.num.z0/1) --
    long.phi.denom.z0 = rep(x = phi.denom.z0, each = m.z0) #### for Z = 0
    long.phi.denom.z1 = rep(x = phi.denom.z1, each = m.z1) #### for Z = 1
    ### Divide them to get phi = E{I(S=s)|Y,Z} ---------------------------------
    phi.z0 = phi.num.z0 / long.phi.denom.z0
    phi.z1 = phi.num.z1 / long.phi.denom.z1
    phi = c(phi.z0, phi.z1)
    #### Add indicators for non-missing rows -----------------------------------
    phi.aug = c(rep(x = 1, times = (n0 + n1)), phi)
    ## M step
    ### Re-fit the linear regression model
    new.fit = lm(formula = Y ~ Z * S,
                 data = cd,
                 weights = phi.aug)
    new.beta = as.numeric(new.fit$coefficients)
    new.sigma = calc.sigma(new.fit) # sigma(new.fit)
    ### Re-estimate empirical probabilities for distribution of S | Z
    #### Among the untreated (Z = 0)
    sum.phi.z0 = rowsum(x = phi.aug[cd$Z == 0], #### Sum over i = 1, ..., n of phi-hats...
                        group = cd$S[cd$Z == 0], #### for each k = 1, ..., m ...
                        reorder = FALSE) #### and keep them in the original order
    lambda.z0 = sum(sum.phi.z0) #### Sum over k = 1, ..., m for constraint
    new.p.z0 = sum.phi.z0 / lambda.z0 #### updated probabilities
    #### Among the treated (Z = 1)
    sum.phi.z1 = rowsum(x = phi.aug[cd$Z == 1], #### Sum over i = 1, ..., n of phi-hats...
                        group = cd$S[cd$Z == 1], #### for each k = 1, ..., m ...
                        reorder = FALSE) #### and keep them in the original order
    lambda.z1 = sum(sum.phi.z1) #### Sum over k = 1, ..., m for constraint
    new.p.z1 = sum.phi.z1 / lambda.z1 #### updated probabilities
    ## Check for convergence
    beta.conv = !any(abs(prev.beta - new.beta) > tol)
    sigma.conv = !(abs(prev.sigma - new.sigma) > tol)
    p.conv = c(!any(abs(prev.p.z0 - new.p.z0) > tol),
               !any(abs(prev.p.z1 - new.p.z1) > tol))
    if (mean(c(beta.conv, sigma.conv, p.conv)) == 1) {
      converged = TRUE ### Success!
    }

    ## If not converged, prepare move on to next iteration
    it = it + 1
    prev.beta = new.beta
    prev.sigma = new.sigma
    prev.p.z0 = new.p.z0
    prev.p.z1 = new.p.z1
  }

  # Return percent of treatment effect explained
  if (converged) {
    ## Define linear regression coefficients
    beta0 = new.beta[1]
    beta1 = new.beta[2]
    beta2 = new.beta[3]
    beta3 = new.beta[4]

    ## Define conditional mean coefficients
    alpha0 = sum(S.z0 * new.p.z0) ### E(S|Z=0)
    ### Complete case mean: mean(long.dat$S[long.dat$Z == 0], na.rm = TRUE)
    alpha1 = sum(S.z1 * new.p.z1) ### E(S|Z=1)
    ### Complete case mean: mean(long.dat$S[long.dat$Z == 1], na.rm = TRUE)

    ## Construct percent treatment effect explained
    delta = beta1 + (beta2 + beta3) * alpha1 - beta2 * alpha0
    delta.S = beta1 + beta3 * alpha0
    R.S = 1 - delta.S / delta

    ## Return
    if (full.output) {
      list(delta = delta,
           delta.s = delta.S,
           R.s = R.S,
           betas = new.beta,
           sigma = new.sigma,
           p0 = new.p.z0,
           p1 = new.p.z1,
           alphas = c(alpha0, alpha1))
    } else {
      list(delta = delta,
           delta.s = delta.S,
           R.s = R.S)
    }
  } else {
    ## Return
    if (full.output) {
      list(delta = NA,
           delta.s = NA,
           R.s = NA,
           betas = rep(NA, length(new.beta)),
           sigma = NA,
           p0 = NA,
           p1 = NA,
           alphas = rep(NA, 2))
    } else {
      list(delta = NA,
           delta.s = NA,
           R.s = NA)
    }
  }
}
