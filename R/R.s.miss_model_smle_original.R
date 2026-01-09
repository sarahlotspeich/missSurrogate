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
  long_dat_z0 = long.dat[long.dat$Z == 0, ] ## control group
  long_dat_z1 = long.dat[long.dat$Z == 1, ] ## treatment group

  # Create initial values for betas
  ## Linear regression of Y ~ Z + S + X x Z
  prev_beta = rep(0, 4)
  prev_sigma = 0.1

  # Create even longer version of *complete* data (without missingness)
  cd_nonmiss = long.dat[long.dat$R == 1, ] ## patients with non-missing surrogate markers, both treatment groups
  cd_nonmiss_z0 = cd_nonmiss[cd_nonmiss$Z == 0, ] ## separate control group patients
  cd_nonmiss_z1 = cd_nonmiss[cd_nonmiss$Z == 1, ] ## separate treatment group patients

  ## Count and save unique non-missing values of the surrogate
  ### Control group
  S_z0 = unique(cd_nonmiss_z0$S) ## unique values of non-missing surrogates (ordered descendingly)
  m_z0 = length(S_z0) ## number of unique non-missing surrogates
  n0 = nrow(cd_nonmiss_z0) ## number of patients with non-missing surrogates
  miss0 = N0 - n0 ## number of patients with missing surrogates

  ### Treatment group
  S_z1 = unique(cd_nonmiss_z1$S) ## unique values of non-missing surrogates (ordered descendingly)
  m_z1 = length(S_z1) ## number of unique non-missing surrogates
  n1 = nrow(cd_nonmiss_z1) ## number of patients with non-missing surrogates
  miss1 = N1 - n1 ## number of patients with missing surrogates

  ## Create even longer version of *complete* data (filling in missingness)
  ### Complete data for Z = 0
  cd_miss_z0 = long_dat_z0[rep(x = (n0 + 1):N0, each = m_z0), ] ### create m0 copies of each patient with missing surrogate
  cd_miss_z0$S = rep(x = S_z0, times = miss0) ## try out different surrogate values

  ### Complete data for Z = 1
  cd_miss_z1 = long_dat_z1[rep(x = (n1 + 1):N1, each = m_z1), ] ### create m1 copies of each patient with missing surrogate
  cd_miss_z1$S = rep(x = S_z1, times = miss1) ### try out different surrogate values

  ## Combined complete dataset
  cd_miss = rbind(cd_miss_z0, cd_miss_z1) ### only those with missing surrogates
  cd = rbind(cd_nonmiss, cd_miss) ### all patients

  # Conditional distribution of S given Z
  ## Empirical probabilities of S | Z = 0
  prev_p_z0 = p0_z0 = matrix(data = 1 / m_z0,
                             nrow = m_z0,
                             ncol = 1)

  ## Empirical probabilities of S | Z = 1
  prev_p_z1 = p0_z1 = matrix(data = 1 / m_z1,
                             nrow = m_z1,
                             ncol = 1)

  # If converged values supplied, use them as initials
  if (!is.null(conv.res)) {
    ## Linear regression of Y ~ Z + S + X x Z
    prev_beta = conv.res$betas
    prev_sigma = conv.res$sigma
    # Conditional distribution of S given Z
    ## Empirical probabilities of S | Z = 0
    prev_p_z0 = p0_z0 = matrix(data = conv.res$p0,
                               ncol = 1)

    ## Empirical probabilities of S | Z = 1
    prev_p_z1 = p0_z1 = matrix(data = conv.res$p1,
                               ncol = 1)
  }

  # EM Algorithm
  converged = FALSE ## initialize as unconverged
  it = 1 ## initialize iteration counter
  while (!converged & it <= max.it) {
    ## E step ------------------------------------------------------------------
    ### Update the phi_ki = P(S=sk|Zi) for patients w/ missing surrogate -------
    #### Outcome model: P(Y|S,Z) -----------------------------------------------
    ##### Calculate mu = beta0 + beta1X + beta2Z + from previous betas ...
    mu_beta = prev_beta[1] + prev_beta[2] * cd_miss$Z +
      prev_beta[3] * cd_miss$S + prev_beta[4] * cd_miss$S * cd_miss$Z
    ##### and use it with previous sigma to compute P(Y|S,Z) from normal PDF ---
    pYgivSZ = dnorm(x = cd_miss$Y,
                    mean = mu_beta,
                    sd = prev_sigma)
    ############################################################################
    ### Conditional distribution of surrogate given treatment: P(S|Z) ----------
    pSgivZ = c(prev_p_z0[rep(x = 1:m_z0, times = miss0)], #### take from p_k0 if Z = 0
               prev_p_z1[rep(x = 1:m_z1, times = miss1)]) #### take from p_k1 if Z = 1
    ############################################################################
    ## Estimate conditional expectations ---------------------------------------
    ### Update numerator -------------------------------------------------------
    #### P(Y|S,Z)P(S|Z) --------------------------------------------------------
    phi_num = pYgivSZ * pSgivZ
    ### Update denominator -----------------------------------------------------
    #### Sum over P(Y|S,Z)P(S|Z) per patient -----------------------------------
    ##### Control group --------------------------------------------------------
    phi_num_z0 = phi_num[cd_miss$Z == 0]
    phi_denom_z0 = rowsum(x = phi_num_z0,
                          group = cd_miss_z0$ID,
                          reorder = FALSE)
    ##### Treatment group ------------------------------------------------------
    phi_num_z1 = phi_num[cd_miss$Z == 1]
    phi_denom_z1 = rowsum(x = phi_num_z1,
                          group = cd_miss_z1$ID,
                          reorder = FALSE)
    #### Avoid NaN resulting from dividing by 0 --------------------------------
    phi_denom_z0[phi_denom_z0 == 0] = 1
    phi_denom_z1[phi_denom_z1 == 0] = 1
    #### Make long version of phi_denom_z0/1 (same dimension as phi_num_z0/1) --
    long_phi_denom_z0 = rep(x = phi_denom_z0, each = m_z0) #### for Z = 0
    long_phi_denom_z1 = rep(x = phi_denom_z1, each = m_z1) #### for Z = 1
    ### Divide them to get phi = E{I(S=s)|Y,Z} ---------------------------------
    phi_z0 = phi_num_z0 / long_phi_denom_z0
    phi_z1 = phi_num_z1 / long_phi_denom_z1
    phi = c(phi_z0, phi_z1)
    #### Add indicators for non-missing rows -----------------------------------
    phi_aug = c(rep(x = 1, times = (n0 + n1)), phi)
    ## M step
    ### Re-fit the linear regression model
    new_fit = lm(formula = Y ~ Z * S,
                 data = cd,
                 weights = phi_aug)
    new_beta = as.numeric(new_fit$coefficients)
    new_sigma = sigma(new_fit)
    ### Re-estimate empirical probabilities for distribution of S | Z
    #### Among the untreated (Z = 0)
    sum_phi_z0 = rowsum(x = phi_aug[cd$Z == 0], #### Sum over i = 1, ..., n of phi-hats...
                        group = cd$S[cd$Z == 0], #### for each k = 1, ..., m ...
                        reorder = FALSE) #### and keep them in the original order
    lambda_z0 = sum(sum_phi_z0) #### Sum over k = 1, ..., m for constraint
    new_p_z0 = sum_phi_z0 / lambda_z0 #### updated probabilities
    #### Among the treated (Z = 1)
    sum_phi_z1 = rowsum(x = phi_aug[cd$Z == 1], #### Sum over i = 1, ..., n of phi-hats...
                        group = cd$S[cd$Z == 1], #### for each k = 1, ..., m ...
                        reorder = FALSE) #### and keep them in the original order
    lambda_z1 = sum(sum_phi_z1) #### Sum over k = 1, ..., m for constraint
    new_p_z1 = sum_phi_z1 / lambda_z1 #### updated probabilities
    ## Check for convergence
    beta_conv = !any(abs(prev_beta - new_beta) > tol)
    sigma_conv = !(abs(prev_sigma - new_sigma) > tol)
    p_conv = c(!any(abs(prev_p_z0 - new_p_z0) > tol),
               !any(abs(prev_p_z1 - new_p_z1) > tol))
    if (mean(c(beta_conv, sigma_conv, p_conv)) == 1) {
      converged = TRUE ### Success!
    }

    ## If not converged, prepare move on to next iteration
    it = it + 1
    prev_beta = new_beta
    prev_sigma = new_sigma
    prev_p_z0 = new_p_z0
    prev_p_z1 = new_p_z1
  }

  # Return percent of treatment effect explained
  if (converged) {
    ## Define linear regression coefficients
    beta0 = new_beta[1]
    beta1 = new_beta[2]
    beta2 = new_beta[3]
    beta3 = new_beta[4]

    ## Define conditional mean coefficients
    alpha0 = sum(S_z0 * new_p_z0) ### E(S|Z=0)
    ### Complete case mean: mean(long.dat$S[long.dat$Z == 0], na.rm = TRUE)
    alpha1 = sum(S_z1 * new_p_z1) ### E(S|Z=1)
    ### Complete case mean: mean(long.dat$S[long.dat$Z == 1], na.rm = TRUE)

    ## Construct percent treatment effect explained
    delta = beta1 + (beta2 + beta3) * alpha1 - beta2 * alpha0
    delta_S = beta1 + beta3 * alpha0
    R_S = 1 - delta_S / delta

    ## Return
    if (full.output) {
      list(delta = delta,
           delta.s = delta_S,
           R.s = R_S,
           betas = new_beta,
           sigma = new_sigma,
           p0 = new_p_z0,
           p1 = new_p_z1,
           alphas = c(alpha0, alpha1))
    } else {
      list(delta = delta,
           delta.s = delta_S,
           R.s = R_S)
    }
  } else {
    ## Return
    if (full.output) {
      list(delta = NA,
           delta.s = NA,
           R.s = NA,
           betas = rep(NA, length(new_beta)),
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
