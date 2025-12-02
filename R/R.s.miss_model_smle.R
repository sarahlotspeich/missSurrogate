# Sieve maximum likelihood estimator -- Only type = "model"
R.s.miss_model_smle <- function(sone, szero, yone, yzero,
                                nonparam = TRUE,
                                conv_res = NULL,
                                max.it = 1e4,
                                tol = 1e-3,
                                full_output = FALSE) {
  # sizes
  N0 <- length(szero)
  N1 <- length(sone)
  n  <- N0 + N1

  # Build long observed data in natural order; create ID column
  long_dat <- data.frame(
    Y = c(yzero, yone),
    Z = c(rep(0, N0), rep(1, N1)),
    S = c(szero, sone),
    R = as.numeric(!is.na(c(szero, sone))),
    ID = seq_len(n)
  )

  # Separate observed and missing
  long_nonmiss <- long_dat[long_dat$R == 1, , drop = FALSE]
  long_miss    <- long_dat[long_dat$R == 0, , drop = FALSE]

  # Support of observed S values by arm
  S_z0 <- sort(unique(na.omit(long_dat$S[long_dat$Z == 0])))
  S_z1 <- sort(unique(na.omit(long_dat$S[long_dat$Z == 1])))
  m_z0 <- length(S_z0)
  m_z1 <- length(S_z1)

  # Initialize parameters (make prev_p numeric vectors)
  prev_beta  <- rep(0, 4)   # intercept, Z, S, Z:S
  prev_sigma <- 0.1
  if (nonparam) {
    prev_p_z0 <- if (m_z0 > 0) as.numeric(rep(1 / m_z0, m_z0)) else numeric(0)
    prev_p_z1 <- if (m_z1 > 0) as.numeric(rep(1 / m_z1, m_z1)) else numeric(0)
  } else {
    prev_gamma <- c(0, 0)
    prev_eta   <- 0.1
  }

  # Use conv_res if provided
  if (!is.null(conv_res)) {
    if (!is.null(conv_res$betas)) prev_beta  <- conv_res$betas
    if (!is.null(conv_res$sigma)) prev_sigma <- conv_res$sigma
    if (nonparam && !is.null(conv_res$p0) && !is.null(conv_res$p1)) {
      prev_p_z0 <- as.numeric(conv_res$p0)
      prev_p_z1 <- as.numeric(conv_res$p1)
    } else if (!nonparam && !is.null(conv_res$gamma) && !is.null(conv_res$eta)) {
      prev_gamma <- conv_res$gamma
      prev_eta   <- conv_res$eta
    }
  }

  # helper for sigma
  calc_sigma <- function(lmobj) {
    if (is.null(lmobj)) return(NA_real_)
    res <- lmobj$residuals
    rdf <- lmobj$df.residual
    if (is.null(rdf) || rdf <= 0) return(NA_real_)
    sqrt(sum(res^2, na.rm = TRUE) / rdf)
  }

  converged <- FALSE
  it <- 1L

  while (!converged && it <= max.it) {
    # ---------- E-step ----------
    cd_nonmiss <- long_nonmiss
    # expand each missing subject to all candidate S in that arm
    expand_rows <- function(df, Svals) {
      if (nrow(df) == 0 || length(Svals) == 0) return(df[FALSE, , drop = FALSE])
      out_list <- lapply(seq_len(nrow(df)), function(i) {
        rowi <- df[i, , drop = FALSE]
        newrows <- rowi[rep(1, length(Svals)), , drop = FALSE]
        newrows$S <- Svals
        newrows
      })
      do.call(rbind, out_list)
    }
    cd_miss_z0 <- expand_rows(long_miss[long_miss$Z == 0, , drop = FALSE], S_z0)
    cd_miss_z1 <- expand_rows(long_miss[long_miss$Z == 1, , drop = FALSE], S_z1)
    cd_miss    <- if (nrow(cd_miss_z0) + nrow(cd_miss_z1) == 0) cd_miss_z0 else rbind(cd_miss_z0, cd_miss_z1)
    cd         <- if (nrow(cd_miss) == 0) cd_nonmiss else rbind(cd_nonmiss, cd_miss)

    # build phi_aug: 1's for observed rows then posterior probs for expanded rows
    phi_aug <- rep(1, nrow(cd_nonmiss))
    if (nrow(cd_miss) > 0) {
      # P(Y | S, Z) on expanded rows
      mu_beta_miss <- prev_beta[1] + prev_beta[2] * cd_miss$Z +
        prev_beta[3] * cd_miss$S + prev_beta[4] * cd_miss$S * cd_miss$Z
      pY_miss <- dnorm(cd_miss$Y, mean = mu_beta_miss, sd = prev_sigma)

      # P(S | Z) on expanded rows -> use numeric prev_p vectors robustly
      if (nonparam) {
        pS_miss <- numeric(nrow(cd_miss))
        # fill block for Z==0 expanded rows
        if (nrow(cd_miss_z0) > 0) {
          n_miss_z0_subj <- nrow(long_miss[long_miss$Z == 0, , drop = FALSE])
          # repeat prev_p_z0 for each missing subject in z0 (order matches expand_rows)
          pS_miss[seq_len(nrow(cd_miss_z0))] <- rep(as.numeric(prev_p_z0), times = n_miss_z0_subj)
        }
        if (nrow(cd_miss_z1) > 0) {
          start1 <- if (nrow(cd_miss_z0) > 0) nrow(cd_miss_z0) + 1 else 1
          n_miss_z1_subj <- nrow(long_miss[long_miss$Z == 1, , drop = FALSE])
          if (n_miss_z1_subj > 0) {
            pS_miss[start1:(start1 + nrow(cd_miss_z1) - 1)] <- rep(as.numeric(prev_p_z1), times = n_miss_z1_subj)
          }
        }
      } else {
        mu_gamma_miss <- prev_gamma[1] + prev_gamma[2] * cd_miss$Z
        pS_miss <- dnorm(cd_miss$S, mean = mu_gamma_miss, sd = prev_eta)
      }

      phi_num <- pY_miss * pS_miss
      # group-sum phi_num by original subject ID
      denom_byID <- rowsum(phi_num, group = cd_miss$ID, reorder = FALSE)
      denom_rep <- as.numeric(denom_byID[as.character(cd_miss$ID)])
      denom_rep[is.na(denom_rep)] <- 1
      denom_rep[denom_rep == 0] <- 1
      phi_expanded <- phi_num / denom_rep
      phi_aug <- c(phi_aug, phi_expanded)
    }

    if (length(phi_aug) != nrow(cd)) stop("phi_aug length mismatch with cd rows")

    # ---------- M-step ----------
    new_fit <- tryCatch(lm(Y ~ Z * S, data = cd, weights = phi_aug), error = function(e) NULL)
    if (is.null(new_fit)) {
      if (full_output) {
        return(list(delta = NA, delta.s = NA, R.s = NA, betas = rep(NA,4), sigma = NA, p0 = NA, p1 = NA, alphas = c(NA,NA)))
      } else {
        return(list(delta = NA, delta.s = NA, R.s = NA))
      }
    }
    new_beta <- coef(new_fit)
    # ensure length 4 and order by "(Intercept)","Z","S","Z:S"
    all_names <- c("(Intercept)", "Z", "S", "Z:S")
    tmp <- setNames(rep(0, 4), all_names)
    tmp[names(new_beta)] <- new_beta
    new_beta <- as.numeric(tmp)
    new_sigma <- calc_sigma(new_fit)

    # Update p_{kz} robustly: sum phi_aug over rows with Z=z and S == s_k
    if (nonparam) {
      # Z == 0
      idx_z0 <- which(cd$Z == 0)
      if (length(S_z0) == 0) {
        new_p_z0 <- numeric(0)
      } else {
        sums_z0 <- sapply(S_z0, function(sk) {
          sum(phi_aug[idx_z0][cd$S[idx_z0] == sk], na.rm = TRUE)
        })
        if (sum(sums_z0) == 0) sums_z0 <- rep(1/length(sums_z0), length(sums_z0))
        new_p_z0 <- as.numeric(sums_z0 / sum(sums_z0))
      }

      # Z == 1
      idx_z1 <- which(cd$Z == 1)
      if (length(S_z1) == 0) {
        new_p_z1 <- numeric(0)
      } else {
        sums_z1 <- sapply(S_z1, function(sk) {
          sum(phi_aug[idx_z1][cd$S[idx_z1] == sk], na.rm = TRUE)
        })
        if (sum(sums_z1) == 0) sums_z1 <- rep(1/length(sums_z1), length(sums_z1))
        new_p_z1 <- as.numeric(sums_z1 / sum(sums_z1))
      }
    } else {
      sfit <- tryCatch(lm(S ~ Z, data = cd, weights = phi_aug), error = function(e) NULL)
      if (is.null(sfit)) {
        new_gamma <- prev_gamma; new_eta <- prev_eta
      } else {
        new_gamma <- coef(sfit); if (length(new_gamma) < 2) new_gamma <- c(new_gamma, 0)[1:2]
        new_eta <- calc_sigma(sfit)
      }
    }

    # ---------- convergence ----------
    beta_conv  <- all(!is.na(prev_beta)) && all(!is.na(new_beta)) && all(abs(prev_beta - new_beta) < tol)
    sigma_conv <- !is.na(prev_sigma) && !is.na(new_sigma) && (abs(prev_sigma - new_sigma) < tol)
    if (nonparam) {
      # check shape and non-NA
      ok0 <- length(prev_p_z0) == length(new_p_z0) && length(new_p_z0) > 0
      ok1 <- length(prev_p_z1) == length(new_p_z1) && length(new_p_z1) > 0
      if (ok0 && ok1) {
        p_conv <- all(abs(prev_p_z0 - new_p_z0) < tol) && all(abs(prev_p_z1 - new_p_z1) < tol)
      } else {
        # if no missing subjects, treat p_conv as TRUE (nothing to update)
        if (nrow(long_miss) == 0) p_conv <- TRUE else p_conv <- FALSE
      }
      gamma_conv <- TRUE; eta_conv <- TRUE
    } else {
      p_conv <- TRUE
      gamma_conv <- all(!is.na(prev_gamma)) && all(!is.na(new_gamma)) && all(abs(prev_gamma - new_gamma) < tol)
      eta_conv   <- !is.na(prev_eta) && !is.na(new_eta) && (abs(prev_eta - new_eta) < tol)
    }

    converged <- isTRUE(all(c(beta_conv, sigma_conv, p_conv, gamma_conv, eta_conv)))

    # update prev for next iter
    prev_beta  <- new_beta; prev_sigma <- new_sigma
    if (nonparam) { prev_p_z0 <- new_p_z0; prev_p_z1 <- new_p_z1 } else { prev_gamma <- new_gamma; prev_eta <- new_eta }
    it <- it + 1L
  } # end EM loop

  if (!converged) {
    if (full_output) {
      return(list(delta = NA, delta.s = NA, R.s = NA, betas = rep(NA,4), sigma = NA, p0 = NA, p1 = NA, alphas = c(NA,NA)))
    } else {
      return(list(delta = NA, delta.s = NA, R.s = NA))
    }
  }

  # compute alpha0, alpha1 and estimands
  beta0 <- prev_beta[1]; beta1 <- prev_beta[2]; beta2 <- prev_beta[3]; beta3 <- prev_beta[4]
  if (nonparam) {
    alpha0 <- if (length(S_z0) > 0 && length(prev_p_z0) > 0) sum(S_z0 * as.numeric(prev_p_z0)) else NA_real_
    alpha1 <- if (length(S_z1) > 0 && length(prev_p_z1) > 0) sum(S_z1 * as.numeric(prev_p_z1)) else NA_real_
  } else {
    alpha0 <- prev_gamma[1]; alpha1 <- alpha0 + prev_gamma[2]
  }

  delta   <- beta1 + (beta2 + beta3) * alpha1 - beta2 * alpha0
  delta.s <- beta1 + beta3 * alpha0
  R_S     <- 1 - delta.s / delta


  if (full_output) {
    return(list(delta = delta, delta.s = delta.s, R.s = R_S,
                betas = prev_beta, sigma = prev_sigma, p0 = prev_p_z0, p1 = prev_p_z1, alphas = c(alpha0, alpha1)))
  } else {
    return(list(delta = delta, delta.s = delta.s, R.s = R_S))
  }
}
