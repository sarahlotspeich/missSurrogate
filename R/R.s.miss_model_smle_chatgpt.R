# Sieve maximum likelihood estimator -- Only type = "model"
R.s.miss_model_smle_chatgpt <- function(sone, szero, yone, yzero,
                                nonparam = TRUE,
                                conv.res = NULL,
                                max.it = 1e4,
                                tol = 1e-3,
                                full.output = FALSE) {
  # sizes
  N0 <- length(szero)
  N1 <- length(sone)
  n  <- N0 + N1

  # Build long observed data in natural order; create ID column
  long.dat <- data.frame(
    Y = c(yzero, yone),
    Z = c(rep(0, N0), rep(1, N1)),
    S = c(szero, sone),
    R = as.numeric(!is.na(c(szero, sone))),
    ID = seq_len(n)
  )

  # Separate observed and missing
  long.nonmiss <- long.dat[long.dat$R == 1, , drop = FALSE]
  long.miss    <- long.dat[long.dat$R == 0, , drop = FALSE]

  # Support of observed S values by arm
  S.z0 <- sort(unique(na.omit(long.dat$S[long.dat$Z == 0])))
  S.z1 <- sort(unique(na.omit(long.dat$S[long.dat$Z == 1])))
  m.z0 <- length(S.z0)
  m.z1 <- length(S.z1)

  # Initialize parameters (make prev.p numeric vectors)
  prev.beta  <- rep(0, 4) ## intercept, Z, S, Z:S
  prev.sigma <- 0.1 ## residual standard deviation
  if (nonparam) {
    prev.p.z0 <- if (m.z0 > 0) as.numeric(rep(1 / m.z0, m.z0)) else numeric(0) ## P(S|Z=0)
    prev.p.z1 <- if (m.z1 > 0) as.numeric(rep(1 / m.z1, m.z1)) else numeric(0) ## P(S|Z=1)
  } else {
    prev.gamma <- c(0, 0) ## intercept, Z
    prev.eta   <- 0.1 ## residual standard deviation
  }

  # Use conv.res if provided
  if (!is.null(conv.res)) {
    if (!is.null(conv.res$betas)) prev.beta  <- conv.res$betas
    if (!is.null(conv.res$sigma)) prev.sigma <- conv.res$sigma
    if (nonparam && !is.null(conv.res$p0) && !is.null(conv.res$p1)) {
      prev.p.z0 <- as.numeric(conv.res$p0)
      prev.p.z1 <- as.numeric(conv.res$p1)
    } else if (!nonparam && !is.null(conv.res$gamma) && !is.null(conv.res$eta)) {
      prev.gamma <- conv.res$gamma
      prev.eta   <- conv.res$eta
    }
  }

  # helper for sigma
  calc.sigma <- function(lmobj) {
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
    cd.nonmiss <- long.nonmiss
    # expand each missing subject to all candidate S in that arm
    expand.rows <- function(df, Svals) {
      if (nrow(df) == 0 || length(Svals) == 0) return(df[FALSE, , drop = FALSE])
      out.list <- lapply(seq_len(nrow(df)), function(i) {
        rowi <- df[i, , drop = FALSE]
        newrows <- rowi[rep(1, length(Svals)), , drop = FALSE]
        newrows$S <- Svals
        newrows
      })
      do.call(rbind, out.list)
    }
    cd.miss.z0 <- expand.rows(long.miss[long.miss$Z == 0, , drop = FALSE], S.z0)
    cd.miss.z1 <- expand.rows(long.miss[long.miss$Z == 1, , drop = FALSE], S.z1)
    cd.miss    <- if (nrow(cd.miss.z0) + nrow(cd.miss.z1) == 0) cd.miss.z0 else rbind(cd.miss.z0, cd.miss.z1)
    cd         <- if (nrow(cd.miss) == 0) cd.nonmiss else rbind(cd.nonmiss, cd.miss)

    # build phi.aug: 1's for observed rows then posterior probs for expanded rows
    phi.aug <- rep(1, nrow(cd.nonmiss))
    if (nrow(cd.miss) > 0) {
      # P(Y | S, Z) on expanded rows
      mu.beta.miss <- prev.beta[1] + prev.beta[2] * cd.miss$Z +
        prev.beta[3] * cd.miss$S + prev.beta[4] * cd.miss$S * cd.miss$Z
      pY.miss <- dnorm(cd.miss$Y, mean = mu.beta.miss, sd = prev.sigma)

      # P(S | Z) on expanded rows -> use numeric prev.p vectors robustly
      if (nonparam) {
        pS.miss <- numeric(nrow(cd.miss))
        # fill block for Z==0 expanded rows
        if (nrow(cd.miss.z0) > 0) {
          n.miss.z0.subj <- nrow(long.miss[long.miss$Z == 0, , drop = FALSE])
          # repeat prev.p.z0 for each missing subject in z0 (order matches expand.rows)
          pS.miss[seq_len(nrow(cd.miss.z0))] <- rep(as.numeric(prev.p.z0), times = n.miss.z0.subj)
        }
        if (nrow(cd.miss.z1) > 0) {
          start1 <- if (nrow(cd.miss.z0) > 0) nrow(cd.miss.z0) + 1 else 1
          n.miss.z1.subj <- nrow(long.miss[long.miss$Z == 1, , drop = FALSE])
          if (n.miss.z1.subj > 0) {
            pS.miss[start1:(start1 + nrow(cd.miss.z1) - 1)] <- rep(as.numeric(prev.p.z1), times = n.miss.z1.subj)
          }
        }
      } else {
        mu.gamma.miss <- prev.gamma[1] + prev.gamma[2] * cd.miss$Z
        pS.miss <- dnorm(cd.miss$S, mean = mu.gamma.miss, sd = prev.eta)
      }

      phi.num <- pY.miss * pS.miss
      # group-sum phi.num by original subject ID
      denom.byID <- rowsum(phi.num, group = cd.miss$ID, reorder = FALSE)
      denom.rep <- as.numeric(denom.byID[as.character(cd.miss$ID)])
      denom.rep[is.na(denom.rep)] <- 1
      denom.rep[denom.rep == 0] <- 1
      phi.expanded <- phi.num / denom.rep
      phi.aug <- c(phi.aug, phi.expanded)
    }

    if (length(phi.aug) != nrow(cd)) stop("phi.aug length mismatch with cd rows")

    # ---------- M-step ----------
    new.fit <- tryCatch(lm(Y ~ Z * S, data = cd, weights = phi.aug), error = function(e) NULL)
    if (is.null(new.fit)) {
      if (full.output) {
        return(list(delta = NA, delta.s = NA, R.s = NA, betas = rep(NA,4), sigma = NA, p0 = NA, p1 = NA, alphas = c(NA,NA)))
      } else {
        return(list(delta = NA, delta.s = NA, R.s = NA))
      }
    }
    new.beta <- coef(new.fit)
    # ensure length 4 and order by "(Intercept)","Z","S","Z:S"
    all.names <- c("(Intercept)", "Z", "S", "Z:S")
    tmp <- setNames(rep(0, 4), all.names)
    tmp[names(new.beta)] <- new.beta
    new.beta <- as.numeric(tmp)
    new.sigma <- calc.sigma(new.fit)

    # Update p.{kz} robustly: sum phi.aug over rows with Z=z and S == s.k
    if (nonparam) {
      # Z == 0
      idx.z0 <- which(cd$Z == 0)
      if (length(S.z0) == 0) {
        new.p.z0 <- numeric(0)
      } else {
        sums.z0 <- sapply(S.z0, function(sk) {
          sum(phi.aug[idx.z0][cd$S[idx.z0] == sk], na.rm = TRUE)
        })
        if (sum(sums.z0) == 0) sums.z0 <- rep(1/length(sums.z0), length(sums.z0))
        new.p.z0 <- as.numeric(sums.z0 / sum(sums.z0))
      }

      # Z == 1
      idx.z1 <- which(cd$Z == 1)
      if (length(S.z1) == 0) {
        new.p.z1 <- numeric(0)
      } else {
        sums.z1 <- sapply(S.z1, function(sk) {
          sum(phi.aug[idx.z1][cd$S[idx.z1] == sk], na.rm = TRUE)
        })
        if (sum(sums.z1) == 0) sums.z1 <- rep(1/length(sums.z1), length(sums.z1))
        new.p.z1 <- as.numeric(sums.z1 / sum(sums.z1))
      }
    } else {
      sfit <- tryCatch(lm(S ~ Z, data = cd, weights = phi.aug), error = function(e) NULL)
      if (is.null(sfit)) {
        new.gamma <- prev.gamma; new.eta <- prev.eta
      } else {
        new.gamma <- coef(sfit); if (length(new.gamma) < 2) new.gamma <- c(new.gamma, 0)[1:2]
        new.eta <- calc.sigma(sfit)
      }
    }

    # ---------- convergence ----------
    beta.conv  <- all(!is.na(prev.beta)) && all(!is.na(new.beta)) && all(abs(prev.beta - new.beta) < tol)
    sigma.conv <- !is.na(prev.sigma) && !is.na(new.sigma) && (abs(prev.sigma - new.sigma) < tol)
    if (nonparam) {
      # check shape and non-NA
      ok0 <- length(prev.p.z0) == length(new.p.z0) && length(new.p.z0) > 0
      ok1 <- length(prev.p.z1) == length(new.p.z1) && length(new.p.z1) > 0
      if (ok0 && ok1) {
        p.conv <- all(abs(prev.p.z0 - new.p.z0) < tol) && all(abs(prev.p.z1 - new.p.z1) < tol)
      } else {
        # if no missing subjects, treat p.conv as TRUE (nothing to update)
        if (nrow(long.miss) == 0) p.conv <- TRUE else p.conv <- FALSE
      }
      gamma.conv <- TRUE; eta.conv <- TRUE
    } else {
      p.conv <- TRUE
      gamma.conv <- all(!is.na(prev.gamma)) && all(!is.na(new.gamma)) && all(abs(prev.gamma - new.gamma) < tol)
      eta.conv   <- !is.na(prev.eta) && !is.na(new.eta) && (abs(prev.eta - new.eta) < tol)
    }

    converged <- isTRUE(all(c(beta.conv, sigma.conv, p.conv, gamma.conv, eta.conv)))

    # update prev for next iter
    prev.beta  <- new.beta; prev.sigma <- new.sigma
    if (nonparam) { prev.p.z0 <- new.p.z0; prev.p.z1 <- new.p.z1 } else { prev.gamma <- new.gamma; prev.eta <- new.eta }
    it <- it + 1L
  } # end EM loop

  if (!converged) {
    if (full.output) {
      return(list(delta = NA, delta.s = NA, R.s = NA, betas = rep(NA,4), sigma = NA, p0 = NA, p1 = NA, alphas = c(NA,NA)))
    } else {
      return(list(delta = NA, delta.s = NA, R.s = NA))
    }
  }

  # compute alpha0, alpha1 and estimands
  beta0 <- prev.beta[1]; beta1 <- prev.beta[2]; beta2 <- prev.beta[3]; beta3 <- prev.beta[4]
  if (nonparam) {
    alpha0 <- if (length(S.z0) > 0 && length(prev.p.z0) > 0) sum(S.z0 * as.numeric(prev.p.z0)) else NA_real_
    alpha1 <- if (length(S.z1) > 0 && length(prev.p.z1) > 0) sum(S.z1 * as.numeric(prev.p.z1)) else NA_real_
  } else {
    alpha0 <- prev.gamma[1]; alpha1 <- alpha0 + prev.gamma[2]
  }

  delta   <- beta1 + (beta2 + beta3) * alpha1 - beta2 * alpha0
  delta.s <- beta1 + beta3 * alpha0
  R.S     <- 1 - delta.s / delta


  if (full.output) {
    return(list(delta = delta, delta.s = delta.s, R.s = R.S,
                betas = prev.beta, sigma = prev.sigma, p0 = prev.p.z0, p1 = prev.p.z1, alphas = c(alpha0, alpha1)))
  } else {
    return(list(delta = delta, delta.s = delta.s, R.s = R.S))
  }
}
