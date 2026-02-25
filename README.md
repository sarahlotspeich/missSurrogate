Surrogate marker evaluation with missing data
================

<p style="display:inline-block;">
  <img src="miss_surrogate_hex.png" width="200" title="a woman with a tiara and pageant sash that says 'Miss Surrogate'">
  <h1>Surrogate marker validation with missing data</h1>
</p>


## Package Installation

Installation of the `missSurrogate` from GitHub requires the
[`devtools`](https://www.r-project.org/nosvn/pandoc/devtools.html)
package and can be done in the following way.

``` r
devtools::install_github(repo = "sarahlotspeich/missSurrogate")
```

Then, you are ready to load the package as follows.

``` r
library(missSurrogate)
```

*Art for the hex sticker was made with the help of Dall-E.*

## Setup

To illustrate the various functions in the `missSurrogate` package, we
simulate the following dataset. It is designed to mimic a clinical trial
with a continuous primary outcome $Y$, continuous surrogate marker for
it $S$, and random treatment assignment $Z$. However, $S$ is missing for
some patients in the trial. Specific details of how the data are
generated can be found in the table below.

| Variable: Description | Simulated Distribution |
|----|----|
| $Z$: Treatment assignment | $50/50$ assignment to control ($Z=0$) vs treatment ($Z=1$) |
| $S$: Surrogate marker | $S\|Z \sim \textrm{Normal}(\mu = 5 + Z, \sigma^2 = 1 + 3Z)$ |
| $Y$: Primary outcome | $Y = 2 + Z + 5S + Z \times S + \epsilon$, where $\epsilon \overset{\textrm{iid}}{\sim} \textrm{Normal}(\mu = 0, \sigma^2 = 1)$ |
| $O$: Observation indicator for $S$ | $O\|Y \sim \textrm{Bernoulli}(\pi = 1 / \left\{1 + \exp(-0.015Y)\right\})$ |

Since the probability that surrogate marker $S$ is non-missing depends
on $Y$, this setting is **missingness at random (MAR)**. For any
simulated patients with $O = 0$, their surrogate marker $S$ will be
redacted/treated as missing.

The following code chunk simulates data for a trial of $n = 2000$
patients as described.

``` r
# Simulate data 
## Treatment assignment
Z = rep(x = c(0, 1), each = 1000)
## Surrogate marker 
S = rnorm(n = 2000, mean = 5 + Z, sd = sqrt(1 + 3 * Z))
## Primary outcome 
eps = rnorm(n = 2000, mean = 0, sd = 1)
Y = 2 + Z + 5 * S + Z * S + eps
## Observation indicator for surrogate marker
O = rbinom(n = 2000, size = 1, prob = 1 / (1 + exp(- 0.015 * Y)))
## Redact surrogate markers if O = 0
S[O == 0] = NA 
## Build dataframe
dat = data.frame(Z, S, Y, O)
head(dat) ### view first few rows 
```

    ##   Z        S        Y O
    ## 1 0 5.498679 29.46668 1
    ## 2 0 4.143396 22.26453 1
    ## 3 0       NA 33.11936 0
    ## 4 0 4.123199 21.84686 1
    ## 5 0       NA 28.18187 0
    ## 6 0 4.594432 26.41873 1

## Functionality

The `missSurrogate` package contains a single function, `R.s.miss()`,
that can estimate the proportion of treatment effect explained (PTE) of
a surrogate marker subject to missing data. The arguments to
`R.s.miss()` are as follows:

- `sone`: vector of surrogate markers for the treatment group
- `szero`: vector of surrogate markers for the control group
- `yone`: vector of primary outcomes for the treatment group
- `yzero`: vector of primary outcomes for the control group
- `wone`: vector of inverse probability weights for the treatment group
  (only applicable for certain missing data approaches)
- `wzero`: vector of inverse probability weights for the control group
  (only applicable for certain missing data approaches)
- `type`: character specifying whether `"robust"` (nonparametric) or
  `"model"` (parametric) PTE estimation is desired
- `max.it`: maximum number of iterations attempted before quitting
  without converging (default is `10000`)
- `tol`: numeric tolerance used to define convergence (default is
  `0.001`)
- `ipw.formula`: model formula for the inverse probability weights (only
  applicable for certain missing data approaches when confidence
  intervals are requested)
- `conf.int`: whether bootstrapped confidence intervals are needed
  (default is `FALSE`)

The user can choose between **two estimation approaches**, nonparametric
and parametric, and **three missing data corrections**, complete case,
inverse probability weighting (IPW), and semiparametric maximum
likelihood estimation (SMLE).

In the subsections that follow, we illustrate how this function can be
used to apply these different estimation methods and missing data
corrections using the simulated `dat`. We define the following separate
vectors for the surrogate markers and primary outcomes in the two
treatment groups, in alignment with the arguments. We also define
separate vectors of the observed indicators.

``` r
# Define separate vectors of Y,S,O for treatment/control
s1 = dat$S[dat$Z == 1]
s0 = dat$S[dat$Z == 0]
y1 = dat$Y[dat$Z == 1]
y0 = dat$Y[dat$Z == 0]
o1 = dat$O[dat$Z == 1]
o0 = dat$O[dat$Z == 0]
```

## Complete Case Estimation

To obtain the complete case estimates, we subset the surrogate marker
and primary outcome vectors to the rows corresponding to patients with
non-missing surrogates (i.e., with $O = 1$). We begin with the complete
case **parametric PTE estimates**.

``` r
# Complete case analysis 
## Parametric PTE 
param_cc = R.s.miss(sone = s1[o1 == 1], 
                    szero = s0[o0 == 1], 
                    yone = y1[o1 == 1], 
                    yzero = y0[o0 == 1], 
                    type = "model", 
                    conf.int = TRUE)
param_cc
```

    ## $delta
    ## [1] 12.19665
    ## 
    ## $delta.s
    ## [1] 6.043017
    ## 
    ## $R.s
    ## [1] 0.5045345
    ## 
    ## $delta.var
    ## [1] 0.2505104
    ## 
    ## $delta.s.var
    ## [1] 0.005418184
    ## 
    ## $R.s.var
    ## [1] 0.0004841111
    ## 
    ## $conf.int.normal.delta
    ## [1] 11.21565 13.17765
    ## 
    ## $conf.int.quantile.delta
    ## [1] 11.18660 13.15946
    ## 
    ## $conf.int.normal.delta.s
    ## [1] 5.898745 6.187289
    ## 
    ## $conf.int.quantile.delta.s
    ## [1] 5.905278 6.188543
    ## 
    ## $conf.int.normal.R.s
    ## [1] 0.4614095 0.5476594
    ## 
    ## $conf.int.quantile.R.s
    ## [1] 0.4570221 0.5418250

The returned list, saved as `param_cc`, contains the following:

- `delta`: the estimated overall average treatment effect on $Y$
- `delta.s`: the estimated residual treatment effect
- `R.s`: the estimated proportion of treatment effect explained by $S$
- `delta.var`, `delta.s.var`, and `R.s.var`: bootstrap-estimated
  variances for the corresponding quantities
- `conf.int.normal.delta`, `conf.int.normal.delta.s`, and
  `conf.int.normal.R.s`: Wald-type 95% confidence interval for the
  corresponding quantities based on the Normal distribution and
  bootstrap-estimated variances
- `conf.int.quantile.delta`, `conf.int.quantile.delta.s`, and
  `conf.int.quantile.R.s`: quantile-based 95% confidence interval for
  the corresponding quantities based on the central 95% of the
  bootstrapped values

To instead obtain the complete case **nonparametric PTE estimates**, we
only need to switch from `type = "model"` to `type = "robust"`. The
output list follows the same structure as above.

``` r
## Nonparametric PTE 
nonparam_cc = R.s.miss(sone = s1[o1 == 1], 
                       szero = s0[o0 == 1], 
                       yone = y1[o1 == 1], 
                       yzero = y0[o0 == 1], 
                       type = "robust", 
                       conf.int = TRUE)
```

## Inverse Probability Weighting

Before using `R.s.miss()` for IPW, we need to calculate the appropriate
weights, defined as the inverse of the probability of a patient having a
nonmissing surrogate marker (potentially conditional upon other
patient-level information). For this simulation, we know that whether or
not a patient’s surrogate $S$ is observed (i.e., $O = 1$) depends on
their primary outcome $Y$. Thus, we fit a logistic regression model for
$O \sim Y$ to estimate these weights.

``` r
## Calculate weights for IPW approaches
o = c(o1, o0)
y = c(y1, y0)
ipw_fit = glm(formula = o ~ y, 
              family = "binomial")
p1 = 1 / (1 + exp(-(ipw_fit$coefficients[1] +
                      ipw_fit$coefficients[2] * y1)))
w1 = 1 / p1 ### vector of weights for treatment group
p0 = 1 / (1 + exp(-(ipw_fit$coefficients[1] + ipw_fit$coefficients[2] * y0)))
w0 = 1 / p0 ### vector of weights for control group
```

Once we have the weight vectors `w1` and `w0` for the treatment and
control groups, we can proceed with obtain the IPW-corrected parametric
or nonparametric PTE estimates as follows. In addition, we need to
provide the `ipw_formula` for our weighting model.

``` r
# IPW analysis 
## Parametric PTE 
param_ipw = R.s.miss(sone = s1, 
                     szero = s0, 
                     yone = y1, 
                     yzero = y0, 
                     wone = w1, 
                     wzero = w0,
                     ipw.formula = o ~ y,
                     type = "model", 
                     conf.int = TRUE)
## Nonparametric PTE 
nonparam_ipw = R.s.miss(sone = s1, 
                        szero = s0, 
                        yone = y1, 
                        yzero = y0, 
                        wone = w1, 
                        wzero = w0,
                        ipw.formula = o ~ y,
                        type = "robust", 
                        conf.int = TRUE)
```

Based on our data generation setup, the following weighting models are
also *possible* (but only the first one is correctly specified):

1.  If we thought that the surrogate was missing completely at random
    (i.e., independently of all other variables), then we would use
    `ipw.formula = o ~ 1`, since our weights model only needs an
    intercept.
2.  If we thought that the surrogate was MAR given treatment group $Z$,
    rather than $Y$, it would be `ipw.formula = o ~ z`.
3.  If we thought that the surrogate was MAR given primary outcome $Y$
    and treatment group $Z$, it would be `ipw.formula = o ~ y + z`.
4.  If we thought that the surrogate was MAR given primary outcome $Y$
    and treatment group $Z$ *and* their interaction, it would be
    `ipw.formula = o ~ y * z`.

Notice that in all cases, `ipw.formula` is in terms of general `s`, `y`,
or `z`, rather than the subgroups (e.g., `szero` and `sone`). Currently,
the `missSurrogate` package cannot accommodate additional variables in
this model.

## Semiparametric Maximum Likelihood Estimation

Unlike complete case estimation and IPW, the SMLE does not require any
subsetting or data pre-processing. Currently, the SMLE missing
correction is only available for the parametric PTE estimation approach,
and it can be applied as follows.

``` r
# SMLE analysis 
## Parametric PTE 
Rparam_miss_smle = R.s.miss(sone = s1, 
                            szero = s0,
                            yone = y1,
                            yzero = y0, 
                            type = "model", 
                            conf.int = TRUE) 
```
