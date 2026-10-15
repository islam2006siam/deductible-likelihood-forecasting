# Packages and settings ----
suppressPackageStartupMessages({
  library(data.table); library(rstan); library(loo); library(MASS)
  library(scoringRules); library(ggplot2); library(patchwork); library(scales)
})
options(mc.cores = 4)
rstan_options(auto_write = TRUE)

SEED <- 20260703
data_dir <- "data"
out_dir <- "outputs"
fit_dir <- file.path(out_dir, "fits")
fig_dir <- file.path(out_dir, "figures")
for (d in c(out_dir, fit_dir, fig_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# Stan programs ----
sev_functions <- "
functions {
  real gb2_lpdf(real y, real a, real p, real q, real b) {
    real z = a * (log(y) - log(b));
    return log(a) + p * z - log(y) - lbeta(p, q) - (p + q) * log1p_exp(z);
  }
  real gb2_cdf_val(real y, real a, real p, real q, real b) {
    real z = a * (log(y) - log(b));
    return inc_beta(p, q, inv_logit(z));
  }
  real burr_lpdf(real y, real a, real q, real b) {
    real z = a * (log(y) - log(b));
    return log(a) + log(q) + z - log(y) - (q + 1) * log1p_exp(z);
  }
  real burr_lccdf_val(real y, real a, real q, real b) {
    real z = a * (log(y) - log(b));
    return -q * log1p_exp(z);
  }
  real sev_ld(real y, int family, real loc, real sigma, real alpha,
              real a, real p, real q) {
    if (family == 1) return lognormal_lpdf(y | loc, sigma);
    else if (family == 2) return gamma_lpdf(y | alpha, alpha / exp(loc));
    else if (family == 4) return burr_lpdf(y | a, q, exp(loc));
    else return gb2_lpdf(y | a, p, q, exp(loc));
  }
  real sev_lcc(real y, int family, real loc, real sigma, real alpha,
               real a, real p, real q) {
    if (family == 1) return lognormal_lccdf(y | loc, sigma);
    else if (family == 2) return gamma_lccdf(y | alpha, alpha / exp(loc));
    else if (family == 4) return burr_lccdf_val(y, a, q, exp(loc));
    else return log1m(gb2_cdf_val(y, a, p, q, exp(loc)));
  }
}
"

joint_code <- paste0(sev_functions, "
data {
  int<lower=1> n_py;
  array[n_py] int<lower=0> N;
  int<lower=1> K_f;
  matrix[n_py, K_f] Xf;
  int<lower=1> K_s;
  matrix[n_py, K_s] Xs_py;
  vector<lower=0>[n_py] d_py;
  int<lower=1> n_pid;
  array[n_py] int<lower=1, upper=n_pid> pid;
  int<lower=0> n_cl;
  vector<lower=0>[n_cl] y;
  array[n_cl] int<lower=0, upper=1> cens;
  vector<lower=0>[n_cl] d_cl;
  vector<lower=0>[n_cl] u_cl;
  matrix[n_cl, K_s] Xs_cl;
  array[n_cl] int<lower=1, upper=n_py> py_of_cl;
  int<lower=1, upper=4> family;
  int<lower=0, upper=1> ltrc;
  int<lower=0, upper=2> dep;
  int<lower=0, upper=1> nb;
}
parameters {
  vector[K_f] beta_f;
  vector[K_s] beta_s;
  real<lower=0> sigma_b;
  vector[n_pid] b_raw;
  real<lower=0> sigma_ln;
  real<lower=0> alpha_g;
  real<lower=0> gb2_a;
  real<lower=0> gb2_p;
  real<lower=0> gb2_q;
  real<lower=0> phi_nb;
  real gamma_dep;
}
transformed parameters {
  vector[n_pid] b = sigma_b * b_raw;
  vector[n_py] eta_f;
  vector[n_py] log_thin;
  vector[n_py] mu_obs;
  vector[n_py] resid_N;
  {
    vector[n_py] loc_py = Xs_py * beta_s;
    for (i in 1:n_py) {
      real loc_i = loc_py[i];
      if (dep == 2) loc_i += gamma_dep * b[pid[i]];
      eta_f[i] = dot_product(Xf[i], beta_f) + b[pid[i]];
      if (ltrc == 1)
        log_thin[i] = sev_lcc(d_py[i], family, loc_i, sigma_ln, alpha_g,
                              gb2_a, gb2_p, gb2_q);
      else
        log_thin[i] = 0;
      mu_obs[i] = exp(eta_f[i] + log_thin[i]);
      resid_N[i] = (N[i] - mu_obs[i]) / sqrt(mu_obs[i] + 1e-8);
    }
  }
}
model {
  beta_f ~ normal(0, 2.5);
  beta_s ~ normal(0, 2.5);
  sigma_b ~ exponential(1);
  b_raw ~ std_normal();
  sigma_ln ~ normal(0, 2);
  alpha_g ~ lognormal(0, 1);
  gb2_a ~ lognormal(0, 1);
  gb2_p ~ lognormal(0, 1);
  gb2_q ~ lognormal(0, 1);
  phi_nb ~ lognormal(1, 1);
  gamma_dep ~ normal(0, 1);
  if (nb == 0) {
    for (i in 1:n_py) target += poisson_log_lpmf(N[i] | eta_f[i] + log_thin[i]);
  } else {
    for (i in 1:n_py) target += neg_binomial_2_lpmf(N[i] | mu_obs[i] + 1e-10, phi_nb);
  }
  for (k in 1:n_cl) {
    int i = py_of_cl[k];
    real loc_k = dot_product(Xs_cl[k], beta_s);
    if (dep == 1) loc_k += gamma_dep * resid_N[i];
    if (dep == 2) loc_k += gamma_dep * b[pid[i]];
    if (cens[k] == 1 && ltrc == 1) {
      target += sev_lcc(u_cl[k], family, loc_k, sigma_ln, alpha_g, gb2_a, gb2_p, gb2_q)
              - sev_lcc(d_cl[k], family, loc_k, sigma_ln, alpha_g, gb2_a, gb2_p, gb2_q);
    } else if (ltrc == 1) {
      target += sev_ld(y[k], family, loc_k, sigma_ln, alpha_g, gb2_a, gb2_p, gb2_q)
              - sev_lcc(d_cl[k], family, loc_k, sigma_ln, alpha_g, gb2_a, gb2_p, gb2_q);
    } else {
      target += sev_ld(y[k], family, loc_k, sigma_ln, alpha_g, gb2_a, gb2_p, gb2_q);
    }
  }
}
generated quantities {
  vector[n_py] log_lik;
  for (i in 1:n_py) {
    if (nb == 0)
      log_lik[i] = poisson_log_lpmf(N[i] | eta_f[i] + log_thin[i]);
    else
      log_lik[i] = neg_binomial_2_lpmf(N[i] | mu_obs[i] + 1e-10, phi_nb);
  }
  for (k in 1:n_cl) {
    int i = py_of_cl[k];
    real loc_k = dot_product(Xs_cl[k], beta_s);
    if (dep == 1) loc_k += gamma_dep * resid_N[i];
    if (dep == 2) loc_k += gamma_dep * b[pid[i]];
    if (cens[k] == 1 && ltrc == 1) {
      log_lik[i] += sev_lcc(u_cl[k], family, loc_k, sigma_ln, alpha_g, gb2_a, gb2_p, gb2_q)
                  - sev_lcc(d_cl[k], family, loc_k, sigma_ln, alpha_g, gb2_a, gb2_p, gb2_q);
    } else if (ltrc == 1) {
      log_lik[i] += sev_ld(y[k], family, loc_k, sigma_ln, alpha_g, gb2_a, gb2_p, gb2_q)
                  - sev_lcc(d_cl[k], family, loc_k, sigma_ln, alpha_g, gb2_a, gb2_p, gb2_q);
    } else {
      log_lik[i] += sev_ld(y[k], family, loc_k, sigma_ln, alpha_g, gb2_a, gb2_p, gb2_q);
    }
  }
}
")

sev_code <- paste0(sev_functions, "
data {
  int<lower=1> n_cl;
  vector<lower=0>[n_cl] y;
  array[n_cl] int<lower=0, upper=1> cens;
  vector<lower=0>[n_cl] d_cl;
  vector<lower=0>[n_cl] u_cl;
  int<lower=1> K_s;
  matrix[n_cl, K_s] Xs;
  int<lower=1, upper=4> family;
  int<lower=0, upper=1> ltrc;
}
parameters {
  vector[K_s] beta_s;
  real<lower=0> sigma_ln;
  real<lower=0> alpha_g;
  real<lower=0> gb2_a;
  real<lower=0> gb2_p;
  real<lower=0> gb2_q;
}
model {
  beta_s ~ normal(0, 2.5);
  sigma_ln ~ normal(0, 2);
  alpha_g ~ lognormal(0, 1);
  gb2_a ~ lognormal(0, 1);
  gb2_p ~ lognormal(0, 1);
  gb2_q ~ lognormal(0, 1);
  for (k in 1:n_cl) {
    real loc_k = dot_product(Xs[k], beta_s);
    if (cens[k] == 1 && ltrc == 1)
      target += sev_lcc(u_cl[k], family, loc_k, sigma_ln, alpha_g, gb2_a, gb2_p, gb2_q)
              - sev_lcc(d_cl[k], family, loc_k, sigma_ln, alpha_g, gb2_a, gb2_p, gb2_q);
    else if (ltrc == 1)
      target += sev_ld(y[k], family, loc_k, sigma_ln, alpha_g, gb2_a, gb2_p, gb2_q)
              - sev_lcc(d_cl[k], family, loc_k, sigma_ln, alpha_g, gb2_a, gb2_p, gb2_q);
    else
      target += sev_ld(y[k], family, loc_k, sigma_ln, alpha_g, gb2_a, gb2_p, gb2_q);
  }
}
generated quantities {
  vector[n_cl] log_lik;
  for (k in 1:n_cl) {
    real loc_k = dot_product(Xs[k], beta_s);
    if (cens[k] == 1 && ltrc == 1)
      log_lik[k] = sev_lcc(u_cl[k], family, loc_k, sigma_ln, alpha_g, gb2_a, gb2_p, gb2_q)
                 - sev_lcc(d_cl[k], family, loc_k, sigma_ln, alpha_g, gb2_a, gb2_p, gb2_q);
    else if (ltrc == 1)
      log_lik[k] = sev_ld(y[k], family, loc_k, sigma_ln, alpha_g, gb2_a, gb2_p, gb2_q)
                 - sev_lcc(d_cl[k], family, loc_k, sigma_ln, alpha_g, gb2_a, gb2_p, gb2_q);
    else
      log_lik[k] = sev_ld(y[k], family, loc_k, sigma_ln, alpha_g, gb2_a, gb2_p, gb2_q);
  }
}
")

sim_dep2_code <- "
data {
  int<lower=1> n;
  array[n] int<lower=0> N;
  vector[n] x1;
  vector[n] x2;
  vector<lower=0>[n] d;
  vector<lower=0>[n] u;
  int<lower=0> n_cl;
  vector<lower=0>[n_cl] y;
  array[n_cl] int<lower=0, upper=1> cens;
  array[n_cl] int<lower=1, upper=n> py_of_cl;
  int<lower=0, upper=1> ltrc;
  int<lower=0, upper=1> depon;
}
parameters {
  vector[3] beta_f;
  vector[3] beta_s;
  real<lower=0> delta;
  real gam;
  real<lower=0> sigma;
  vector[n] theta_raw;
}
model {
  beta_f ~ normal(0, 2.5);
  beta_s ~ normal(0, 2.5);
  delta ~ normal(0, 1);
  gam ~ normal(0, 1);
  sigma ~ normal(0, 2);
  theta_raw ~ std_normal();
  {
    vector[n] eta = beta_f[1] + beta_f[2] * x1 + beta_f[3] * x2 + delta * theta_raw;
    vector[n] loc = beta_s[1] + beta_s[2] * x1 + beta_s[3] * x2
                    + (depon == 1 ? gam : 0) * theta_raw;
    vector[n] lthin;
    for (i in 1:n) {
      lthin[i] = ltrc == 1 ? lognormal_lccdf(d[i] | loc[i], sigma) : 0;
      target += poisson_log_lpmf(N[i] | eta[i] + lthin[i]);
    }
    for (k in 1:n_cl) {
      int i = py_of_cl[k];
      if (cens[k] == 1 && ltrc == 1)
        target += lognormal_lccdf(u[i] | loc[i], sigma) - lognormal_lccdf(d[i] | loc[i], sigma);
      else if (ltrc == 1)
        target += lognormal_lpdf(y[k] | loc[i], sigma) - lognormal_lccdf(d[i] | loc[i], sigma);
      else
        target += lognormal_lpdf(y[k] | loc[i], sigma);
    }
  }
}
"

# Functions ----
lgpif_prepare <- function(dat) {
  tp <- copy(dat$train_pol); hp <- copy(dat$hold_pol)
  tc <- copy(dat$train_clm); hc <- copy(dat$hold_clm)
  tc <- tc[!is.na(BCcov)]; hc <- hc[!is.na(BCcov)]
  tp[, py := .I]; hp[, py := .I]
  tc <- merge(tc, tp[, .(PolicyNum, Year, py)], by = c("PolicyNum", "Year"))
  hc <- merge(hc, hp[, .(PolicyNum, Year, py)], by = c("PolicyNum", "Year"))
  Xf_cols <- function(d) cbind(1, d$zLnCov, d$TypeCity, d$TypeCounty, d$TypeSchool,
                               d$TypeTown, d$TypeVillage, d$NoClaimCredit, d$Fire5)
  Xs_cols <- function(d) cbind(1, d$zLnCov, d$TypeCity, d$TypeCounty, d$TypeSchool,
                               d$TypeTown, d$TypeVillage)
  setorder(tc, py); setorder(hc, py)
  list(tp = tp, hp = hp, tc = tc, hc = hc,
       Xf_tr = Xf_cols(tp), Xs_tr_py = Xs_cols(tp), Xs_tr_cl = Xs_cols(tc),
       Xf_ho = Xf_cols(hp), Xs_ho_py = Xs_cols(hp), Xs_ho_cl = Xs_cols(hc),
       n_pid = max(dat$pol$pid), std = dat$std)
}

lgpif_stan_data <- function(pp, family, ltrc, dep, nb) {
  tp <- pp$tp; tc <- pp$tc
  list(n_py = nrow(tp), N = tp$FreqA, K_f = ncol(pp$Xf_tr), Xf = pp$Xf_tr,
       K_s = ncol(pp$Xs_tr_py), Xs_py = pp$Xs_tr_py, d_py = tp$Deduct / 1000,
       n_pid = pp$n_pid, pid = tp$pid,
       n_cl = nrow(tc), y = pmin(tc$Claim, tc$BCcov) / 1000, cens = tc$censored,
       d_cl = tc$Deduct / 1000, u_cl = tc$BCcov / 1000, Xs_cl = pp$Xs_tr_cl,
       py_of_cl = tc$py, family = family, ltrc = ltrc, dep = dep, nb = nb)
}

lgpif_sev_data <- function(pp, family, ltrc) {
  tc <- pp$tc
  list(n_cl = nrow(tc), y = pmin(tc$Claim, tc$BCcov) / 1000, cens = tc$censored,
       d_cl = tc$Deduct / 1000, u_cl = tc$BCcov / 1000,
       K_s = ncol(pp$Xs_tr_cl), Xs = pp$Xs_tr_cl, family = family, ltrc = ltrc)
}

joint_inits <- function(K_f, K_s) function() list(
  beta_f = c(-0.5, rep(0, K_f - 1)), beta_s = c(1.0, rep(0, K_s - 1)),
  sigma_b = 0.5, gb2_a = 1.5, gb2_p = 1, gb2_q = 1.5,
  sigma_ln = 1.5, alpha_g = 0.8, phi_nb = 3, gamma_dep = 0)

sev_logS_v <- function(x, loc, fam, sp) {
  if (fam == 1) plnorm(x, loc, sp$sigma, lower.tail = FALSE, log.p = TRUE)
  else if (fam == 2) pgamma(x, shape = sp$alpha, rate = sp$alpha / exp(loc),
                            lower.tail = FALSE, log.p = TRUE)
  else if (fam == 4) -sp$q * log1p(exp(pmin(sp$a * (log(x) - loc), 700)))
  else pbeta(plogis(sp$a * (log(x) - loc)), sp$p, sp$q,
             lower.tail = FALSE, log.p = TRUE)
}
sev_logf_v <- function(x, loc, fam, sp) {
  if (fam == 1) dlnorm(x, loc, sp$sigma, log = TRUE)
  else if (fam == 2) dgamma(x, shape = sp$alpha, rate = sp$alpha / exp(loc), log = TRUE)
  else if (fam == 4) {
    z <- sp$a * (log(x) - loc)
    log(sp$a) + log(sp$q) + z - log(x) - (sp$q + 1) * log1p(exp(pmin(z, 700)))
  } else {
    z <- sp$a * (log(x) - loc)
    log(sp$a) + sp$p * z - log(x) - lbeta(sp$p, sp$q) - (sp$p + sp$q) * log1p(exp(pmin(z, 700)))
  }
}
sev_cdf_v <- function(x, loc, fam, sp) {
  if (fam == 1) plnorm(x, loc, sp$sigma)
  else if (fam == 2) pgamma(x, shape = sp$alpha, rate = sp$alpha / exp(loc))
  else if (fam == 4) 1 - exp(-sp$q * log1p(exp(pmin(sp$a * (log(x) - loc), 700))))
  else pbeta(plogis(sp$a * (log(x) - loc)), sp$p, sp$q)
}
sev_q_v <- function(u, loc, fam, sp) {
  if (fam == 1) qlnorm(u, loc, sp$sigma)
  else if (fam == 2) qgamma(u, shape = sp$alpha, rate = sp$alpha / exp(loc))
  else if (fam == 4) exp(loc) * ((1 - u)^(-1 / sp$q) - 1)^(1 / sp$a)
  else { qb <- qbeta(u, sp$p, sp$q); exp(loc) * (qb / (1 - qb))^(1 / sp$a) }
}

get_sp <- function(fit, s) {
  f <- fit$family
  if (f == 1) list(sigma = fit$draws$sigma_ln[s])
  else if (f == 2) list(alpha = fit$draws$alpha_g[s])
  else if (f == 4) list(a = fit$draws$gb2_a[s], q = fit$draws$gb2_q[s])
  else list(a = fit$draws$gb2_a[s], p = fit$draws$gb2_p[s], q = fit$draws$gb2_q[s])
}

predict_model <- function(fit, Xf, Xs_py, d_py, u_py, pid_v, N_obs, claims,
                          seen_pids, S_use = 1000, seed = 1) {
  set.seed(seed)
  dr <- fit$draws
  S_tot <- length(dr$sigma_b)
  ss <- round(seq(1, S_tot, length.out = S_use))
  n <- nrow(Xf); fam <- fit$family
  dep <- fit$spec$dep; ltrc <- fit$spec$ltrc; nb <- fit$nb
  Ymat <- matrix(0, n, S_use)
  lpN <- matrix(NA_real_, n, S_use)
  n_cl <- length(claims$y)
  lf_cl <- if (n_cl > 0) matrix(NA_real_, n_cl, S_use) else NULL
  for (js in seq_along(ss)) {
    s <- ss[js]
    sp <- get_sp(fit, s)
    b_v <- dr$b[s, ][pid_v]
    new_pol <- !(pid_v %in% seen_pids)
    if (any(new_pol)) b_v[new_pol] <- rnorm(sum(new_pol), 0, dr$sigma_b[s])
    eta_f <- as.vector(Xf %*% dr$beta_f[s, ]) + b_v
    loc0 <- as.vector(Xs_py %*% dr$beta_s[s, ])
    loc_th <- loc0 + if (dep == 2) dr$gamma_dep[s] * b_v else 0
    lthin <- if (ltrc == 1) sev_logS_v(d_py, loc_th, fam, sp) else 0
    mu <- exp(eta_f + lthin)
    lpN[, js] <- if (nb == 0) dpois(N_obs, mu, log = TRUE)
                 else dnbinom(N_obs, mu = mu, size = dr$phi_nb[s], log = TRUE)
    Nsim <- if (nb == 0) rpois(n, mu) else rnbinom(n, mu = mu, size = dr$phi_nb[s])
    tot <- rep(0, n)
    idx <- which(Nsim > 0)
    if (length(idx) > 0) {
      py <- rep(idx, Nsim[idx])
      loc_k <- loc_th[py]
      if (dep == 1) {
        resid <- (Nsim - mu) / sqrt(mu)
        loc_k <- loc_k + dr$gamma_dep[s] * resid[py]
      }
      if (ltrc == 1) {
        Fd <- sev_cdf_v(d_py[py], loc_k, fam, sp)
        X <- sev_q_v(Fd + runif(length(py)) * (1 - Fd), loc_k, fam, sp)
      } else {
        X <- sev_q_v(runif(length(py)), loc_k, fam, sp)
      }
      pay <- pmin(X, u_py[py])
      agg <- rowsum(pay, py)
      tot[as.integer(rownames(agg))] <- agg
    }
    Ymat[, js] <- tot
    if (n_cl > 0) {
      loc_c <- as.vector(claims$Xs %*% dr$beta_s[s, ])
      i_py <- claims$py
      if (dep == 2) loc_c <- loc_c + dr$gamma_dep[s] * b_v[i_py]
      if (dep == 1) {
        resid_o <- (N_obs - mu) / sqrt(mu)
        loc_c <- loc_c + dr$gamma_dep[s] * resid_o[i_py]
      }
      if (ltrc == 1) {
        lS_d <- sev_logS_v(claims$d, loc_c, fam, sp)
        lf <- ifelse(claims$cens == 1,
                     sev_logS_v(claims$u, loc_c, fam, sp),
                     sev_logf_v(claims$y, loc_c, fam, sp))
        lf_cl[, js] <- lf - lS_d
      } else {
        lf_cl[, js] <- sev_logf_v(claims$y, loc_c, fam, sp)
      }
    }
  }
  list(Ymat = Ymat, lpN = lpN, lf_cl = lf_cl)
}

log_mean_exp_rows <- function(M) {
  mx <- apply(M, 1, max)
  mx + log(rowMeans(exp(M - mx)))
}

pit_from_samples <- function(y, Ymat, seed = 3) {
  set.seed(seed)
  S <- ncol(Ymat)
  lt <- rowSums(Ymat < y)
  eq <- rowSums(Ymat == y)
  (lt + runif(length(y)) * (eq + 1)) / (S + 1)
}

gini_lift <- function(pred, y) {
  o <- order(pred)
  L <- cumsum(y[o]) / sum(y)
  Lprev <- c(0, L[-length(L)])
  auc <- sum((L + Lprev) / 2) / length(y)
  1 - 2 * auc
}

cqr_correct <- function(y_cal, Ymat_cal, alphas = c(0.2, 0.1)) {
  out <- list()
  for (a in alphas) {
    ql <- apply(Ymat_cal, 1, quantile, probs = a / 2, names = FALSE)
    qh <- apply(Ymat_cal, 1, quantile, probs = 1 - a / 2, names = FALSE)
    E <- pmax(ql - y_cal, y_cal - qh)
    n <- length(E)
    k <- ceiling((1 - a) * (n + 1))
    out[[paste0("a", a * 100)]] <- sort(E)[min(k, n)]
  }
  out
}

prem_curve <- function(fit, xs, xf, d_grid, u_grid, u_for_d, d_for_u,
                       S_use = 500, M = 4000, seed = 11) {
  set.seed(seed)
  dr <- fit$draws; fam <- fit$family; dep <- fit$spec$dep
  S_tot <- length(dr$sigma_b)
  ss <- round(seq(1, S_tot, length.out = S_use))
  PPd <- matrix(0, length(d_grid), S_use)
  PPu <- matrix(0, length(u_grid), S_use)
  for (js in seq_along(ss)) {
    s <- ss[js]
    sp <- get_sp(fit, s)
    b <- rnorm(M, 0, dr$sigma_b[s])
    lam <- exp(sum(xf * dr$beta_f[s, ]) + b)
    loc <- rep(sum(xs * dr$beta_s[s, ]), M) + if (dep == 2) dr$gamma_dep[s] * b else 0
    Ns <- rpois(M, lam)
    idx <- which(Ns > 0)
    if (!length(idx)) next
    py <- rep(idx, Ns[idx])
    loc_k <- loc[py]
    if (dep == 1)
      loc_k <- loc_k + dr$gamma_dep[s] * (Ns[py] - lam[py]) / sqrt(lam[py])
    X <- sev_q_v(runif(length(py)), loc_k, fam, sp)
    for (gi in seq_along(d_grid))
      PPd[gi, js] <- sum(pmax(pmin(X, u_for_d) - d_grid[gi], 0)) / M
    for (gi in seq_along(u_grid))
      PPu[gi, js] <- sum(pmax(pmin(X, u_grid[gi]) - d_for_u, 0)) / M
  }
  list(PPd = PPd, PPu = PPu)
}

dgb2 <- function(y, a, p, q, b, log = FALSE) {
  z <- a * (log(y) - log(b))
  lf <- log(a) + p * z - log(y) - lbeta(p, q) - (p + q) * log1p(exp(pmin(z, 700)))
  if (log) lf else exp(lf)
}
pgb2 <- function(y, a, p, q, b) pbeta(plogis(a * (log(y) - log(b))), p, q)
qgb2 <- function(u, a, p, q, b) b * (qbeta(u, p, q) / (1 - qbeta(u, p, q)))^(1 / a)
rgb2 <- function(n, a, p, q, b) qgb2(runif(n), a, p, q, b)

gauss_hermite_prob <- function(k) {
  i <- seq_len(k - 1)
  J <- matrix(0, k, k)
  J[cbind(i, i + 1)] <- J[cbind(i + 1, i)] <- sqrt(i / 2)
  e <- eigen(J, symmetric = TRUE)
  idx <- order(e$values)
  list(nodes = sqrt(2) * e$values[idx], weights = (e$vectors[1, idx])^2)
}

logsumexp_rows <- function(M) {
  mx <- apply(M, 1, max)
  mx + log(rowSums(exp(M - mx)))
}

`%||%` <- function(a, b) if (is.null(a)) b else a

TRUE_PARS <- list(
  beta_f = c(-0.3, 0.3, 0.2),
  beta_s = c(1.0, 0.4, -0.3),
  sigma  = 1.5,
  gb2    = c(a = 1.4, p = 1.0, q = 1.5),
  delta  = 0.4
)
N_POL <- 3000; N_YR <- 5
GH <- gauss_hermite_prob(7)

.SQ_CACHE <- new.env(parent = emptyenv())
sev_quantiles <- function(family, probs) {
  key <- paste0(family, "|", paste(probs, collapse = "_"))
  if (!is.null(.SQ_CACHE[[key]])) return(.SQ_CACHE[[key]])
  had_seed <- exists(".Random.seed", envir = globalenv())
  if (had_seed) old <- get(".Random.seed", envir = globalenv())
  set.seed(991)
  n <- 4e5
  x1 <- rnorm(n); x2 <- rbinom(n, 1, 0.35)
  loc <- TRUE_PARS$beta_s[1] + TRUE_PARS$beta_s[2] * x1 + TRUE_PARS$beta_s[3] * x2
  X <- if (family == "ln") rlnorm(n, loc, TRUE_PARS$sigma)
       else rgb2(n, TRUE_PARS$gb2["a"], TRUE_PARS$gb2["p"], TRUE_PARS$gb2["q"], exp(loc))
  out <- quantile(X, probs)
  if (had_seed) assign(".Random.seed", old, envir = globalenv())
  .SQ_CACHE[[key]] <- out
  out
}

gen_data <- function(scen, seed) {
  set.seed(seed)
  n <- N_POL * N_YR
  x1 <- rnorm(n); x2 <- rbinom(n, 1, 0.35)
  eta_f <- TRUE_PARS$beta_f[1] + TRUE_PARS$beta_f[2] * x1 + TRUE_PARS$beta_f[3] * x2
  loc0  <- TRUE_PARS$beta_s[1] + TRUE_PARS$beta_s[2] * x1 + TRUE_PARS$beta_s[3] * x2
  if (scen$mech == "dep2") {
    th <- rnorm(n)
    lam <- exp(eta_f + TRUE_PARS$delta * th)
    Nstar <- rpois(n, lam)
    loc_i <- loc0 + scen$gam * th
  } else {
    lam <- exp(eta_f)
    Nstar <- rpois(n, lam)
    loc_i <- loc0 + scen$gam * (Nstar - lam) / sqrt(lam)
  }
  qs <- sev_quantiles(scen$family, c(.30, .35, .45, .55, .65, .80, .90, .999))
  if (scen$dist == "lgpif") {
    dlev <- qs[c(".35" = 2, ".45" = 3, ".55" = 4, ".65" = 5, ".80" = 6)]
    d <- sample(dlev, n, replace = TRUE, prob = c(.47, .19, .12, .17, .05))
    u <- rep(qs[8], n)
  } else {
    d <- rep(qs[1], n); u <- rep(qs[7], n)
  }
  py <- rep(seq_len(n), Nstar)
  Xg <- if (scen$family == "ln") rlnorm(length(py), loc_i[py], TRUE_PARS$sigma)
        else rgb2(length(py), TRUE_PARS$gb2["a"], TRUE_PARS$gb2["p"], TRUE_PARS$gb2["q"],
                  exp(loc_i[py]))
  keep <- Xg > d[py]
  py_o <- py[keep]
  y_o  <- pmin(Xg[keep], u[py_o])
  cens <- as.integer(Xg[keep] >= u[py_o])
  Nobs <- tabulate(py_o, nbins = n)
  list(n = n, x1 = x1, x2 = x2, d = d, u = u, Nobs = Nobs,
       py_o = py_o, y = y_o, cens = cens, scen = scen)
}

sev_logf <- function(y, loc, sp, family) {
  if (family == "ln") dlnorm(y, loc, sp$sigma, log = TRUE)
  else dgb2(y, sp$a, sp$p, sp$q, exp(loc), log = TRUE)
}
sev_logS <- function(y, loc, sp, family) {
  if (family == "ln") pnorm((log(y) - loc) / sp$sigma, lower.tail = FALSE, log.p = TRUE)
  else log(pmax(1 - pgb2(y, sp$a, sp$p, sp$q, exp(loc)), 1e-300))
}

unpack <- function(par, model, mech, family) {
  out <- list(beta_f = par[1:3], beta_s = par[4:6])
  i <- 7
  out$gam <- 0; out$delta <- 0
  if (model %in% c("MC", "MD")) { out$gam <- par[i]; i <- i + 1 }
  if (mech == "dep2") { out$delta <- exp(par[i]); i <- i + 1 }
  if (family == "ln") { out$sp <- list(sigma = exp(par[i])) }
  else out$sp <- list(a = exp(par[i]), p = 1, q = exp(par[i + 1]))
  out
}

neg_log_prior <- function(p, family) {
  lp <- sum(dnorm(c(p$beta_f, p$beta_s), 0, 2.5, log = TRUE)) +
        dnorm(p$gam, 0, 1, log = TRUE) +
        dnorm(p$delta, 0, 1, log = TRUE)
  if (family == "ln") lp <- lp + dnorm(p$sp$sigma, 0, 2, log = TRUE)
  else lp <- lp + dlnorm(p$sp$a, 0, 1, log = TRUE) + dlnorm(p$sp$q, 0, 1, log = TRUE)
  -lp
}

nll_dep1 <- function(par, dat, model, family, ltrc) {
  p <- unpack(par, model, "dep1", dat$scen$family_fit %||% family)
  eta <- p$beta_f[1] + p$beta_f[2] * dat$x1 + p$beta_f[3] * dat$x2
  loc0 <- p$beta_s[1] + p$beta_s[2] * dat$x1 + p$beta_s[3] * dat$x2
  lam <- exp(pmin(eta, 30))
  if (ltrc) {
    if (abs(p$gam) < 1e-12) {
      mu <- lam * exp(sev_logS(dat$d, loc0, p$sp, family))
    } else {
      mu <- rep(0, length(lam))
      for (nn in 1:12) {
        locn <- loc0 + p$gam * (nn - lam) / sqrt(lam)
        mu <- mu + dpois(nn, lam) * nn * exp(sev_logS(dat$d, locn, p$sp, family))
      }
    }
    resid <- (dat$Nobs - mu) / sqrt(lam)
  } else {
    mu <- lam
    resid <- (dat$Nobs - mu) / sqrt(mu + 1e-12)
  }
  ll_f <- sum(dpois(dat$Nobs, pmax(mu, 1e-12), log = TRUE))
  k <- dat$py_o
  loc_k <- loc0[k] + p$gam * resid[k]
  if (ltrc) {
    lden <- ifelse(dat$cens == 1,
                   sev_logS(dat$u[k], loc_k, p$sp, family),
                   sev_logf(dat$y, loc_k, p$sp, family))
    ll_s <- sum(lden - sev_logS(dat$d[k], loc_k, p$sp, family))
  } else {
    ll_s <- sum(sev_logf(dat$y, loc_k, p$sp, family))
  }
  if (!is.finite(ll_f + ll_s)) return(1e10)
  -(ll_f + ll_s) + neg_log_prior(p, family)
}

nll_dep2 <- function(par, dat, model, family, ltrc) {
  p <- unpack(par, model, "dep2", family)
  eta <- p$beta_f[1] + p$beta_f[2] * dat$x1 + p$beta_f[3] * dat$x2
  loc0 <- p$beta_s[1] + p$beta_s[2] * dat$x1 + p$beta_s[3] * dat$x2
  n <- dat$n; K <- length(GH$nodes)
  k <- dat$py_o; m <- length(k)
  lam <- exp(pmin(outer(eta, p$delta * GH$nodes, "+"), 30))
  loc_th <- outer(loc0, p$gam * GH$nodes, "+")
  if (ltrc) {
    lS0 <- matrix(sev_logS(rep(dat$d, K), as.vector(loc_th), p$sp, family), n, K)
    mu <- lam * exp(lS0)
  } else {
    mu <- lam
  }
  ll_mat <- matrix(dpois(dat$Nobs, pmax(mu, 1e-12), log = TRUE), n, K)
  if (m > 0) {
    lock <- loc_th[k, , drop = FALSE]
    if (ltrc) {
      unc <- dat$cens == 0
      lden <- matrix(0, m, K)
      if (any(unc))
        lden[unc, ] <- matrix(sev_logf(rep(dat$y[unc], K),
                                       as.vector(lock[unc, , drop = FALSE]),
                                       p$sp, family), sum(unc), K)
      if (any(!unc))
        lden[!unc, ] <- matrix(sev_logS(rep(dat$u[k[!unc]], K),
                                        as.vector(lock[!unc, , drop = FALSE]),
                                        p$sp, family), sum(!unc), K)
      ls <- lden - lS0[k, , drop = FALSE]
    } else {
      ls <- matrix(sev_logf(rep(dat$y, K), as.vector(lock), p$sp, family), m, K)
    }
    agg <- rowsum(ls, k)
    gi <- as.integer(rownames(agg))
    ll_mat[gi, ] <- ll_mat[gi, ] + agg
  }
  ll <- sum(logsumexp_rows(sweep(ll_mat, 2, log(GH$weights), "+")))
  if (!is.finite(ll)) return(1e10)
  -ll + neg_log_prior(p, family)
}

start_vals <- function(dat, model, mech, family, warm = NULL) {
  lam_bar <- mean(dat$Nobs)
  s0 <- c(log(lam_bar), 0, 0, mean(log(dat$y)), 0, 0)
  if (model %in% c("MC", "MD")) s0 <- c(s0, 0.05)
  if (mech == "dep2") s0 <- c(s0, log(0.3))
  if (family == "ln") s0 <- c(s0, log(sd(log(dat$y))))
  else s0 <- c(s0, log(1.5), log(1.5))
  if (model == "MB" && !is.null(warm$MA)) s0 <- warm$MA
  if (model == "MC" && !is.null(warm$MA)) s0 <- append(warm$MA, 0.05, after = 6)
  if (model == "MD") {
    if (!is.null(warm$MB)) s0 <- append(warm$MB, 0.1, after = 6)
    else if (!is.null(warm$MC)) s0 <- warm$MC
  }
  s0
}

fit_one <- function(dat, model, family) {
  mech <- dat$scen$mech
  ltrc <- model %in% c("MB", "MD")
  nll <- if (mech == "dep1") nll_dep1 else nll_dep2
  s0 <- start_vals(dat, model, mech, family, warm = attr(dat, "warm"))
  o <- try(nlminb(s0, nll, dat = dat, model = model, family = family, ltrc = ltrc,
                  control = list(eval.max = 1200, iter.max = 250, rel.tol = 1e-8)),
           silent = TRUE)
  if (inherits(o, "try-error")) return(list(ok = FALSE))
  list(ok = o$convergence == 0, par = o$par, nll = o$objective,
       pars = unpack(o$par, model, mech, family))
}

premium_mc <- function(pars, mech, family, model, contracts, nsim = 1e5, seed = 7) {
  set.seed(seed)
  eta <- pars$beta_f[1]; loc0 <- pars$beta_s[1]
  gam <- pars$gam; delta <- pars$delta
  if (mech == "dep2") {
    th <- rnorm(nsim)
    lam <- exp(eta + delta * th)
    Ns <- rpois(nsim, lam)
    locv <- loc0 + gam * th
  } else {
    lam <- exp(eta)
    Ns <- rpois(nsim, lam)
    locv <- loc0 + gam * (Ns - lam) / sqrt(lam)
  }
  py <- rep(seq_len(nsim), Ns)
  X <- if (family == "ln") rlnorm(length(py), locv[py], pars$sp$sigma)
       else rgb2(length(py), pars$sp$a, pars$sp$p, pars$sp$q, exp(locv[py]))
  sapply(contracts, function(cc)
    sum(pmax(pmin(X, cc["u"]) - cc["d"], 0)) / nsim)
}

true_premium <- function(scen, contracts, nsim = 2e6) {
  pars <- list(beta_f = TRUE_PARS$beta_f, beta_s = TRUE_PARS$beta_s,
               gam = scen$gam, delta = ifelse(scen$mech == "dep2", TRUE_PARS$delta, 0),
               sp = if (scen$family == "ln") list(sigma = TRUE_PARS$sigma)
                    else list(a = TRUE_PARS$gb2["a"], p = TRUE_PARS$gb2["p"],
                              q = TRUE_PARS$gb2["q"]))
  premium_mc(pars, scen$mech, scen$family, "TRUE", contracts, nsim = nsim, seed = 4242)
}

# Data preparation ----
set.seed(SEED + 1)
pol <- fread(file.path(data_dir, "Insample.csv"))
clm <- fread(file.path(data_dir, "CLAIMLEVEL.csv"))
stopifnot(all(pol$Deduct > 0))

clm2 <- merge(clm, pol[, .(PolicyNum, Year, BCcov, LnCoverage, NoClaimCredit,
                           TypeCity, TypeCounty, TypeMisc, TypeSchool, TypeTown, TypeVillage,
                           AC00, AC05, AC10, AC15, lnDeduct, Deduct_pol = Deduct,
                           Fire5_pol = Fire5)],
              by = c("PolicyNum", "Year"), all.x = TRUE)

sub_ded  <- clm2[Claim <= Deduct, .N]
zero_clm <- clm2[Claim <= 0, .N]
cens     <- clm2[!is.na(BCcov) & Claim >= BCcov, .N]
cat("Recorded claims:", nrow(clm2), "| at or below deductible:", sub_ded,
    "| at or above coverage:", cens, "\n")

clm2[, censored := as.integer(!is.na(BCcov) & Claim >= BCcov)]
clma <- clm2[Claim > Deduct]

cnta <- clma[, .(FreqA = .N), by = .(PolicyNum, Year)]
polA <- merge(pol, cnta, by = c("PolicyNum", "Year"), all.x = TRUE)
polA[is.na(FreqA), FreqA := 0L]

polA[, EntityType := fifelse(TypeCity == 1, "City",
                     fifelse(TypeCounty == 1, "County",
                     fifelse(TypeSchool == 1, "School",
                     fifelse(TypeTown == 1, "Town",
                     fifelse(TypeVillage == 1, "Village", "Misc")))))]
clma[, EntityTypeF := EntityType]
mu_lc <- polA[Year <= 2009, mean(LnCoverage)]
sd_lc <- polA[Year <= 2009, sd(LnCoverage)]
polA[, zLnCov := (LnCoverage - mu_lc) / sd_lc]
clma <- merge(clma, polA[, .(PolicyNum, Year, zLnCov, EntityType, FreqA)],
              by = c("PolicyNum", "Year"), all.x = TRUE)
polA[, pid := as.integer(factor(PolicyNum))]
clma <- merge(clma, unique(polA[, .(PolicyNum, pid)]), by = "PolicyNum", all.x = TRUE)

dat <- list(pol = polA, clm = clma,
            train_pol = polA[Year <= 2009], hold_pol = polA[Year == 2010],
            train_clm = clma[Year <= 2009], hold_clm = clma[Year == 2010],
            std = list(mu_lc = mu_lc, sd_lc = sd_lc),
            gate = list(sub_ded = sub_ded, zero_clm = zero_clm, cens = cens,
                        n_claims_raw = nrow(clm2), n_claims_analysis = nrow(clma)))
saveRDS(dat, file.path(out_dir, "lgpif.rds"))
cat("Train policy-years:", nrow(dat$train_pol), "| holdout policy-years:", nrow(dat$hold_pol),
    "| train claims:", nrow(dat$train_clm), "| holdout claims:", nrow(dat$hold_clm), "\n")

pp <- lgpif_prepare(dat)

# Model compilation and density checks ----
set.seed(SEED + 2)
m_joint <- stan_model(model_code = joint_code, model_name = "joint_fs")
m_sev <- stan_model(model_code = sev_code, model_name = "sev_margin")
expose_stan_functions(m_sev)

dgb2R <- function(y, a, p, q, b) {
  z <- a * (log(y) - log(b))
  exp(log(a) + p * z - log(y) - lbeta(p, q) - (p + q) * log1p(exp(z)))
}
pgb2R <- function(y, a, p, q, b) pbeta(plogis(a * (log(y) - log(b))), p, q)
grid <- expand.grid(y = c(0.05, 0.5, 1, 5, 50, 500),
                    a = c(0.7, 1.5, 3), p = c(0.5, 1.2), q = c(0.8, 2.5), b = c(0.5, 4))
e_lpdf <- max(abs(mapply(function(y, a, p, q, b) gb2_lpdf(y, a, p, q, b),
                         grid$y, grid$a, grid$p, grid$q, grid$b) -
                  log(dgb2R(grid$y, grid$a, grid$p, grid$q, grid$b))))
e_cdf <- max(abs(mapply(function(y, a, p, q, b) gb2_cdf_val(y, a, p, q, b),
                        grid$y, grid$a, grid$p, grid$q, grid$b) -
                 pgb2R(grid$y, grid$a, grid$p, grid$q, grid$b)))
stopifnot(e_lpdf < 1e-8, e_cdf < 1e-8)
stopifnot(abs(integrate(function(t) dgb2R(t, 1.5, 1.2, 2.5, 4), 0, Inf,
                        rel.tol = 1e-10)$value - 1) < 1e-6)

dburrR <- function(y, a, q, b) a * q * (y / b)^(a - 1) / b / (1 + (y / b)^a)^(q + 1)
sburrR <- function(y, a, q, b) (1 + (y / b)^a)^(-q)
g <- expand.grid(y = c(0.1, 1, 10, 200), a = c(0.8, 1.6), q = c(0.9, 2.2), b = c(0.6, 5))
e1 <- max(abs(mapply(function(y, a, q, b) burr_lpdf(y, a, q, b), g$y, g$a, g$q, g$b) -
              log(dburrR(g$y, g$a, g$q, g$b))))
e2 <- max(abs(mapply(function(y, a, q, b) burr_lccdf_val(y, a, q, b), g$y, g$a, g$q, g$b) -
              log(sburrR(g$y, g$a, g$q, g$b))))
stopifnot(e1 < 1e-10, e2 < 1e-10)
cat("GB2 and Burr density checks passed\n")

# Severity margin selection ----
set.seed(SEED + 3)
fams <- c(lognormal = 1, gamma = 2, burr = 4)
sev_fits <- list(); sev_loos <- list()
for (fn in names(fams)) {
  sd_ <- lgpif_sev_data(pp, family = fams[[fn]], ltrc = 1)
  init_f <- function() list(beta_s = c(1.0, rep(0, sd_$K_s - 1)), sigma_ln = 1.5,
                            alpha_g = 0.8, gb2_a = 1.5, gb2_p = 1, gb2_q = 1.5)
  t0 <- Sys.time()
  fit <- sampling(m_sev, data = sd_, chains = 4, iter = 2000, warmup = 1000,
                  seed = SEED, refresh = 250, init = init_f,
                  control = list(adapt_delta = 0.9))
  rt <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
  ll <- extract_log_lik(fit, merge_chains = FALSE)
  lo <- loo(ll, r_eff = relative_eff(exp(ll)))
  sm <- summary(fit)$summary
  used <- c(paste0("beta_s[", 1:sd_$K_s, "]"),
            switch(fn, lognormal = "sigma_ln", gamma = "alpha_g", burr = c("gb2_a", "gb2_q")))
  print(round(sm[used, c("mean", "sd", "2.5%", "97.5%", "n_eff", "Rhat")], 3))
  sev_loos[[fn]] <- lo
  sev_fits[[fn]] <- list(summary = sm[used, ], loo = lo, runtime_min = rt,
                         divergences = get_num_divergent(fit),
                         draws = extract(fit, pars = c("beta_s", "sigma_ln", "alpha_g",
                                                       "gb2_a", "gb2_p", "gb2_q")))
}
sev_cmp <- loo_compare(sev_loos)
print(sev_cmp)
saveRDS(list(fits = sev_fits, loos = sev_loos, cmp = sev_cmp),
        file.path(out_dir, "sev_margins.rds"))

set.seed(SEED + 34)
tc <- pp$tc
y_s <- pmin(tc$Claim, tc$BCcov) / 1000
d_s <- tc$Deduct / 1000; u_s <- tc$BCcov / 1000; c_s <- tc$censored
X_s <- pp$Xs_tr_cl
K <- ncol(X_s)
nll_gb2 <- function(par, fix_p1 = FALSE) {
  beta <- par[1:K]
  a <- exp(par[K + 1])
  p <- if (fix_p1) 1 else exp(par[K + 2])
  q <- exp(par[K + 2 + !fix_p1])
  b <- exp(as.vector(X_s %*% beta))
  lS_d <- log(pmax(1 - pgb2(d_s, a, p, q, b), 1e-300))
  ll <- ifelse(c_s == 1,
               log(pmax(1 - pgb2(u_s, a, p, q, b), 1e-300)),
               dgb2(y_s, a, p, q, b, log = TRUE)) - lS_d
  pen <- sum(dnorm(beta, 0, 2.5, log = TRUE)) + dlnorm(a, 0, 1, log = TRUE) +
         (if (!fix_p1) dlnorm(p, 0, 1, log = TRUE) else 0) + dlnorm(q, 0, 1, log = TRUE)
  out <- -(sum(ll) + pen)
  if (!is.finite(out)) 1e10 else out
}
ob <- nlminb(c(1, rep(0, K - 1), log(1.6), log(0.6)), nll_gb2, fix_p1 = TRUE,
             control = list(iter.max = 500))
og <- nlminb(append(ob$par, 0, after = K + 1), nll_gb2, fix_p1 = FALSE,
             control = list(iter.max = 500))
cat(sprintf("Burr (p = 1): %.2f | GB2: %.2f | p = %.3f | improvement %.2f\n",
            ob$objective, og$objective, exp(og$par[K + 2]), ob$objective - og$objective))
saveRDS(list(burr = ob, gb2 = og, K = K, improvement = ob$objective - og$objective,
             gb2_pars = c(a = exp(og$par[K + 1]), p = exp(og$par[K + 2]),
                          q = exp(og$par[K + 3]))),
        file.path(out_dir, "gb2_ml.rds"))

# Joint frequency-severity models ----
fit_plan <- data.table(
  tag     = c("MA", "MB", "MC1", "MC2", "MD1", "MD2", "MB_NB", "MD2_NB"),
  iter    = c(2000, 3000, 2000, 3000, 3000, 3000, 2000, 3000),
  warmup  = c(1000, 1500, 1000, 1500, 1500, 1500, 1000, 1500),
  mapinit = c(FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE))
family_joint <- 4
adapt_delta <- 0.9

fit_joint <- function(tag, family, iter, warmup, adapt_delta, mapinit) {
  set.seed(SEED + 4)
  nb <- as.integer(grepl("_NB$", tag))
  spec <- switch(sub("_NB$", "", tag),
    MA  = list(dep = 0, ltrc = 0), MB  = list(dep = 0, ltrc = 1),
    MC1 = list(dep = 1, ltrc = 0), MC2 = list(dep = 2, ltrc = 0),
    MD1 = list(dep = 1, ltrc = 1), MD2 = list(dep = 2, ltrc = 1))
  sdat <- lgpif_stan_data(pp, family, spec$ltrc, spec$dep, nb)
  init_use <- joint_inits(sdat$K_f, sdat$K_s)
  if (mapinit) {
    op <- optimizing(m_joint, data = sdat, seed = SEED, as_vector = FALSE,
                     init = init_use(), iter = 4000)
    mp <- op$par
    init_use <- function() list(
      beta_f = mp$beta_f + rnorm(length(mp$beta_f), 0, 0.02),
      beta_s = mp$beta_s + rnorm(length(mp$beta_s), 0, 0.02),
      sigma_b = max(mp$sigma_b, 0.05) * exp(rnorm(1, 0, 0.05)),
      b_raw = mp$b_raw + rnorm(length(mp$b_raw), 0, 0.05),
      sigma_ln = max(mp$sigma_ln, 0.05) * exp(rnorm(1, 0, 0.05)),
      alpha_g = max(mp$alpha_g, 0.05) * exp(rnorm(1, 0, 0.05)),
      gb2_a = max(mp$gb2_a, 0.05) * exp(rnorm(1, 0, 0.05)),
      gb2_p = max(mp$gb2_p, 0.05) * exp(rnorm(1, 0, 0.05)),
      gb2_q = max(mp$gb2_q, 0.05) * exp(rnorm(1, 0, 0.05)),
      phi_nb = max(mp$phi_nb, 0.05) * exp(rnorm(1, 0, 0.05)),
      gamma_dep = mp$gamma_dep + rnorm(1, 0, 0.02))
  }
  t0 <- Sys.time()
  fit <- sampling(m_joint, data = sdat, chains = 4, iter = iter, warmup = warmup,
                  seed = SEED, refresh = 100, init = init_use,
                  control = list(adapt_delta = adapt_delta, max_treedepth = 11),
                  pars = c("b_raw", "eta_f", "log_thin", "mu_obs", "resid_N"),
                  include = FALSE)
  rt <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
  sm <- summary(fit)$summary
  key <- c(paste0("beta_f[", 1:sdat$K_f, "]"), paste0("beta_s[", 1:sdat$K_s, "]"),
           "sigma_b",
           if (family == 1) "sigma_ln" else if (family == 2) "alpha_g"
           else if (family == 4) c("gb2_a", "gb2_q") else c("gb2_a", "gb2_p", "gb2_q"),
           if (nb == 1) "phi_nb", if (spec$dep > 0) "gamma_dep")
  print(round(sm[key, c("mean", "sd", "2.5%", "97.5%", "n_eff", "Rhat")], 3))
  div <- get_num_divergent(fit)
  diag_pars <- sm[grep("^(beta_|sigma_b|b\\[|gamma_dep|sigma_ln|alpha_g|gb2_|phi_nb)",
                       rownames(sm)), ]
  drop_unused <- c(
    if (family != 1) "sigma_ln", if (family != 2) "alpha_g",
    if (family == 1 || family == 2) c("gb2_a", "gb2_q"),
    if (family != 3) "gb2_p",
    if (nb == 0) "phi_nb", if (spec$dep == 0) "gamma_dep")
  diag_pars <- diag_pars[!rownames(diag_pars) %in% drop_unused, ]
  cat(sprintf("%s: max Rhat %.4f, min ESS %.0f, divergences %d, %.1f min\n",
              tag, max(diag_pars[, "Rhat"], na.rm = TRUE),
              min(diag_pars[, "n_eff"], na.rm = TRUE), div, rt))
  dr <- extract(fit, pars = c("beta_f", "beta_s", "sigma_b", "b", "sigma_ln", "alpha_g",
                              "gb2_a", "gb2_p", "gb2_q", "phi_nb", "gamma_dep"))
  out <- list(tag = tag, family = family, spec = spec, nb = nb, draws = dr,
              log_lik = extract_log_lik(fit, merge_chains = FALSE),
              summary_key = sm[key, ], runtime_min = rt, divergences = div,
              diag = list(max_rhat = max(diag_pars[, "Rhat"], na.rm = TRUE),
                          min_neff = min(diag_pars[, "n_eff"], na.rm = TRUE)))
  saveRDS(out, file.path(fit_dir, paste0("fit_", tag, ".rds")))
  invisible(out)
}

for (i in seq_len(nrow(fit_plan))) {
  fp <- fit_plan[i]
  if (file.exists(file.path(fit_dir, paste0("fit_", fp$tag, ".rds")))) next
  fit_joint(fp$tag, family_joint, fp$iter, fp$warmup, adapt_delta, fp$mapinit)
}

fits <- list(); joint_loo <- list()
for (tg in fit_plan$tag) {
  f <- readRDS(file.path(fit_dir, paste0("fit_", tg, ".rds")))
  joint_loo[[tg]] <- loo(f$log_lik, r_eff = relative_eff(exp(f$log_lik), cores = 1), cores = 1)
  f$log_lik <- NULL
  fits[[tg]] <- f
}
joint_cmp <- loo_compare(joint_loo)
print(joint_cmp)
saveRDS(list(loos = joint_loo, cmp = joint_cmp), file.path(out_dir, "joint_loo.rds"))

# Two-part GLM benchmark ----
set.seed(SEED + 5)
tp <- pp$tp; hp <- pp$hp; tc <- copy(pp$tc); hc <- pp$hc
me_freq <- glm.nb(FreqA ~ zLnCov + TypeCity + TypeCounty + TypeSchool + TypeTown +
                    TypeVillage + NoClaimCredit + Fire5 + lnDeduct, data = tp)
tc[, ymod := pmin(Claim, BCcov) / 1000]
me_sev <- glm(ymod ~ zLnCov + TypeCity + TypeCounty + TypeSchool + TypeTown + TypeVillage,
              family = Gamma(link = "log"), data = tc)
me_disp <- summary(me_sev)$dispersion
print(summary(me_freq)$coefficients)
print(summary(me_sev)$coefficients)
saveRDS(list(freq = me_freq, sev = me_sev, disp = me_disp), file.path(fit_dir, "fit_ME.rds"))

# Out-of-sample prediction (2010) ----
set.seed(SEED + 6)
seen_pids <- unique(tp$pid)
y_hold <- {
  agg <- hc[, .(Y = sum(pmin(Claim, BCcov) / 1000)), by = py]
  v <- rep(0, nrow(hp)); v[agg$py] <- agg$Y; v
}
y_cal <- {
  cal_idx <- which(tp$Year == 2009)
  cc <- pp$tc[Year == 2009]
  agg <- cc[, .(Y = sum(pmin(Claim, BCcov) / 1000)), by = py]
  v <- rep(0, nrow(tp)); v[agg$py] <- agg$Y
  list(idx = cal_idx, y = v[cal_idx])
}
hold_claims <- list(py = hc$py, y = pmin(hc$Claim, hc$BCcov) / 1000,
                    cens = hc$censored, d = hc$Deduct / 1000, u = hc$BCcov / 1000,
                    Xs = pp$Xs_ho_cl)
cal_claims <- {
  cc <- pp$tc[Year == 2009]
  list(py = match(cc$py, y_cal$idx), y = pmin(cc$Claim, cc$BCcov) / 1000,
       cens = cc$censored, d = cc$Deduct / 1000, u = cc$BCcov / 1000,
       Xs = cbind(1, cc$zLnCov, cc$TypeCity, cc$TypeCounty, cc$TypeSchool,
                  cc$TypeTown, cc$TypeVillage))
}

results <- list()
intervals <- list()
for (tg in fit_plan$tag) {
  fit <- fits[[tg]]
  pr <- predict_model(fit, pp$Xf_ho, pp$Xs_ho_py, hp$Deduct / 1000, hp$BCcov / 1000,
                      hp$pid, hp$FreqA, hold_claims, seen_pids, S_use = 1000, seed = 601)
  logs_freq <- log_mean_exp_rows(pr$lpN)
  logs_sev <- if (!is.null(pr$lf_cl)) log_mean_exp_rows(pr$lf_cl) else numeric(0)
  crps <- crps_sample(y_hold, pr$Ymat)
  pit <- pit_from_samples(y_hold, pr$Ymat, seed = 602)
  pmean <- rowMeans(pr$Ymat)
  port <- colSums(pr$Ymat)
  prc <- predict_model(fit, pp$Xf_tr[y_cal$idx, , drop = FALSE],
                       pp$Xs_tr_py[y_cal$idx, , drop = FALSE],
                       tp$Deduct[y_cal$idx] / 1000, tp$BCcov[y_cal$idx] / 1000,
                       tp$pid[y_cal$idx], tp$FreqA[y_cal$idx], cal_claims, seen_pids,
                       S_use = 1000, seed = 603)
  Qs <- cqr_correct(y_cal$y, prc$Ymat, alphas = c(0.2, 0.1))
  cov_tab <- list()
  for (a in c(0.2, 0.1)) {
    ql <- apply(pr$Ymat, 1, quantile, probs = a / 2, names = FALSE)
    qh <- apply(pr$Ymat, 1, quantile, probs = 1 - a / 2, names = FALSE)
    Q <- Qs[[paste0("a", a * 100)]]
    cov_tab[[paste0("a", a * 100)]] <- data.table(
      alpha = a,
      cov_raw = mean(y_hold >= ql & y_hold <= qh),
      cov_conf = mean(y_hold >= pmax(ql - Q, 0) & y_hold <= qh + Q),
      wid_raw = mean(qh - ql), wid_conf = mean((qh + Q) - pmax(ql - Q, 0)),
      Q = Q)
  }
  ql <- apply(pr$Ymat, 1, quantile, probs = 0.05, names = FALSE)
  qh <- apply(pr$Ymat, 1, quantile, probs = 0.95, names = FALSE)
  Q <- Qs[["a10"]]
  intervals[[tg]] <- list(y = y_hold, ql = ql, qh = qh, Q = Q, pmean = pmean)
  bytype <- data.table(EntityType = hp$EntityType,
                       raw = (y_hold >= ql & y_hold <= qh),
                       conf = (y_hold >= pmax(ql - Q, 0) & y_hold <= qh + Q))[,
                .(cov_raw = mean(raw), cov_conf = mean(conf), n = .N), by = EntityType]
  results[[tg]] <- list(
    tag = tg,
    logs_total = sum(logs_freq) + sum(logs_sev),
    logs_freq = sum(logs_freq), logs_sev = sum(logs_sev),
    logs_mean_py = (sum(logs_freq) + sum(logs_sev)) / nrow(hp),
    crps_mean = mean(crps), pit = pit, gini = gini_lift(pmean, y_hold),
    pred_total = mean(port), realized_total = sum(y_hold),
    port_draws = port,
    var95 = quantile(port, .95), var99 = quantile(port, .99),
    tvar95 = mean(port[port >= quantile(port, .95)]),
    tvar99 = mean(port[port >= quantile(port, .99)]),
    cov = rbindlist(cov_tab), cov_bytype = bytype, pred_mean_py = pmean)
  cat(sprintf("%s: logS %.1f (freq %.1f, sev %.1f), CRPS %.3f, Gini %.3f, pred/real %.3f\n",
              tg, results[[tg]]$logs_total, results[[tg]]$logs_freq,
              results[[tg]]$logs_sev, results[[tg]]$crps_mean, results[[tg]]$gini,
              results[[tg]]$pred_total / results[[tg]]$realized_total))
}

mu_f <- predict(me_freq, newdata = hp, type = "response")
th <- me_freq$theta
loc_s <- predict(me_sev, newdata = hc, type = "link")
loc_py <- predict(me_sev, newdata = hp, type = "link")
shape <- 1 / me_disp
S_use <- 1000
set.seed(701)
Ymat <- matrix(0, nrow(hp), S_use)
for (js in seq_len(S_use)) {
  Nsim <- rnbinom(nrow(hp), mu = mu_f, size = th)
  idx <- which(Nsim > 0)
  if (length(idx)) {
    py <- rep(idx, Nsim[idx])
    X <- rgamma(length(py), shape = shape, rate = shape / exp(loc_py[py]))
    agg <- rowsum(pmin(X, hp$BCcov[py] / 1000), py)
    Ymat[as.integer(rownames(agg)), js] <- agg
  }
}
logs_freq_me <- dnbinom(hp$FreqA, mu = mu_f, size = th, log = TRUE)
logs_sev_me <- dgamma(pmin(hc$Claim, hc$BCcov) / 1000, shape = shape,
                      rate = shape / exp(loc_s), log = TRUE)
port_me <- colSums(Ymat)
results[["ME"]] <- list(tag = "ME",
  logs_total = sum(logs_freq_me) + sum(logs_sev_me),
  logs_freq = sum(logs_freq_me), logs_sev = sum(logs_sev_me),
  logs_mean_py = (sum(logs_freq_me) + sum(logs_sev_me)) / nrow(hp),
  crps_mean = mean(crps_sample(y_hold, Ymat)), pit = pit_from_samples(y_hold, Ymat, 702),
  gini = gini_lift(rowMeans(Ymat), y_hold),
  pred_total = mean(port_me), realized_total = sum(y_hold),
  port_draws = port_me,
  var95 = quantile(port_me, .95), var99 = quantile(port_me, .99),
  tvar95 = mean(port_me[port_me >= quantile(port_me, .95)]),
  tvar99 = mean(port_me[port_me >= quantile(port_me, .99)]),
  pred_mean_py = rowMeans(Ymat))
cat(sprintf("ME: logS %.1f, CRPS %.3f, pred/real %.3f\n", results$ME$logs_total,
            results$ME$crps_mean, results$ME$pred_total / results$ME$realized_total))
saveRDS(list(results = results, intervals = intervals, y_hold = y_hold,
             hp_entity = hp$EntityType),
        file.path(out_dir, "predictions.rds"))

# Pricing, ground-up recovery and dependence ----
set.seed(SEED + 7)
prof <- list(
  School  = list(cov_med = hp[TypeSchool == 1, median(BCcov)],
                 x = function(z) c(1, z, 0, 0, 1, 0, 0)),
  Village = list(cov_med = hp[TypeVillage == 1, median(BCcov)],
                 x = function(z) c(1, z, 0, 0, 0, 0, 1)))
zval <- function(cov) (log(cov / 1e6) - pp$std$mu_lc) / pp$std$sd_lc

d_grid <- c(0.5, 1, 2.5, 5, 10, 15, 25, 50, 100)
u_grid <- c(50, 100, 250, 500, 1000, 5000)
rel <- list()
for (pn in names(prof)) {
  z <- zval(prof[[pn]]$cov_med)
  xs <- prof[[pn]]$x(z)
  xf <- c(xs, 0, 0)
  for (tg in c("MB", "MD2")) {
    pc <- prem_curve(fits[[tg]], xs, xf, d_grid, u_grid,
                     u_for_d = prof[[pn]]$cov_med / 1000, d_for_u = 1,
                     seed = 11 + match(tg, c("MB", "MD2")))
    nbad <- sum(!is.finite(pc$PPd)) + sum(!is.finite(pc$PPu)) +
            sum(pc$PPd[1, ] <= 0, na.rm = TRUE) + sum(pc$PPu[3, ] <= 0, na.rm = TRUE)
    pc$PPd[!is.finite(pc$PPd)] <- NA
    pc$PPu[!is.finite(pc$PPu)] <- NA
    Rd <- sweep(pc$PPd, 2, ifelse(pc$PPd[1, ] > 0, pc$PPd[1, ], NA), "/")
    Ru <- sweep(pc$PPu, 2, ifelse(pc$PPu[3, ] > 0, pc$PPu[3, ], NA), "/")
    qq <- function(M, p) apply(M, 1, quantile, p, na.rm = TRUE)
    rel[[paste(pn, tg, sep = "_")]] <- list(
      profile = pn, tag = tg, d_grid = d_grid, u_grid = u_grid, n_bad = nbad,
      PPd_mean = rowMeans(pc$PPd, na.rm = TRUE), PPd_lo = qq(pc$PPd, .025),
      PPd_hi = qq(pc$PPd, .975),
      Rd_mean = rowMeans(Rd, na.rm = TRUE), Rd_lo = qq(Rd, .025), Rd_hi = qq(Rd, .975),
      Ru_mean = rowMeans(Ru, na.rm = TRUE), Ru_lo = qq(Ru, .025), Ru_hi = qq(Ru, .975))
  }
}

etypes <- c("City", "County", "Misc", "School", "Town", "Village")
xrow <- function(et, z) c(1, z, et == "City", et == "County", et == "School",
                          et == "Town", et == "Village")
ent_prem <- list()
for (et in etypes) {
  cov_med <- hp[EntityType == et, median(BCcov)]
  xs <- xrow(et, zval(cov_med)); xf <- c(xs, 0, 0)
  for (tg in c("MB", "MD2")) {
    pc <- prem_curve(fits[[tg]], xs, xf, d_grid = 1, u_grid = 250,
                     u_for_d = cov_med / 1000, d_for_u = 1, S_use = 500,
                     M = 4000, seed = 77)
    pd <- pc$PPd[1, ]
    pd[!is.finite(pd)] <- NA
    ent_prem[[paste(et, tg, sep = "_")]] <-
      data.table(EntityType = et, model = tg, cov_med = cov_med,
                 pp_mean = mean(pd, na.rm = TRUE), pp_lo = quantile(pd, .025, na.rm = TRUE),
                 pp_hi = quantile(pd, .975, na.rm = TRUE))
  }
}
ent_prem <- rbindlist(ent_prem)
print(ent_prem)

fit <- fits$MD2
dr <- fit$draws; fam <- fit$family
ss <- round(seq(1, length(dr$sigma_b), length.out = 500))
S_thin <- matrix(0, nrow(tp), length(ss))
lam_gu <- matrix(0, nrow(tp), length(ss))
for (js in seq_along(ss)) {
  s <- ss[js]
  sp <- get_sp(fit, s)
  b_v <- dr$b[s, tp$pid]
  loc <- as.vector(pp$Xs_tr_py %*% dr$beta_s[s, ]) + dr$gamma_dep[s] * b_v
  S_thin[, js] <- exp(sev_logS_v(tp$Deduct / 1000, loc, fam, sp))
  lam_gu[, js] <- exp(as.vector(pp$Xf_tr %*% dr$beta_f[s, ]) + b_v)
}
gu_total_draws <- colSums(lam_gu)
obs_total_draws <- colSums(lam_gu * S_thin)
actual_all <- clm[Year <= 2009, .N]
actual_obs <- sum(tp$FreqA)
gu_check <- list(
  post_gu_mean = mean(gu_total_draws), post_gu_lo = quantile(gu_total_draws, .025),
  post_gu_hi = quantile(gu_total_draws, .975),
  post_obs_mean = mean(obs_total_draws),
  actual_all_recorded = actual_all, actual_observed = actual_obs,
  subded_share_pred = 1 - mean(obs_total_draws) / mean(gu_total_draws),
  subded_share_actual = 1 - actual_obs / actual_all)
cat(sprintf("Ground-up claims 2006-2009: %.0f [%.0f, %.0f]; recorded %d; expected observed %.0f; observed %d\n",
            gu_check$post_gu_mean, gu_check$post_gu_lo, gu_check$post_gu_hi,
            actual_all, gu_check$post_obs_mean, actual_obs))

dep_post <- lapply(fits[c("MC1", "MD1", "MC2", "MD2")], function(f) as.vector(f$draws$gamma_dep))
saveRDS(list(rel = rel, ent_prem = ent_prem, thin_mean = rowMeans(S_thin),
             lam_mean = rowMeans(lam_gu), gu_check = gu_check, dep_post = dep_post,
             d_grid = d_grid, u_grid = u_grid),
        file.path(out_dir, "pricing.rds"))

# Simulation study ----
SCEN_GRID <- expand.grid(mech = c("dep1", "dep2"), family = c("ln", "gb2"),
                         gam_lvl = c("none", "mod", "strong"),
                         dist = c("lgpif", "heavy"), stringsAsFactors = FALSE)
SCEN_GRID$gam <- with(SCEN_GRID, ifelse(mech == "dep1",
                      c(none = 0, mod = 0.15, strong = 0.3)[gam_lvl],
                      c(none = 0, mod = 0.3, strong = 0.6)[gam_lvl]))
SCEN_GRID$R <- with(SCEN_GRID, ifelse(mech == "dep1",
                    ifelse(family == "ln", ifelse(gam_lvl == "none", 100, 60),
                                           ifelse(gam_lvl == "none", 60, 40)),
                    ifelse(family == "ln", 25, 15)))
SCEN_GRID$cell <- seq_len(nrow(SCEN_GRID))
sim_dir <- file.path(out_dir, "simulation")
dir.create(sim_dir, showWarnings = FALSE)

for (ci in SCEN_GRID$cell) {
  cell_file <- file.path(sim_dir, sprintf("cell_%02d.rds", ci))
  if (file.exists(cell_file)) next
  sc <- SCEN_GRID[ci, ]
  scen <- list(mech = sc$mech, family = sc$family, gam = sc$gam, dist = sc$dist)
  qs <- sev_quantiles(scen$family, c(.25, .30, .50, .90, .95, .99))
  contracts <- list(C1 = c(d = qs[[1]], u = qs[[6]]),
                    C2 = c(d = qs[[3]], u = qs[[5]]),
                    C3 = c(d = qs[[2]], u = qs[[4]]))
  pp_true <- true_premium(scen, contracts)
  out <- vector("list", sc$R)
  for (r in seq_len(sc$R)) {
    seed <- SEED + ci * 1000 + r
    sdat <- gen_data(scen, seed)
    warm <- list()
    res <- list(cell = ci, rep = r, seed = seed,
                mech = sc$mech, family = sc$family, gam_lvl = sc$gam_lvl,
                gam_true = sc$gam, dist = sc$dist,
                n_obs_claims = length(sdat$y), cens_share = mean(sdat$cens),
                pp_true = pp_true)
    for (m in c("MA", "MB", "MC", "MD")) {
      attr(sdat, "warm") <- warm
      t1 <- Sys.time()
      f <- fit_one(sdat, m, scen$family)
      el <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
      if (f$ok) warm[[m]] <- f$par
      pr <- if (f$ok) premium_mc(f$pars, sc$mech, scen$family, m, contracts,
                                 nsim = 1e5, seed = seed + 17) else rep(NA_real_, 3)
      res[[m]] <- list(ok = f$ok, t = el,
                       gam = if (f$ok) f$pars$gam else NA,
                       delta = if (f$ok) f$pars$delta else NA,
                       bf = if (f$ok) f$pars$beta_f else rep(NA, 3),
                       bs = if (f$ok) f$pars$beta_s else rep(NA, 3),
                       tail = if (f$ok) { if (scen$family == "ln") f$pars$sp$sigma
                                          else f$pars$sp$q } else NA,
                       pp = pr)
    }
    out[[r]] <- res
  }
  saveRDS(out, cell_file)
  cat(sprintf("Cell %d (%s, %s, %s, %s) finished, R = %d\n", ci, sc$mech, sc$family,
              sc$gam_lvl, sc$dist, sc$R))
}

rows <- list()
for (f in list.files(sim_dir, pattern = "^cell_\\d+\\.rds$", full.names = TRUE)) {
  for (res in readRDS(f)) {
    for (m in c("MA", "MB", "MC", "MD")) {
      fm <- res[[m]]
      if (is.null(fm)) next
      rows[[length(rows) + 1]] <- data.table(
        cell = res$cell, rep = res$rep, mech = res$mech, family = res$family,
        gam_lvl = res$gam_lvl, gam_true = res$gam_true, dist = res$dist,
        model = m, ok = fm$ok, t = fm$t,
        gam = fm$gam, delta = fm$delta, tail = fm$tail,
        bf1 = fm$bf[2], bs1 = fm$bs[2],
        pp1 = fm$pp[1], pp2 = fm$pp[2], pp3 = fm$pp[3],
        pt1 = res$pp_true[1], pt2 = res$pp_true[2], pt3 = res$pp_true[3])
    }
  }
}
sim <- rbindlist(rows)
cat("Simulation fits:", nrow(sim), "| convergence rate:", round(mean(sim$ok), 4), "\n")
sim <- sim[ok == TRUE]
tail_true <- 1.5
sim_summ <- sim[, .(
  R = .N,
  gam_bias = mean(gam) - gam_true[1], gam_rmse = sqrt(mean((gam - gam_true[1])^2)),
  gam_se = sd(gam) / sqrt(.N),
  tail_bias = mean(tail) - tail_true, tail_rmse = sqrt(mean((tail - tail_true)^2)),
  bf1_bias = mean(bf1) - 0.3, bf1_rmse = sqrt(mean((bf1 - 0.3)^2)),
  bs1_bias = mean(bs1) - 0.4,
  pp1_rb = mean(pp1 / pt1 - 1), pp2_rb = mean(pp2 / pt2 - 1),
  pp3_rb = mean(pp3 / pt3 - 1),
  pp1_rb_se = sd(pp1 / pt1) / sqrt(.N), pp2_rb_se = sd(pp2 / pt2) / sqrt(.N),
  pp3_rb_se = sd(pp3 / pt3) / sqrt(.N)),
  by = .(mech, family, gam_lvl, gam_true, dist, model)]
setorder(sim_summ, mech, family, dist, gam_lvl, model)
print(sim_summ)
saveRDS(list(raw = sim, summ = sim_summ), file.path(out_dir, "sim_summary.rds"))

# Simulation check against MCMC ----
set.seed(SEED + 8)
m_sim <- stan_model(model_code = sim_dep2_code, model_name = "sim_dep2")
scen <- list(mech = "dep2", family = "ln", gam = 0.6, dist = "heavy")
CELL_ID <- 22
mc_res <- list()
for (r in 1:5) {
  seed <- SEED + CELL_ID * 1000 + r
  sdat <- gen_data(scen, seed)
  ml <- fit_one(sdat, "MD", "ln")
  stan_dat <- list(n = sdat$n, N = sdat$Nobs, x1 = sdat$x1, x2 = sdat$x2,
                   d = sdat$d, u = sdat$u, n_cl = length(sdat$y), y = sdat$y,
                   cens = sdat$cens, py_of_cl = sdat$py_o, ltrc = 1, depon = 1)
  t0 <- Sys.time()
  fit <- sampling(m_sim, data = stan_dat, chains = 4, iter = 800, warmup = 400,
                  seed = seed, refresh = 0,
                  init = function() list(beta_f = c(-0.3, 0.3, 0.2),
                                         beta_s = c(1, 0.4, -0.3), delta = 0.4,
                                         gam = 0.3, sigma = 1.5,
                                         theta_raw = rep(0, stan_dat$n)),
                  control = list(adapt_delta = 0.9))
  rt <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  sm <- summary(fit, pars = c("beta_f", "beta_s", "delta", "gam", "sigma"))$summary
  mc_res[[r]] <- data.frame(rep = r, ml_gam = ml$pars$gam, ml_sigma = ml$pars$sp$sigma,
                            ml_delta = ml$pars$delta, ml_bf1 = ml$pars$beta_f[2],
                            mc_gam = sm["gam", "mean"], mc_gam_sd = sm["gam", "sd"],
                            mc_sigma = sm["sigma", "mean"], mc_delta = sm["delta", "mean"],
                            mc_bf1 = sm["beta_f[2]", "mean"],
                            max_rhat = max(sm[, "Rhat"]), min_neff = min(sm[, "n_eff"]),
                            div = get_num_divergent(fit), rt_min = rt)
}
mc_res <- do.call(rbind, mc_res)
print(mc_res)
cat(sprintf("Mean |ML - MCMC| for gamma: %.4f (posterior sd %.4f)\n",
            mean(abs(mc_res$ml_gam - mc_res$mc_gam)), mean(mc_res$mc_gam_sd)))
saveRDS(mc_res, file.path(out_dir, "sim_mcmc_check.rds"))

# Figures ----
set.seed(SEED + 9)
thm <- theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        strip.background = element_rect(fill = "grey93"))
save_fig <- function(p, name, w = 6.5, h = 4) {
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = w, height = h,
         device = cairo_pdf)
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
}
PAL <- c(`M-A` = "#E69F00", `M-B` = "#56B4E9", `M-C` = "#009E73", `M-D` = "#0072B2")
status_cols <- c("Below deductible (fund-recorded only)" = "#E69F00",
                 "Insurer-observed" = "#0072B2",
                 "At/above coverage (censored)" = "#D55E00")
entity_cols <- c(City = "#E69F00", County = "#56B4E9", Misc = "#999999",
                 School = "#009E73", Town = "#D55E00", Village = "#CC79A7")

polf <- dat$pol
clm_all <- merge(copy(clm), polf[, .(PolicyNum, Year, BCcov, EntityType, Deduct_p = Deduct)],
                 by = c("PolicyNum", "Year"), all.x = FALSE)
clm_all[, status := fifelse(Claim <= Deduct, "Below deductible (fund-recorded only)",
                    fifelse(Claim >= BCcov, "At/above coverage (censored)",
                            "Insurer-observed"))]

pA <- ggplot(clm_all, aes(Claim, fill = status)) +
  geom_histogram(bins = 60, color = NA) +
  scale_x_log10(labels = label_dollar(scale_cut = cut_short_scale())) +
  scale_fill_manual(values = status_cols, name = NULL) +
  labs(x = "Recorded ground-up claim amount (log scale)", y = "Claims",
       subtitle = "(a) All fund-recorded claims, 2006-2010") +
  thm + guides(fill = guide_legend(nrow = 3))
set.seed(42); samp <- clm_all[sample(.N, 3000)]
pB <- ggplot(samp, aes(Deduct, Claim, color = status)) +
  geom_jitter(size = 0.5, alpha = 0.5, width = 0.04, height = 0) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.4) +
  scale_x_log10(labels = label_dollar()) +
  scale_y_log10(labels = label_dollar(scale_cut = cut_short_scale())) +
  scale_color_manual(values = status_cols, guide = "none") +
  labs(x = "Policy deductible (log scale)", y = "Claim amount (log scale)",
       subtitle = "(b) Truncation frontier: claim vs deductible") +
  thm
save_fig(pA + pB + plot_layout(guides = "collect") & theme(legend.position = "bottom"),
         "fig_data_distortion", 6.5, 3.4)

pA <- ggplot(polf, aes(factor(Deduct / 1000), fill = EntityType)) +
  geom_bar(position = "stack") +
  scale_fill_manual(values = entity_cols, name = NULL) +
  labs(x = "Deductible ($000)", y = "Policy-years", subtitle = "(a) Deductibles") + thm +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
pB <- ggplot(polf, aes(EntityType, BCcov / 1e6, fill = EntityType)) +
  geom_boxplot(outlier.size = 0.4, linewidth = 0.3, show.legend = FALSE) +
  scale_y_log10(labels = label_dollar(suffix = "M")) +
  scale_fill_manual(values = entity_cols) +
  labs(x = NULL, y = "Coverage amount (log scale)", subtitle = "(b) Coverage amounts") + thm +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(pA + pB + plot_layout(guides = "collect") & theme(legend.position = "bottom"),
         "fig_deduct_coverage_dist", 6.5, 3.6)

eb <- sim[model %in% c("MC", "MD") & gam_lvl != "none"]
eb[, model_lab := fifelse(model == "MC", "M-C (uncorrected)", "M-D (LTRC)")]
eb[, mech_lab := fifelse(mech == "dep1", "DEP-1 truth", "DEP-2 truth")]
eb[, cell_lab := paste0(fifelse(family == "ln", "LN", "Burr"), ", ",
                        fifelse(dist == "lgpif", "LGPIF-matched", "heavy"))]
eb[, glab := factor(gam_lvl, c("mod", "strong"), c("moderate", "strong"))]
tru <- unique(eb[, .(mech_lab, cell_lab, glab, gam_true)])
p <- ggplot(eb, aes(glab, gam, fill = model_lab)) +
  geom_boxplot(outlier.size = 0.3, linewidth = 0.3, position = position_dodge(0.75)) +
  geom_segment(data = tru, aes(x = as.integer(glab) - 0.38, xend = as.integer(glab) + 0.38,
                               y = gam_true, yend = gam_true),
               inherit.aes = FALSE, linetype = 2, linewidth = 0.5) +
  facet_grid(mech_lab ~ cell_lab, scales = "free_y") +
  scale_fill_manual(values = c("M-C (uncorrected)" = "#009E73", "M-D (LTRC)" = "#0072B2"),
                    name = NULL) +
  labs(x = "True dependence strength", y = "Estimated dependence parameter") + thm
save_fig(p, "fig_simA_dependence_bias", 6.5, 4.6)

e5 <- sim_summ[gam_lvl == "strong"]
long <- rbindlist(lapply(1:3, function(ci)
  data.table(mech = e5$mech, family = e5$family, dist = e5$dist, model = e5$model,
             contract = paste0("C", ci), rb = e5[[paste0("pp", ci, "_rb")]],
             se = e5[[paste0("pp", ci, "_rb_se")]])))
long[, model_lab := paste0("M-", substr(model, 2, 2))]
long[, mech_lab := fifelse(mech == "dep1", "DEP-1 truth", "DEP-2 truth")]
long[, cell_lab := paste0(fifelse(family == "ln", "LN", "Burr"), ", ",
                          fifelse(dist == "lgpif", "LGPIF-matched", "heavy"))]
p <- ggplot(long, aes(contract, 100 * rb, color = model_lab)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4) +
  geom_pointrange(aes(ymin = 100 * (rb - 2 * se), ymax = 100 * (rb + 2 * se)),
                  position = position_dodge(0.6), size = 0.25, linewidth = 0.5) +
  facet_grid(mech_lab ~ cell_lab, scales = "free_y") +
  scale_color_manual(values = PAL, name = NULL) +
  labs(x = "Contract (d, u)", y = "Pure-premium relative bias (%)") + thm
save_fig(p, "fig_simA_premium_bias", 6.5, 4.6)

pc <- readRDS(file.path(out_dir, "pricing.rds"))
dd <- data.table(thin = pc$thin_mean)
pA <- ggplot(dd, aes(thin)) +
  geom_histogram(bins = 50, fill = "#0072B2", color = NA) +
  labs(x = expression(paste("Posterior mean ", S(d ~ "|" ~ x))), y = "Policy-years",
       subtitle = "(a) Estimated claim-observation probability") + thm
gc <- pc$gu_check
lv <- c("Recorded (all)", "Model ground-up", "Observed (> d)", "Model observed")
bars <- data.table(
  what = factor(lv, lv),
  val = c(gc$actual_all_recorded, gc$post_gu_mean, gc$actual_observed, gc$post_obs_mean),
  lo = c(NA, gc$post_gu_lo, NA, NA), hi = c(NA, gc$post_gu_hi, NA, NA),
  src = c("Actual", "M-D2 posterior", "Actual", "M-D2 posterior"))
pB <- ggplot(bars, aes(what, val, fill = src)) +
  geom_col(width = 0.65) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2, linewidth = 0.4, na.rm = TRUE) +
  scale_fill_manual(values = c(Actual = "#999999", `M-D2 posterior` = "#0072B2"), name = NULL) +
  labs(x = NULL, y = "Claim count, 2006-2009",
       subtitle = "(b) Ground-up recovery vs fund records") + thm +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_fig(pA + pB, "fig_thinning_effect", 6.5, 3.4)

dp <- rbindlist(lapply(names(pc$dep_post), function(tg)
  data.table(tag = tg, gamma = pc$dep_post[[tg]])))
dp[, mech := fifelse(grepl("1$", tag), "DEP-1 (count residual)", "DEP-2 (shared random effect)")]
dp[, corr := fifelse(grepl("^MC", tag), "M-C (uncorrected)", "M-D (LTRC)")]
p <- ggplot(dp, aes(gamma, color = corr, fill = corr)) +
  geom_density(alpha = 0.25, linewidth = 0.6) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.4) +
  facet_wrap(~mech, scales = "free") +
  scale_color_manual(values = c("M-C (uncorrected)" = "#009E73", "M-D (LTRC)" = "#0072B2"),
                     name = NULL) +
  scale_fill_manual(values = c("M-C (uncorrected)" = "#009E73", "M-D (LTRC)" = "#0072B2"),
                    name = NULL) +
  labs(x = expression(paste("Dependence parameter ", gamma)), y = "Posterior density") + thm
save_fig(p, "fig_lgpif_dependence_posterior", 6.5, 3.2)

prd <- readRDS(file.path(out_dir, "predictions.rds"))
pit <- rbindlist(lapply(c("MA", "MC2", "MD2"), function(tg)
  data.table(model = sub("M(.)(\\d?)", "M-\\1\\2", tg), pit = prd$results[[tg]]$pit)))
pit[, model := factor(model, c("M-A", "M-C2", "M-D2"))]
p <- ggplot(pit, aes(pit)) +
  geom_histogram(breaks = seq(0, 1, by = 0.05), aes(y = after_stat(density)),
                 fill = "#0072B2", color = "white", linewidth = 0.2) +
  geom_hline(yintercept = 1, linetype = 2, linewidth = 0.4) +
  facet_wrap(~model) +
  labs(x = "Randomized PIT of 2010 aggregate policy-year loss", y = "Density") + thm
save_fig(p, "fig_pit_2010", 6.5, 2.6)

mk <- function(key, what) {
  r <- pc$rel[[key]]
  data.table(x = if (what == "d") r$d_grid else r$u_grid,
             m = if (what == "d") r$Rd_mean else r$Ru_mean,
             lo = if (what == "d") r$Rd_lo else r$Ru_lo,
             hi = if (what == "d") r$Rd_hi else r$Ru_hi,
             model = fifelse(r$tag == "MB", "M-B", "M-D2"), profile = r$profile)
}
rel_scales <- list(
  scale_linetype_manual(values = c(`M-B` = "22", `M-D2` = "solid"), name = NULL),
  scale_color_manual(values = c(`M-B` = "#56B4E9", `M-D2` = "#0072B2"), name = NULL),
  scale_fill_manual(values = c(`M-B` = "#56B4E9", `M-D2` = "#0072B2"), name = NULL))
dd <- rbindlist(lapply(names(pc$rel), mk, what = "d"))
du <- rbindlist(lapply(names(pc$rel), mk, what = "u"))
pA <- ggplot(dd, aes(x, m, color = model, fill = model, linetype = model)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, color = NA) +
  geom_line(linewidth = 0.6) + facet_wrap(~profile) + scale_x_log10() + rel_scales +
  labs(x = "Deductible ($000, log scale)", y = "Relativity vs d = $500",
       subtitle = "(a) Deductible relativities") + thm
pB <- ggplot(du, aes(x, m, color = model, fill = model, linetype = model)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, color = NA) +
  geom_line(linewidth = 0.6) + facet_wrap(~profile) +
  scale_x_log10(labels = label_comma()) + rel_scales +
  labs(x = "Per-occurrence limit ($000, log scale)", y = "Relativity vs u = $250k",
       subtitle = "(b) Limit relativities") + thm
save_fig(pA / pB + plot_layout(guides = "collect") & theme(legend.position = "bottom"),
         "fig_relativity_curves", 6.5, 5.6)

ii <- prd$intervals$MD2
set.seed(SEED)
pick <- sort(sample(seq_along(ii$y), 50))
dd <- data.table(y = ii$y[pick], ql = ii$ql[pick], qh = ii$qh[pick],
                 qlc = pmax(ii$ql[pick] - ii$Q, 0), qhc = ii$qh[pick] + ii$Q,
                 pm = ii$pmean[pick])
dd <- dd[order(pm)]; dd[, idx := .I]
dd[, covered := y >= qlc & y <= qhc]
p <- ggplot(dd, aes(idx)) +
  geom_linerange(aes(ymin = qlc + 0.05, ymax = qhc + 0.05), color = "#56B4E9",
                 linewidth = 1.6, alpha = 0.5) +
  geom_linerange(aes(ymin = ql + 0.05, ymax = qh + 0.05), color = "#0072B2", linewidth = 0.7) +
  geom_point(aes(y = y + 0.05, shape = covered), size = 1.1, color = "black") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 4),
                     labels = c(`TRUE` = "covered", `FALSE` = "missed"), name = NULL) +
  scale_y_log10(labels = label_dollar(suffix = "k")) +
  labs(x = "50 random 2010 policy-years (sorted by predicted mean)",
       y = "Aggregate loss + $50 (log scale)",
       subtitle = "Dark: raw 90% posterior predictive. Light: conformal-calibrated.") + thm
save_fig(p, "fig_predictive_intervals_2010", 6.5, 3.6)

sel <- c("MA", "MB", "MC2", "MD2")
dd <- rbindlist(lapply(sel, function(tg)
  data.table(model = sub("M(.)(\\d?)", "M-\\1\\2", tg),
             tot = prd$results[[tg]]$port_draws / 1000)))
dd[, model := factor(model, c("M-A", "M-B", "M-C2", "M-D2"))]
realized <- prd$results[[sel[1]]]$realized_total / 1000
p <- ggplot(dd, aes(tot, color = model)) +
  geom_density(linewidth = 0.7) +
  geom_vline(xintercept = realized, linetype = 2, linewidth = 0.5) +
  annotate("text", x = realized, y = Inf, label = "  realized 2010",
           hjust = 0, vjust = 1.5, size = 2.8) +
  scale_color_manual(values = c(`M-A` = "#E69F00", `M-B` = "#56B4E9",
                                `M-C2` = "#009E73", `M-D2` = "#0072B2"), name = NULL) +
  labs(x = "2010 portfolio aggregate observed loss ($M)",
       y = "Posterior predictive density") + thm
save_fig(p, "fig_portfolio_risk", 6.5, 3.4)

# Session information ----
writeLines(c(paste("Stan (via rstan):", as.character(stan_version())),
             capture.output(sessionInfo())),
           file.path(out_dir, "sessionInfo.txt"))
