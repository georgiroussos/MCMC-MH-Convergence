
## MCMC convergence: MH sampler & convergence diagnostics ------------------

# Example from Lesaffre & Lawson's osteoporosis regression (page 97)
# Total body bone mineral content (TBBMC) ~ BMI 

# Model: Bayesian simple linear regression   y_i = b0 + b1*x_i + e_i,
# [e_i ~ N(0, sigma^2)]
# Parameters sampled: b0, b1, log(sigma)   [keeps sigma > 0]

# Example shows low mixing for the regression coefficients, fast mixing for the 
# variance 

library(coda)
set.seed(17082026)

# Simulate data -----------------------------------------------------------

n <- 100
true_b0 <- 2.0
true_b1 <- 1.5
true_sig <- 1.0

x <- runif(n, 0, 10)
y <- true_b0 + true_b1 * x + rnorm(n, 0, true_sig)



# Calc log-likelihood -----------------------------------------------------
# log-likelihood + log-prior
# theta = c(b0, b1, log_sigma)

log_posterior <- function(theta, x, y) {
  b0 <- theta[1]
  b1 <- theta[2]
  log_sigma <- theta[3]
  sigma <- exp(log_sigma)
  
  mu <- b0 + b1 * x
  loglik <- sum(dnorm(y, mean = mu, sd = sigma, log = TRUE))
  
  # weekly informative priors:
  # b0, b1 ~ N(0, 10^2), and, 
  # log(sigma) ~ N(0, 2^2)
  logprior <- dnorm(b0, 0, 10, log = TRUE) +
    dnorm(b1, 0, 10, log = TRUE) +
    dnorm(log_sigma, 0, 2, log = TRUE)
  
  loglik + logprior
}



# Metropolis-Hastings (MH) sampler  ---------------------------------------
# theta* = theta_{t-1} + rnorm(0, step) 
# [remember theta = c(b0, b1, log_sigma)]

run_mh_chain <- function(start, n_iter, step_sd, x, y) {
  n_par <- length(start)
  chain <- matrix(NA_real_, nrow = n_iter, ncol = n_par)
  chain[1, ] <- start
  
  cur_lp <- log_posterior(start, x, y)
  n_accept <- 0
  
  for (t in 2:n_iter) {
    proposal <- chain[t - 1, ] + rnorm(n_par, 0, step_sd)
    prop_lp <- log_posterior(proposal, x, y)
    
    log_alpha <- prop_lp - cur_lp
    if (log(runif(1)) < log_alpha) {
      chain[t, ] <- proposal
      cur_lp <- prop_lp
      n_accept <- n_accept + 1
    } else {
      chain[t, ] <- chain[t - 1, ]
    }
  }
  attr(chain, "accept_rate") <- n_accept / (n_iter - 1)
  chain
}


# Is your chain running? (You better catch it!) ---------------------------
# Lol I'm so funny :D

# Run 4 chains with 8000 iterations each and a burnin of the first 2000 
# iterations; we want the starting points to be deliberately overdispersed 

n_iter  <- 8000
burnin  <- 2000
step_sd <- c(0.16, 0.04, 0.06)   # step size per parameter (b0, b1, log_sigma)

start_points <- list(
  c(b0 = -5, b1 = -1, log_sigma = 2),
  c(b0 = 10, b1 =  4, log_sigma = -2),
  c(b0 =  0, b1 =  0, log_sigma = 0),
  c(b0 =  15, b1 = -3, log_sigma = 1.5)
)

raw_chains <- lapply(start_points, function(s0) {
  run_mh_chain(start = s0, n_iter = n_iter, step_sd = step_sd, x = x, y = y)
})

# Acceptance rates by chain:
print(sapply(raw_chains, function(ch) round(attr(ch, "accept_rate"), 3)))
# (Target ~20-40% for random-walk MH -> retune step_sd if far outside that)



# CODA  -------------------------------------------------------------------
# raw chains as coda mcmc objects (discarding burnin), then an mcmc.list 
# for the multi-chain diagnostic


to_mcmc <- function(chain_matrix) {
  post <- chain_matrix[(burnin + 1):n_iter, , drop = FALSE]
  post <- cbind(post, sigma = exp(post[, 3]))
  colnames(post) <- c("b0", "b1", "log_sigma", "sigma")
  coda::mcmc(post, start = burnin + 1, end = n_iter, thin = 1)
}

mcmc_chains <- lapply(raw_chains, to_mcmc)
mcmc_multi <- coda::mcmc.list(mcmc_chains)
chain1 <- mcmc_chains[[1]] 


# Graphical diagnostics ---------------------------------------------------

#### Trace plots (all 4 chains) ---------------------------------------------

# the coda default plot() on an mcmc.list gives both side by side, one row 
# per parameter
plot(mcmc_multi, col = c("#61223B", "#B79961", "#4F9C82", "#CE3F27"), ask = FALSE)

traceplot(mcmc_multi, col = c("#61223B", "#B79961", "#4F9C82", "#CE3F27"), ask = FALSE)


#### Running mean (ergodic mean) plot ---------------------------------------------
# should stabilise once stationary

cumuplot(chain1, probs = c(0.025, 0.5, 0.975), ask = FALSE)


#### Autocorrelation plots ---------------------------------------------
# the ACF/"lag" plot, one panel per parameter

autocorr.plot(chain1, ask = FALSE)


#### Cross-correlation plots ---------------------------------------------

crosscorr(chain1)
crosscorr.plot(chain1)


# Formal diagnostics ------------------------------------------------------

#### Geweke diagnostic (single chain) ---------------------------------------------

# compares the mean of the first 10% to the mean of the last 50%, 
# in a Z-score that accounts for autocorrelation; 
# |Z| > ~1.96 flags nonstationarity

geweke.diag(chain1, frac1 = 0.1, frac2 = 0.5)

geweke.plot(chain1, ask = FALSE)


#### Heidelberger-Welch diagnostic (single chain) ---------------------------------
print(heidel.diag(chain1))


#### Raftery-Lewis diagnostic (single chain) ---------------------------------
# number of iterations needed for an accurate estimate of a target quantile 
# (default: the 2.5% quantile, i.e. 95% credible interval's lower bound

print(raftery.diag(chain1, q = 0.025, r = 0.005, s = 0.95))


#### Gelman-Rubin / BGR diagnostic (all 4 chains) ---------------------------------

print(gelman.diag(mcmc_multi, autoburnin = FALSE))

gelman.plot(mcmc_multi, autoburnin = FALSE, ask = FALSE)




# ESS and MCSE ------------------------------------------------------------
# effective sample size (ESS) and Monte Carlo standard error (MCSE)

print(t(sapply(mcmc_chains, effectiveSize)))

print(effectiveSize(mcmc_multi))


mcse_timeseries <- function(x) {
  sqrt(coda::spectrum0.ar(x)$spec / length(x))
}
combined <- as.matrix(mcmc_multi)


mcse <- apply(combined, 2, mcse_timeseries)
post_sd <- apply(combined, 2, sd)
print(data.frame(
  MCSE = round(mcse, 5),
  posterior_SD = round(post_sd, 4),
  `MCSE_as_percent_of_SD` = round(100 * mcse / post_sd, 2)
))



# Posterior dist. summary -------------------------------------------------
# all chains:
print(summary(mcmc_multi))

# vs. true vales:
print(c(b0 = true_b0, b1 = true_b1, sigma = true_sig))

