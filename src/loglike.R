library(MCMCpack)

# Specify likelihoods of parameters for individual models ---------------------

# Independant parallel self terminating
ll_IST <- function(x, data) { # nolint
  ydat <- data[data$response == 2, ]
  ndat <- data[data$response == 1, ]
  A <- x["A"] # nolint
  t0 <- x["t0"]

  yes <- 2 *
    dlba_norm(ydat$rt, A, x["b_pos"], t0, x[ydat$v_pos], 1) *
    (1 - plba_norm(ydat$rt, A, x["b_pos"], t0, x[ydat$v_pos], 1)) *
    (1 - plba_norm(ydat$rt, A, x["b_neg"], t0, x[ydat$v_neg], 1)**2)
  no <- 2 *
    dlba_norm(ndat$rt, A, x["b_neg"], t0, x[ndat$v_neg], 1) *
    plba_norm(ndat$rt, A, x["b_neg"], t0, x[ndat$v_neg], 1) *
    (1 - plba_norm(ndat$rt, A, x["b_pos"], t0, x[ndat$v_pos], 1))**2
  c(yes, no)
}

sample_ll_IST <- function(x, data) {
  data$response <- NA
  data$rt <- NA

  yesrt <- rlba_norm(nrow(data), x["A"], x["b_pos"], x["t0"], x[data$v_pos], 1)
  nort <- rlba_norm(nrow(data), x["A"], x["b_neg"], x["t0"], x[data$v_neg], 1)
  data$rt <- pmin(yesrt, nort)
  data$response <- ifelse(yesrt < nort, 1, 2)
}

ll_IEX <- function(x, data) { # nolint
  ydat <- data[data$response == 2, ]
  ndat <- data[data$response == 1, ]
  A <- x["A"] # nolint
  t0 <- x["t0"]
  yes <- 2 *
    dlba_norm(ydat$rt, A, x["b_pos"], t0, x[ydat$v_pos], 1) *
    plba_norm(ydat$rt, A, x["b_pos"], t0, x[ydat$v_pos], 1) *
    (1 - plba_norm(ydat$rt, A, x["b_neg"], t0, x[ydat$v_neg], 1))**2
  no <- 2 *
    dlba_norm(ndat$rt, A, x["b_neg"], t0, x[ndat$v_neg], 1) *
    (1 - plba_norm(ndat$rt, A, x["b_neg"], t0, x[ndat$v_neg], 1)) *
    (1 - plba_norm(ndat$rt, A, x["b_pos"], t0, x[ndat$v_pos], 1)**2)
  c(yes, no)
}

ll_CYST <- function(x, data) { # nolint
  ydat <- data[data$response == 2, ]
  ndat <- data[data$response == 1, ]
  A <- x["A"] # nolint
  t0 <- x["t0"]

  yes <- dlba_norm(ydat$rt, 2 * A, 2 * x["b_pos"], t0, 2 * x[ydat$v_pos], sqrt(2)) * # nolint
    (1 - plba_norm(ydat$rt, A, x["b_neg"], t0, x[ydat$v_neg], 1)**2)
  no <- 2 *
    dlba_norm(ndat$rt, A, x["b_neg"], t0, x[ndat$v_neg], 1) *
    plba_norm(ndat$rt, A, x["b_neg"], t0, x[ndat$v_neg], 1) *
    (1 - plba_norm(ndat$rt, 2 * A, 2 * x["b_pos"], t0, 2 * x[ndat$v_pos], sqrt(2))) # nolint
  c(yes, no)
}

ll_CYEX <- function(x, data) { # nolint
  ydat <- data[data$response == 2, ]
  ndat <- data[data$response == 1, ]
  A <- x["A"] # nolint
  t0 <- x["t0"]

  yes <- dlba_norm(ydat$rt, 2 * A, 2 * x["b_pos"], t0, 2 * x[ydat$v_pos], sqrt(2)) * # nolint
    (1 - plba_norm(ydat$rt, A, x["b_neg"], t0, x[ydat$v_neg], 1))**2
  no <- 2 *
    dlba_norm(ndat$rt, A, x["b_neg"], t0, x[ndat$v_neg], 1) *
    (1 - plba_norm(ndat$rt, A, x["b_neg"], t0, x[ndat$v_neg], 1)) *
    (1 - plba_norm(ndat$rt, 2 * A, 2 * x["b_pos"], t0, 2 * x[ndat$v_pos], sqrt(2))) # nolint
  c(yes, no)
}

ll_CNST <- function(x, data) { # nolint
  ydat <- data[data$response == 2, ]
  ndat <- data[data$response == 1, ]
  A <- x["A"] # nolint
  t0 <- x["t0"]

  yes <- 2 *
    dlba_norm(ydat$rt, A, x["b_pos"], t0, x[ydat$v_pos], 1) *
    (1 - plba_norm(ydat$rt, A, x["b_pos"], t0, x[ydat$v_pos], 1)) *
    (1 - plba_norm(ydat$rt, 2 * A, 2 * x["b_neg"], t0, 2 * x[ydat$v_neg], sqrt(2))) # nolint
  no <- dlba_norm(ndat$rt, 2 * A, 2 * x["b_neg"], t0, 2 * x[ndat$v_neg], sqrt(2)) * # nolint
    (1 - plba_norm(ndat$rt, A, x["b_pos"], t0, x[ndat$v_pos], 1))**2
  c(yes, no)
}

ll_CNEX <- function(x, data) { # nolint
  ydat <- data[data$response == 2, ]
  ndat <- data[data$response == 1, ]
  A <- x["A"] # nolint
  t0 <- x["t0"]

  yes <- 2 *
    dlba_norm(ydat$rt, A, x["b_pos"], t0, x[ydat$v_pos], 1) *
    plba_norm(ydat$rt, A, x["b_pos"], t0, x[ydat$v_pos], 1) *
    (1 - plba_norm(ydat$rt, 2 * A, 2 * x["b_neg"], t0, 2 * x[ydat$v_neg], sqrt(2))) # nolint
  no <- dlba_norm(ndat$rt, 2 * A, 2 * x["b_neg"], t0, 2 * x[ndat$v_neg], sqrt(2)) * # nolint
    (1 - plba_norm(ndat$rt, A, x["b_pos"], t0, x[ndat$v_pos], 1)**2)
  c(yes, no)
}

ll_CB <- function(x, data) { # nolint
  ydat <- data[data$response == 2, ]
  ndat <- data[data$response == 1, ]
  A <- x["A"] # nolint
  t0 <- x["t0"]

  yes <- dlba_norm(ydat$rt, 2 * A, 2 * x["b_pos"], t0, 2 * x[ydat$v_pos], sqrt(2)) * # nolint
    (1 - plba_norm(ydat$rt, 2 * A, 2 * x["b_neg"], t0, 2 * x[ydat$v_neg], sqrt(2))) # nolint
  no <- dlba_norm(ndat$rt, 2 * A, 2 * x["b_neg"], t0, 2 * x[ndat$v_neg], sqrt(2)) * # nolint
    (1 - plba_norm(ndat$rt, 2 * A, 2 * x["b_pos"], t0, 2 * x[ndat$v_pos], sqrt(2))) # nolint
  c(yes, no)
}

ll_funcs <- c(ll_IST, ll_IEX, ll_CYST, ll_CYEX, ll_CNST, ll_CNEX, ll_CB)

# Specify the log likelihood function -----------------------------------------
dirichlet_mix_ll <- function(x, data) {
  x <- exp(x)

  # Enforce alphas to be greater than 0.01 and less than 100
  if (any(x[mix_counts] < 0.01) || any(x[mix_counts] > 100)) {
    return(-1e10)
  }

  # Enforces b cannot be less than A, b parameter is thus threshold - A
  x["b_pos"] <- x["b_pos"] + x["A"]
  x["b_neg"] <- x["b_neg"] + x["A"]

  # all decision rules
  rdev <- rdirichlet(1, x[mix_counts])
  func_idx <- sample(mix_counts, 1, prob = rdev)
  ll_func <- ll_funcs[[func_idx]]
  trial_ll <- ll_func(x, data)
  new_like <- (1 - p_contam) * trial_ll +
    p_contam * (dunif(data$rt, min_rt, max_rt) / 2)
  sum(log(pmax(new_like, 1e-10)))
}

test.likelihood <- function(x, num.values = 9, fake.data, dimensions, ll_func, server = FALSE, ...) {
  pars <- names(x)
  x.tmp <- x
  sequence <- seq(from = -.22, to = .22, length.out = num.values)
  op <- par(mfrow = dimensions)
  par_likelihoods <- lapply(setNames(pars, pars), function(p) {
    testvalues <- sequence + true.x[[p]]
    tmp <- unlist(lapply(testvalues, function(i) { # for each test vlaue, it will apply the function where i = testvalue.
      x.tmp[[p]] <- i
      ll_func(x = x.tmp, data = fake.data)
    }))

    plot(
      x = testvalues,
      y = tmp,
      type = "b",
      main = p,
      xlab = "log par values",
      ylab = "loglikelihood",
      ...
    )
    abline(v = true.x[[p]], col = "red")
    return(tmp)
  })
  par(op)
  return(par_likelihoods)
}

# n.posterior <- 20 # Number of samples from posterior distribution for each parameter.
# pp.data <- list()
# S <- unique(wgnmks2008$subject)
# data <- split(x = wgnmks2008, f = wgnmks2008$subject)
# for (s in S) {
#   cat(s, " ")
#   iterations <- round(seq(from = 1051, to = sampled$samples$idx, length.out = n.posterior))
#   for (i in 1:length(iterations)) {
#     x <- sampled$samples$alpha[, s, iterations[i]]
#     names(x) <- pars
#     tmp <- SDT_ll_fast(x = x, data = wgnmks2008[wgnmks2008$subject == s, ], sample = TRUE)
#     if (i == 1) {
#       pp.data[[s]] <- cbind(i, tmp)
#     } else {
#       pp.data[[s]] <- rbind(pp.data[[s]], cbind(i, tmp))
#     }
#   }
# }
