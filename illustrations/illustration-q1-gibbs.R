library(purrr)
library(coda)
library(mqfa)

load("data/illustration-q1-data.RData")

R <- 5e3
burn_in <- 5e3

start_time <- Sys.time()
out <- mifa_fgld(X = data$X, K = K, q = q, R = R, burn_in = burn_in,
                 params_init = params_init)
end_time <- Sys.time()
elapsed_time <- difftime(end_time, start_time, units = "secs") |>
  as.numeric()

# classification
table(out$z_map, data$params$z)

# out$mu |> as.mcmc() |> plot()
# out$theta |> as.mcmc() |> plot()
# out$prop |> as.mcmc() |> plot()

# indices of parameters excluded from summary stats
idx_excl_params <- map_dbl(c("Y", "omega", "z", "z_map", "waic"), \(x) which(names(out) == x))

# bias
# posterior means
pm <- map(out[-idx_excl_params], colMeans)

perm_true <- ifa2::best_col_permutation(X = matrix(pm$mu, nrow = p, ncol = K),
                                            Xhat = data$params$mu)$order
true_params <- data$params[names(pm)]
true_params[-5] <- map(true_params[-5], \(x) ifa2::lastInd(x, perm_true))
true_params$prop <- true_params$prop[perm_true]

plot(true_params$L, pm$L)
drop(true_params$L)
matrix(pm$L, ncol = 3) |> round(1)
# sign switched for L1
pm_L <- matrix(pm$L, ncol = K)
pm_L[,c(1)] <- -pm_L[,c(1)]
pm_L[,c(3)] <- -pm_L[,c(3)]
pm$L <- pm_L

pm_theta <- matrix(pm$theta, ncol = K)
pm_theta[, c(1)] <- apply(pm_theta[, c(1), drop = FALSE], 2, theta_minus)
pm_theta[, c(3)] <- apply(pm_theta[, c(3), drop = FALSE], 2, theta_minus)
pm$theta <- pm_theta

# bias
true_params <- map(true_params, c)
bias <- map2(pm, true_params, \(x, y) abs(x - y)) |>
  map_dbl(mean)

# ESS
ess <- map(out[-idx_excl_params], effectiveSize) |>
  map_dbl(mean)

out$summary <- list("ess" = ess,
                    "bias" = bias,
                    "time" = elapsed_time)

save(list = "out", file = "data/illustration-q1-out-gibbs.RData")

