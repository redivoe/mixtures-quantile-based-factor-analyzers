library(cmdstanr)
library(purrr)
library(bayesplot)
library(ggplot2)
library(factor.switching)
library(coda)

library(mqfa)

load("data/illustration-q1-data.RData")

R <- 5e3
burn_in <- 5e3

compiled <- cmdstan_model("../stan/mifa-fgld.stan")

data_list <- list("X" = data$X,
                  "n" = n,
                  "p" = p,
                  "q" = q,
                  "K" = K)

params_init_stan <- list()
params_init_stan$psi <- array_branch(params_init$psi, 2)
params_init_stan$Lambda <- array_branch(params_init$L, 3)
params_init_stan$mu <- array_branch(params_init$mu, 2)
params_init_stan$prop <- params_init$prop
params_init_stan$delta <- as.list(params_init$theta[3,,])
params_init_stan$kappa <- as.list(params_init$theta[4,,])

start_time <- Sys.time()
fit <- compiled$sample(data = data_list,
                       chains = 1,
                       parallel_chains = 5,
                       iter_warmup = R,
                       iter_sampling = burn_in,
                       init = list(params_init_stan))
end_time <- Sys.time()
elapsed_time <- difftime(end_time, start_time, units = "secs") |>
  as.numeric()


z_mcmc <- fit$draws("z") |> unclass() |> drop()

out_stan <- list()
out_stan$L <- fit$draws(c("Lambda")) |> 
  unclass() |> 
  drop()
out_stan$theta <- fit$draws(c("alpha", "beta", "delta", "kappa")) |> 
  unclass() |> 
  drop()
out_stan$psi <- fit$draws(c("psi")) |> 
  unclass() |> 
  drop()
out_stan$mu <- fit$draws(c("mu")) |> 
  unclass() |> 
  drop()
out_stan$prop <- fit$draws(c("prop")) |> 
  unclass() |> 
  drop()

mcmc_trace(fit$draws(c("delta", "kappa")))+
  theme_bw()+
  theme(legend.position = "none")
ggsave(filename = "output/illustration-q1-stan-traceplot-delta-kappa.pdf",
       width = 8, height = 4)

mcmc_trace(out_stan$L)+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_blank())
ggsave(filename = "output/illustration-q1-stan-traceplot-lambda.pdf",
       width = 10, height = 6)

ess <- map(out_stan, effectiveSize) |> 
  map_dbl(mean)

pm <- map(out_stan, colMeans)

true_params <- data$params[names(pm)]
matrix(pm$mu, nrow = K) |> round(1)
true_params$mu
perm_true <- c(1, 3, 2)
true_params[-5] <- map(true_params[-5], \(x) ifa2::lastInd(x, perm_true))
true_params$prop <- true_params$prop[perm_true]

true_params$psi <- t(true_params$psi)
true_params$mu <- t(true_params$mu)
true_params$theta <- aperm(true_params$theta, c(3, 1, 2)) |> drop()
true_params$L <- aperm(true_params$L, c(3, 1, 2)) |> drop()

true_params$L
pm$L |> matrix(nrow = K, ncol = p) |> round(1)
true_params$L[c(1, 3), ] <- - true_params$L[c(1, 3), ] 

pm_theta <- matrix(pm$theta, nrow = K)
pm_theta[3, ] <- theta_minus(pm_theta[3, ])
pm_theta[1, ] <- theta_minus(pm_theta[1, ])
pm$theta <- pm_theta

plot(pm$L, true_params$L, asp = 1)
abline(0, 1)
plot(pm$theta, true_params$theta, asp = 1)
abline(0, 1)

bias <- map2(pm, true_params, \(x, y) abs(x - y)) |>
  map_dbl(mean)

out_stan$summary <- list("ess" = ess,
                         "bias" = bias,
                         "time" = elapsed_time)


##---------------------------------
##  Checking for factor switching  
##---------------------------------

# rsp_names <- paste0("Lambda", paste0(rep(paste0("V", 1:p), times = q), "_", rep(1:q, each = p)))
# idx_k <- map(1:K, \(k) grep(paste0("^Lambda\\[", k, ","), x = colnames(out_stan$L)))
# out_stan$L_sep <- map(idx_k, \(x) out_stan$L[, x] |> `colnames<-`(rsp_names))
# out_rsp <- map(out_stan$L_sep, \(x) rsp_exact(lambda_mcmc = x, rotate = FALSE))

save(list = c("out_stan"),
     file = "data/illustration-q1-out-stan.RData")
