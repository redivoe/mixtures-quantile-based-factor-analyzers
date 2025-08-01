library(purrr)
library(furrr)
library(cmdstanr)

library(mqfa)

load("data/datasets.RData")
compiled <- cmdstan_model("../stan/mifa-gk.stan")

reps <- 10
plan(multisession(workers = 16))

out_mfa_gk <- vector(mode = "list", length = length(datasets)) |>
  set_names(names(datasets))

fit_grid <- function(x){
  fit <- compiled$variational(data = list("X" = datasets[[i]]$X,
                                          "n" = nrow(datasets[[i]]$X),
                                          "p" = ncol(datasets[[i]]$X),
                                          "q" = x$q,
                                          "K" = datasets[[i]]$K))
  
  start_time <- Sys.time()
  out <- offline_z_mifa_cmdstanr(fit = fit,
                                 X = datasets[[i]]$X,
                                 distr = "gk",
                                 monte_carlo = TRUE,
                                 B = 200)
  end_time <- Sys.time()
  
  out$time1 <- fit$time()$total
  out$time2 <- end_time - start_time
  
  return(out)
}

for(i in 1:length(datasets)){
  grid <- cross(.l = list("rep" = 1:reps, "q" = 1:datasets[[i]]$q_max))
  out_grid <- future_map(.x = grid,
                         .f = insistently(fit_grid, rate = rate_delay(pause = 0, max_times = 10)),
                         .options = furrr_options(seed = 123))
  
  grid_df <- grid |> transpose() |> map(list_simplify) |> as.data.frame()
  
  grid_df$time1 <- map_dbl(out_grid, "time1")
  grid_df$time2 <- map_dbl(out_grid, "time2")
  grid_df$ari <- map_dbl(out_grid, \(x) mclust::adjustedRandIndex(datasets[[i]]$z, x$z_map))
  grid_df$mr <- map_dbl(out_grid, \(x) misc(datasets[[i]]$z, x$z_map))
  grid_df$waic <- map_dbl(out_grid, "waic")
  grid_df$wbic <- map_dbl(out_grid, "wbic")
  
  out_mfa_gk[[i]] <- grid_df
  cat("\n\n", i, "\n\n")
}

save(list = c("out_mfa_gk"), file = "output/out-mfa-gk-stan.RData")


