library(purrr)
library(furrr)

library(mqfa)

timed_fun <- function(foo){
  function(...){
    start_time <- Sys.time()
    out <- foo(...)
    end_time <- Sys.time()
    elapsed_time <- difftime(end_time, start_time, units = "secs") |>
      as.numeric()
    return(list("out" = out,
                "time" = elapsed_time))
  }
}


load("data/datasets.RData")

reps <- 10
R <- 2e3
burn_in <- 2e3
plan(multisession(workers = 16))

out_mfa <- vector(mode = "list", length = length(datasets)) |>
  set_names(names(datasets))

for(i in 1:length(datasets)){
  grid <- cross(.l = list("rep" = 1:reps, "q" = 1:datasets[[i]]$q_max))
  out_grid <- future_map(grid,
                         \(x) timed_fun(mfa)(X = datasets[[i]]$X,
                                             K = datasets[[i]]$K,
                                             q = x$q,
                                             R = R,
                                             burn_in = burn_in,
                                             z_update = "marginal"),
                          .options = furrr_options(seed = 1))

  grid_df <- grid |> transpose() |> map(list_simplify) |> as.data.frame()
  grid_df$time <- map_dbl(out_grid, "time")
  out_grid <- map(out_grid, "out")
  grid_df$ari <- map_dbl(out_grid, \(x) mclust::adjustedRandIndex(datasets[[i]]$z, x$z_map))
  grid_df$mr <- map_dbl(out_grid, \(x) misc(datasets[[i]]$z, x$z_map))
  grid_df$waic <- map_dbl(out_grid, "waic")
  grid_df$wbic <- map_dbl(out_grid, "wbic")

  out_mfa[[i]] <- grid_df
}

save(list = c("out_mfa"), file = "output/out-mfa.RData")

