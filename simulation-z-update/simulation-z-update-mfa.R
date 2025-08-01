library(purrr)
library(furrr)
library(mqfa)

settings <- list(n = c(100, 200, 500),
                 p = c(10, 20, 50),
                 q = c(1, 2, 3),
                 K = c(2, 3, 4)) |> 
  cross()

n_inits <- 20
R <- 5e3
burn_in <- 5e3


set.seed(1)
data <- map(settings, \(x) r_mfa_data(n = x$n, p = x$p, q = x$q, K = x$K, sigma0 = 5))
inits <- map(settings,
             \(x) map(1:n_inits, \(i) mfa_random_init(n = x$n, p = x$p, q = x$q, K = x$K)))

names <- c("marginal", "conditional")

run_setting <- function(i) {
  out <- vector(mode = "list", length = length(inits[[i]]))
  for (j in 1:length(inits[[i]])) {
    out[[j]] <- map(
      names,
      \(x) mfa(
        X = data[[i]]$X,
        K = settings[[i]]$K,
        q = settings[[i]]$q,
        R = R,
        burn_in = burn_in,
        params_init = inits[[i]][[j]],
        z_update = x
      )$z_map
    ) |>
      set_names(nm = names)
  }
  out <- map_depth(
    .x = out,
    .depth = 2,
    \(x) mclust::adjustedRandIndex(data[[i]]$params$z, x)
  ) |>
    transpose() |>
    map(list_simplify)
  return(out)
}

plan(multisession(workers = 16))
out <- future_map(.x = 1:length(settings),
                  .f = run_setting,
                  .options = furrr_options(seed = 123, chunk_size = 1))

save.image("simulation-z-update-mfa.RData")




