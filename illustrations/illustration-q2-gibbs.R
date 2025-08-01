library(coda)
library(purrr)
library(mqfa)

load("data/illustration-q2-data.RData")

R <- 5e3
burn_in <- 5e3

set.seed(1)
out_mifa_fgld <- mifa_fgld(X = data$X, K, q, R = R, burn_in = burn_in, z_update = "fgld_marginal")
table(data$params$z, out_mifa_fgld$z_map)

save(list = "out_mifa_fgld", file = "data/illustration-q2-out-gibbs.RData")

out_mifa_fgld$theta |> 
  as.mcmc() |> 
  plot()
