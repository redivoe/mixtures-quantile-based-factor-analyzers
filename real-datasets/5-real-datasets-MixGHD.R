library(purrr)
library(furrr)
library(MixGHD)

load("data/datasets.RData")

plan(multisession(workers = 4))

out_mghfa <- list()
for(i in 1:length(datasets)){
  out_mghfa[[i]] <- future_map(1:datasets[[i]]$q_max,
                               \(q) insistently(f = MGHFA,
                                                rate = rate_delay(pause = 0, max_times = 10))(data = datasets[[i]]$X,
                                                               G = datasets[[i]]$K,
                                                               q = q),
                               .options = furrr_options(seed = 2))
  cat("\n\n", i)
}
out_mghfa <- set_names(out_mghfa, names(datasets))

z_map <- map_depth(out_mghfa, 2, \(x) x@map)
bic <- map_depth(out_mghfa, 2, \(x) x@BIC) |> 
  map(list_simplify)

out_mghfa <- imap(datasets, \(x, i) map(z_map[[i]], \(y) c("ari" = mclust::adjustedRandIndex(x = y, x$z),
                                              "mr" = mqfa::misc(classification = y, truth = x$z)))) |> 
  map(bind_rows)

out_mghfa <- map2(out_mghfa, datasets, \(x, y) mutate(x, q = 1:y$q_max))
out_mghfa <- map2(out_mghfa, bic, \(x, y) mutate(x, bic = y))

save(list = c("out_mghfa"), file = "output/out-mfa-skew.RData")


