library(tidyverse)
library(mqfa)

load("data/illustration-q1-out-stan.RData")
load("data/illustration-q1-data.RData")

uncertainty_fun <- function(values_seq, fun, out_mcmc, R_selected = 1000){
  R <- nrow(out_mcmc)
  if(R_selected > R){
    stop("R_selected must be less than R")
  }
  uncertainty_matrix <- sapply(values_seq, \(x) apply(out_mcmc[seq(1, R, len = R_selected), ], MARGIN = 1, \(theta) fun(x, theta)))
  pm <- colMeans(uncertainty_matrix)
  q5 <- apply(uncertainty_matrix, 2, quantile, probs = 0.05)
  q95 <- apply(uncertainty_matrix, 2, quantile, probs = 0.95)
  return(data.frame("values_seq" = values_seq, "pm" = pm, "q5" = q5, "q95" = q95))
}


# these work for Gibbs sampler output
# colnames(out$theta) <- paste0("theta[", rep(1:4, times = q*K),",", rep(1:q, 4*K),",", rep(1:K, each = 4*q),"]")
# idx <- map(1:(q*K), \(i) i:(i+3)+(i-1)*3)

perm <- c(1, 3, 2)
idx <- map(1:(K), \(k) seq(k, k + 9, by = 3))
theta_k <- map(idx, \(x) out_stan$theta[, x])
theta_k <- theta_k[perm]
# theta_k[[3]] <- t(apply(theta_k[[3]], 1, ifa2::theta_minus))
theta_k[[2]] <- t(apply(theta_k[[2]], 1, ifa2::theta_minus))

x_seq <- seq(-4, 4, len = 200)
fx_pm <- imap(theta_k, \(x, i) uncertainty_fun(values_seq = x_seq,
                                               fun = dfgld,
                                               out_mcmc = x,
                                               R_selected = 500) |>
                as_tibble() |>
                mutate(component = i))
fx_pm <- do.call(what = bind_rows, args = fx_pm)

fx_true <- imap(array_branch(data$params$theta, 3),
                \(x, i) tibble("fx" = dfgld(x = x_seq, theta = c(x)),
                               "component" = i))
fx_true <- map(fx_true, \(x) bind_cols(x, "x" = x_seq)) |>
  do.call(args = _, what = bind_rows)

ggplot(fx_pm, aes(x = values_seq)) +
  geom_line(aes(y = pm), alpha = 0.8, lty = 3) +
  geom_ribbon(aes(ymin = q5, ymax = q95), alpha = 0.2) +
  geom_line(data = fx_true, aes(x = x, y = fx), col = "royalblue", lty = 1)+
  facet_wrap(facets = vars(component), labeller = label_both) +
  labs(x = "x", y = "f(x)") +
  theme_bw()+
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5),
        strip.background = element_rect(fill = NA))

ggsave(filename = "output/illustration-mfa-fgld-q1.pdf", width = 6, height = 2.5)
