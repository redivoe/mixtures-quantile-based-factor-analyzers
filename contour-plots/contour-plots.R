library(mqfa)
library(tidyverse)
library(furrr)

mu <- list(rep(0, 2))
L <- list(matrix(c(1, 1), nrow = 2, ncol = 1))
delta <- list(0.5, 0.9)
kappa <- list(0.1, 10)
psi <- list(0.05, 0.2)


params <- cross(list("mu" = mu, "L" = L, "delta" = delta, "kappa" = kappa, "psi" = psi))
theta <- map(params, \(x) matrix(complete_theta_std(x$delta, x$kappa), nrow = 4, ncol = 1))
params <- transpose(params)
params$theta <- theta
params <- transpose(params)

B <- 100
lim <- 2.5
points <- expand.grid("x1" = seq(-lim, lim, len = B),
                      "x2" = seq(-lim, lim, len = B))

plan(multisession(workers = 4))

# Monte Carlo
# fx <- future_map(params,
#                  \(x) fx_fgld_mc(x = points, mu = x$mu, L = x$L, psi = x$psi, theta = x$theta, B = 5000),
#                  .options = furrr_options(seed = 1))

# Cubature
fx <- future_map(params,
                 \(x) fx_fgld_cubature(x = points, mu = x$mu, L = x$L, psi = x$psi, theta = x$theta),
                 .options = furrr_options(seed = 1))

plottable <- transpose(params)[c("delta", "kappa", "psi")] |>
  map(unlist) |>
  as_tibble()
plottable$fx <- map(fx, \(x) bind_cols(points, "fx" = x))
plottable <- unnest(plottable, fx)

(p_contour <- ggplot(plottable, aes(x1, x2, z = fx))+
  geom_contour_filled(alpha = 0.7, bins = 15)+
  facet_grid(rows = vars(kappa), cols = vars(psi, delta),
             labeller = label_bquote(cols = list(delta==.(delta), psi==.(psi)),
                                     rows = kappa==.(kappa))
             )+
  coord_equal()+
  scale_fill_viridis_d(option = "G", direction = -1)+
  labs(x = NULL, y = NULL, col = NULL)+
  guides(fill = "none")+
  theme_bw()+
  theme(panel.grid = element_blank()))

ggsave(filename = "contour-plots.pdf", width = 7.5, height = 4)
