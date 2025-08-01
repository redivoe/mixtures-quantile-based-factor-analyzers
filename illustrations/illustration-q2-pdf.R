library(tidyverse)
library(patchwork)

dmfgld <- function(x, theta){
  if(is.null(dim(x))){
    p <- length(x)
    out <- 0
    for(j in 1:p){
      out <- out + ifa2::dfgld(x[j], theta[, j], log = TRUE)
    }
    return(exp(out))
  }else{
    p <- ncol(x)
    n <- nrow(x)
    out <- numeric(length = n)
    for(j in 1:p){
      out <- out + ifa2::dfgld(x[, j], theta[, j], log = TRUE)
    }
    return(exp(out))
  }
}


load("data/illustration-q2-data.RData")
# load("data/illustration-q2-out-gibbs.RData")
load("data/illustration-q2-out-stan.RData")


ray <- 3
grid <- expand.grid("y1" = seq(-ray, ray, len = 200),
                    "y2" = seq(-ray, ray, len = 200))

fx_true <- tibble("theta" = c(data$params$theta),
                  "component" = rep(1:K, each = 4 * q)) |>
  group_by(component) |>
  reframe(theta = list(matrix(theta, nrow = 4))) |>
  ungroup() |>
  mutate(grid = list(grid)) |>
  rowwise() |>
  mutate(fx = list(dmfgld(grid, theta))) |>
  ungroup() |>
  select(-theta) |>
  unnest(c(grid, fx))

# fitted distribution
theta_pm <- colMeans(do.call(cbind, out_stan$theta_rsp))
# theta_pm <- colMeans(out_mifa_fgld$theta)


fx_pm <- tibble("theta" = theta_pm,
                "component" = rep(1:K, each = 4 * q)) |>
  group_by(component) |>
  reframe(theta = list(matrix(theta, nrow = 4))) |>
  ungroup() |>
  mutate(grid = list(grid)) |>
  rowwise() |>
  mutate(fx = list(dmfgld(grid, theta))) |>
  ungroup() |>
  select(-theta) |>
  unnest(c(grid, fx))

p_fx_pm <- ggplot()+
  geom_contour_filled(data = fx_pm,
                      aes(x = y1, y = y2, z = fx),
                      alpha = 0.7,
                      bins = 7)+
  facet_wrap(facets = vars(component), labeller = label_both)+
  coord_equal()+
  scale_fill_viridis_d(option = "G", direction = -1)+
  labs(x = NULL, y = NULL, col = NULL, title = "Densities at the posterior mean")+
  guides(fill = "none")+
  theme_bw()

p_fx_true <- ggplot()+
  geom_contour_filled(data = fx_true,
                      aes(x = y1, y = y2,z = fx),
                      alpha = 0.7,
                      bins = 7)+
  facet_wrap(facets = vars(component), labeller = label_both)+
  coord_equal()+
  scale_fill_viridis_d(option = "G", direction = -1)+
  labs(x = NULL, y = NULL, col = NULL, title = "True densities")+
  guides(fill = "none")+
  theme_bw()

(p_fx_true / p_fx_pm) & theme(panel.grid = element_blank(),
                              plot.title = element_text(hjust = 0.5))
ggsave("output/illustration-mfa-fgld-q2-stan.pdf", width = 4, height = 5)



##-----------------------
##  Latent factors plot  
##-----------------------

# # latent factors
# Y_pm <- MCMCpstr(object = out_mifa_fgld, params = "Y")$Y |>
#   `colnames<-`(c("y1", "y2"))
# 
# scores_pm <- cbind(Y_pm, "z" = z_map) |>
#   as_tibble()
# 
# scores_true <- data$params$Y |>
#   `colnames<-`(c("y1", "y2")) |>
#   as_tibble() |>
#   mutate(z = data$params$z)
# 
# # offset <- 3
# ggplot()+
#   geom_contour_filled(data = fx_true,
#                       aes(x = y1, y = y2, z = fx),
#                       alpha = 0.7,
#                       bins = 6)+
#   geom_point(data = scores_pm,
#              aes(x = y1, y = y2),
#              pch = 1,
#              alpha = 0.5)+
#   # geom_rug(alpha = 0.1)+
#   facet_wrap(facets = vars(z), labeller = label_both)+
#   # geom_segment(data = loadings,
#   #              aes(x = -offset, y = -offset, xend = y1 - offset, yend = y2 - offset, col = vars),
#   #              arrow = arrow(length = unit(0.5, "lines")))+
#   coord_equal()+
#   scale_fill_viridis_d(option = "G", direction = -1)+
#   labs(x = NULL, y = NULL, col = NULL)+
#   guides(fill = "none")+
#   theme_bw()+
#   theme(panel.grid = element_blank())
