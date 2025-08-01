library(cmdstanr)
library(posterior)
library(tidybayes)
library(tidyverse)
library(furrr)
library(mqfa)
library(patchwork)


n <- 50
K <- 1
q <- 1
p <- 10

set.seed(1)
mu <- array(2 * rt(n = p * K, df = 2), dim = c(p, K))
L <- array(rt(n = p * q * K, df = 2), dim = c(p, q, K))
theta <- array(rtheta_std(p = q*K), dim = c(4, q, K))
psi <- array(rlnorm(n = p * K), dim = c(p, K))
prop <- c(1)

params <- list("mu" = mu,
               "L" = L,
               "theta" = theta,
               "psi" = psi,
               "prop" = prop,
               "common_theta" = FALSE,
               "common_psi" = FALSE)

set.seed(1)
data <- r_mifa_fgld_data(n = n, p = p, q = q, K = K, params = params)


R <- 2e3
burn_in <- 2e3

priors <- c("Standard", "Cauchy", "Diffuse")

m1 <- cmdstan_model("stan/mifa-fgld.stan")
m2 <- cmdstan_model("stan/mifa-fgld-cauchy-priors.stan")
m3 <- cmdstan_model("stan/mifa-fgld-diffuse-priors.stan")
models <- list(m1, m2, m3) |> 
  set_names(priors)

data_list <- list("X" = data$X,
                  "n" = n,
                  "p" = p,
                  "q" = q,
                  "K" = K)

plan(multisession(workers = 3))
out <- future_map(models, \(x) x$sample(data = data_list,
                                                  chains = 1,
                                                  iter_warmup = R,
                                                  iter_sampling = burn_in),
                  .options = furrr_options(seed = 1))

plot_ci <- function(param){
  imap(out, \(x, i) x$draws(out) |>
         spread_draws(psi[component, var]) |> 
         mutate(prior = i)) |> 
    do.call(args = _, what = bind_rows) |> 
    ggplot(aes(y = .data[["psi"]], x = .data[["var"]], col = prior)) +
    stat_pointinterval(.width = c(0.9, 0.99), position = position_dodge(width = 0.5)) +
    theme_bw()
}

p1 <- imap(out, \(x, i) x$draws() |>
       spread_draws(mu[component, var], Lambda[component, var, latent_var], psi[component, var]) |> 
       mutate(Lambda = abs(Lambda)) |> 
       pivot_longer(cols = c(mu, Lambda, psi), names_to = "param", values_to = "value") |> 
       mutate(prior = i)) |> 
  do.call(args = _, what = bind_rows) |> 
  ggplot(aes(y = value, x = factor(var), col = prior)) +
  facet_grid(rows = vars(param), scales = "free",
             labeller = label_parsed)+
  stat_pointinterval(.width = c(0.9, 0.99), position = position_dodge(width = 0.5)) +
  labs(x = NULL, y = NULL, col = "Prior settings")+
  theme_bw()+
  theme(strip.text.y = element_text(angle = 0))


p2 <- imap(out, \(x, i) x$draws() |>
       spread_draws(delta[component, latent_var], kappa[component, latent_var]) |> 
       pivot_longer(cols = c(delta, kappa), names_to = "param", values_to = "value") |> 
       mutate(prior = i)) |> 
  do.call(args = _, what = bind_rows) |> 
  ggplot(aes(y = value, x = factor(latent_var), col = prior)) +
  facet_grid(rows = vars(param), scales = "free",
             labeller = label_parsed)+
  stat_pointinterval(.width = c(0.8, 0.9), position = position_dodge(width = 0.5)) +
  labs(x = NULL, y = NULL, col = "Prior settings")+
  theme_bw()+
  theme(strip.text.y = element_text(angle = 0))

p1 / (plot_spacer() + plot_spacer() + p2) + plot_layout(heights = c(1, 2/3), guides = "collect")
ggsave("prior-sensitivity.pdf", width = 6, height = 7)
  