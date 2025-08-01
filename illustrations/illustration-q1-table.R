library(tidyverse)
library(gt)

load("data/illustration-q1-out-gibbs.RData")
load("data/illustration-q1-out-stan.RData")

summary_tab_gt <-  bind_rows(
  with(out$summary, tibble("param" = names(ess), ess, bias)) |> 
  mutate(method = "gibbs"),
  with(out_stan$summary, tibble("param" = names(ess), ess, bias)) |> 
  mutate(method = "stan")) |> 
  pivot_wider(names_from = method, values_from = c(ess, bias)) |> 
  mutate(ess_time_relative = ess_stan / out_stan$summary$time / (ess_gibbs / out$summary$time)) |> 
  gt() |>
  tab_spanner(label = "ESS", columns = c(2, 3)) |> 
  tab_spanner(label = "|Bias| mean", columns = c(4, 5)) |> 
  cols_label(param ~ "Parameter",
             ess_gibbs ~ "Gibbs",
             ess_stan ~ "Stan",
             bias_gibbs ~ "Gibbs",
             bias_stan ~ "Stan",
             ess_time_relative ~ "Relative ESS/time") |>
  fmt_number(columns = c(ess_gibbs, ess_stan), decimals = 0) |>
  fmt_number(columns = c(bias_gibbs, bias_stan), decimals = 3) |> 
  fmt_number(columns = c(ess_time_relative), decimals = 2)

summary_tab_gt |>
  as_latex() |>
  as.character() |>
  writeLines()




# loadings
# colnames(out$L) <- paste0("L[", rep(1:p, times = q*K),",", rep(1:q, p*K),",", rep(1:K, each = p*q),"]")
# 
# L_sum <- out$L |>
#   tidy_draws() |>
#   spread_draws(L[var, factor, component]) |>
#   summarise(mean = mean(L),
#             q10 = quantile(L, 0.025),
#             q90 = quantile(L, 0.975),
#             .groups = "drop")
# 
# L_sum |>
#   ggplot(aes(y = var, x = mean, xmin = q10, xmax = q90, col = factor(component)))+
#   geom_point(alpha = 0.9)+
#   geom_path()+
#   geom_errorbarh(height = 0, alpha = 0.9)+
#   facet_wrap(facets = vars(component), labeller = label_both)+
#   scale_y_continuous(breaks = 1:p, minor_breaks = NULL)+
#   labs(y = "variable", x = "factor loading", col = "component")+
#   theme_bw()+
#   theme(legend.position = "none")
