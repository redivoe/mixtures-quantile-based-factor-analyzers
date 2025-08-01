library(tidyverse)
library(scales)
library(ggdist)

load("simulation-z-update-mfa.RData")

settings_df <- settings |>
  transpose() |> 
  map(list_simplify) |> 
  as_tibble() |> 
  mutate(marginal = map(out, "marginal"),
         conditional = map(out, "conditional"),
         diff = map2(marginal, conditional, \(x, y) x - y)) |> 
  unnest(diff) |> 
  select(-c(marginal, conditional)) |>
  mutate(across(-diff, as_factor))

settings_df |> 
  ggplot(aes(x = diff, y = after_stat(density)))+
  geom_vline(xintercept = 0, lty = 3)+
  geom_histogram(bins = 30, alpha = 0.6)+
  coord_cartesian(xlim = c(-1, 1))+
  labs(x = "difference in ARI")+
  theme_bw()
ggsave(filename = "diff-ari.pdf", width = 4, height = 3)

mean(settings_df$diff == 0)
mean(settings_df$diff > 0)


settings_df |> 
  ggplot(aes(x = diff, y = after_stat(density)))+
  facet_grid(rows = vars(n), cols = vars(p))+
  geom_vline(xintercept = 0, lty = 3)+
  geom_histogram(bins = 30, alpha = 0.6)+
  coord_cartesian(xlim = c(-1, 1))+
  labs(x = "difference in ARI")+
  theme_bw()


