library(tidyverse)
library(gt)

load("data/datasets.RData")
load("output/out-mfa.RData")
load("output/out-mfa-fgld-stan.RData")
load("output/out-mfa-gk-stan.RData")
load("output/out-mfa-skew.RData")

out <- map(list("MFA" = out_mfa,
         "MFA fgld" = out_mfa_fgld,
         "MFA gk" = out_mfa_gk, 
         "MGHFA" = out_mghfa),
    \(out) imap(out, \(x, i) mutate(x, "dataset" = i)) |> 
      do.call(what = bind_rows, args = _) |> 
      as_tibble()) |> 
  imap(\(x, i) mutate(x, "model" = i)) |> 
  do.call(args = _, what = bind_rows) |> 
  mutate(across(c(model, q), as_factor))

out |> 
  group_by(model, dataset, q) |> 
  slice_min(wbic) |> 
  select(model, q, dataset, ari) |> 
  pivot_wider(values_from = c(ari), names_from = dataset) |> 
  ungroup() |>
  arrange(model, q) |> 
  gt(groupname_col = "model") |> 
  fmt_number(decimals = 2) |> 
  sub_missing() |> 
  tab_options(row_group.as_column = TRUE) |> 
  as_latex() |>
  as.character() |>
  writeLines()

# q selection for each model (italics in Table 3)
out |> 
  filter(model != "MGHFA") |> 
  group_by(model, dataset) |> 
  slice_min(wbic) |> 
  select(dataset, model, q) |> 
  pivot_wider(names_from = dataset, values_from = q)
out |> 
  filter(model == "MGHFA") |> 
  group_by(model, dataset) |> 
  slice_min(bic) |> 
  select(dataset, model, q) |> 
  pivot_wider(names_from = dataset, values_from = q)

out |> 
  filter(!is.na(time1 & time2)) |> 
  group_by(dataset, q) |> 
  summarise(t = mean(time1 + time2), .groups = "drop") |> 
  pivot_wider(values_from = c(t), names_from = dataset) |> 
  gt() |> 
  fmt_number(decimals = 1) |> 
  sub_missing() |> 
  as_latex() |>
  as.character() |>
  writeLines()



##-----------
##  Table 2  
##-----------

library(gt)

tibble(dataset = names(datasets),
       n = map_dbl(datasets, \(x) nrow(x$X)),
       p = map_dbl(datasets, \(x) ncol(x$X)),
       K = map_dbl(datasets, "K"),
       z = c("variety", "region", "type", "variety", "class", "status", "diagnosis", "gp", "ckdmem", "class"),
       Source = c("{pgmm}", "{pgmm}", "{pgmm}", "UCI MLR", "{MoTBs}", "{MixGHD}", "UCI MLR", "CITE", "{teigen}", "{mlbench}")) |> 
  arrange(dataset) |> 
  gt() |> 
  as_latex() |> 
  as.character() |> 
  writeLines()
