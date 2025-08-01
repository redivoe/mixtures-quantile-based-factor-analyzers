library(purrr)

##--------------------
##  Loading datasets
##--------------------

datasets <- list()

data("coffee", package = "pgmm")
datasets$coffee <- list(X = coffee[,-c(1, 2)] |> scale(),
                        z = unclass(coffee$Variety))

data("olive", package = "pgmm")
datasets$olive <- list(X = olive[, -c(1, 2)] |> scale(),
                       z = olive$Region)

data("wine", package = "pgmm")
datasets$wine <- list(X =  wine[, -1] |> scale(),
                      z = wine$Type)

data("seeds", package = "datasetsICR")
datasets$seeds <- list(X = seeds[,-ncol(seeds)] |> scale(),
                       z = unclass(seeds$variety))

data("ecoli", package = "MoTBFs")
datasets$ecoli <- list(X = ecoli[, -c(1, 4, 5, 9)] |> scale(),
                       z = unclass(factor(ecoli$class)))
# lip and chg (cols 4 and 5) only have 2 unique values

data("banknote", package = "MixGHD")
datasets$banknote <- list(X = banknote[, -c(1)] |> scale(),
                          z = unclass(factor(banknote$Status)))


# wbcd <- read.csv(
#     url("https://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/wdbc.data"),
#     header = FALSE,
#     na.strings = "?")
# save(list = "wbcd", file = "wbcd.RData")
load("data/wbcd.RData") # pre-downloaded

datasets$wbcd <- list(X = wbcd[,-c(1, 2)] |> scale(),
                      z = unclass(wbcd[, 2]))

# SOURCE
# http://faculty.washington.edu/kayee/model/
yeast <- read.table("data/norm_cellcycle_384_17.txt", header = TRUE)
datasets$yeast <- list(X = yeast[,-c(1, 2)] |> scale(),
                       z = unclass(yeast$Gp))


data("ckd", package = "teigen")
datasets$ckd <- list(X = ckd[, -1] |> scale(),
                     z = unclass(ckd$ckdmem))


data("Vehicle", package = "mlbench")
datasets$vehicle <- list(X = Vehicle[,-ncol(Vehicle)] |> scale(),
                         z = unclass(Vehicle$Class))


K <- map(datasets, \(x) length(unique(x$z)))
datasets <- map2(datasets, K, \(x, y){x$K <- y; x})

upper_bound_factors <- function(p){
  q <- 1
  while(0.5 * ((p - q)^2 - p - q) > 0){
    q <- q + 1
  }
  return(q-1)
}
p <- map(datasets, \(x) ncol(x$X))
q_max <- map_dbl(p, upper_bound_factors) |>
  pmin(5) |>
  pmax(1)
datasets <- map2(datasets, q_max, \(x, y){x$q_max <- y; x})

save(list = "datasets", file = "data/datasets.RData")
