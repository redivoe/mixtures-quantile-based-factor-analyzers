library(mqfa)

n <- 300
K <- 3
q <- 1
p <- 10

mu <- rbind(diag(1.5, K),
            diag(-1.5, K)[, K:1],
            matrix(0, p - 2*K, K))

L <- array(dim = c(p, q, K))
L[, , 1] <- c(1, 0.5, 0, -0.5, 1, 0, 0, 1, 0, 0)
L[, , 2] <- c(0, 1, 0, 0, 1, 0.5, 1, 0, 0, 0)
L[, , 3] <- c(1, 0, -1, 0, 0, 1, 0, 0, 1, 0)

theta <- array(dim = c(4, q, K))
theta[, , 1] <- complete_theta_std(delta = 0.5, kappa = 0.2)
theta[, , 2] <- complete_theta_std(delta = 0.9, kappa = 0.2)
theta[, , 3] <- complete_theta_std(delta = 0.5, kappa = 5)

psi <- matrix(data = c(0.19, 0.01, 0.19, 0.2, 0.26, 0.2, 0, 0.06, 0.07, 0.29,
                       0.03, 0.01, 0.12, 0.05, 0.46, 0.07, 0.04, 0.33, 0.53, 0.09,
                       0.28, 0.02, 0, 0.06, 0.07, 0.45, 0.1, 0.18, 0.06, 0.09),
              p, K)

prop <- c(0.2, 0.3, 0.5)

params <- list("mu" = mu,
               "L" = L,
               "theta" = theta,
               "psi" = psi,
               "prop" = prop,
               "common_theta" = FALSE,
               "common_psi" = FALSE)

set.seed(321)
data <- r_mifa_fgld_data(n = n, p = p, q = q, K = K, params = params)
# pairs(data$X, col = data$params$z)
# table(data$params$z)

params_init <- mifa_random_init(n, p, q, K)

save(list = c("data", "params_init", "n", "p", "q", "K"),
     file = "data/illustration-q1-data.RData")

