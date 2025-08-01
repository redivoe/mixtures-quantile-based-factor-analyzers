library(purrr)
library(mqfa)

n <- 300
K <- 2
p <- 10
q <- 2

params <- list()
params$mu <- rbind(diag(1.5, K),
                   diag(-1.5, K)[, K:1],
                   matrix(0, p - 2*K, K))

# lambda_values <- c(0, 2, -2, 3, -3)
# lambds_probs <- c(0.5, 0.2, 0.2, 0.05, 0.05)
# get_L <- function(p, q, K, values, probs){
#   array(data = sample(values, size = p*q*K, replace = TRUE, prob = probs),
#         dim = c(p, q, K))
# }
# set.seed(123)
# params$L <- get_L(p, q, K, lambda_values, lambds_probs)

set.seed(1)
params$L <- array(data = rnorm(p * q * K, sd = 2), dim = c(p, q, K))

theta <- array(dim = c(4, q, K))
tt2 <- complete_theta_std(delta = 0.99, kappa = 0.1)
tt3 <- complete_theta_std(delta = 0.5, kappa = 0.1)
tt4 <- complete_theta_std(delta = 0.5, kappa = 10)
tt5 <- complete_theta_std(delta = 0.5, kappa = 0.1)


theta[, , 1] <- cbind(tt2, tt3)
theta[, , 2] <- cbind(tt4, tt5)
params$theta <- theta
params$common_theta <- FALSE

params$prop <- c(0.5, 0.5)

params$psi <- matrix(data = rexp(p, 5),
                     ncol = K)
params$common_psi <- FALSE


set.seed(1)
data <- r_mifa_fgld_data(n, p, q, K, params = params)
# table(data$params$z)
# plot(data$params$Y[data$params$z == 1, ])
# plot(data$params$Y[data$params$z == 2, ])

save(list = c("data", "n", "p", "q", "K"), file = "data/illustration-q2-data.RData")
