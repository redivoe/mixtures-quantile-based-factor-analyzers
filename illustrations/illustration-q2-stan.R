library(cmdstanr)
library(purrr)
library(bayesplot)
library(ggplot2)
library(factor.switching)
library(coda)
library(patchwork)

library(mqfa)

load("data/illustration-q2-data.RData")

R <- 5e3
burn_in <- 5e3

compiled <- cmdstan_model("../../stan/mifa-fgld-no-z.stan")

data_list <- list("X" = data$X,
                  "n" = n,
                  "p" = p,
                  "q" = q,
                  "K" = K)

start_time <- Sys.time()
fit <- compiled$sample(data = data_list,
                       chains = 1,
                       iter_warmup = R,
                       iter_sampling = burn_in,
                       seed = 1)
end_time <- Sys.time()
elapsed_time <- difftime(end_time, start_time, units = "secs") |>
  as.numeric()


out_stan <- list()
out_stan$L <- mcmc_stan(fit, "Lambda")
out_stan$theta <- mcmc_stan(fit, c("alpha", "beta", "delta", "kappa"))
out_stan$psi <- mcmc_stan(fit, "psi")
out_stan$mu <- mcmc_stan(fit, "mu")
out_stan$prop <- mcmc_stan(fit, "prop")


mcmc_trace(fit$draws(c("delta", "kappa")), facet_args = list(nrow = 2))+
  theme_bw()+
  theme(legend.position = "none")
ggsave(filename = "output/illustration-q2-stan-traceplot-delta-kappa.pdf",
       width = 8, height = 4)


mcmc_trace(out_stan$L, facet_args = list(nrow = 4))+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_blank())
ggsave(filename = "output/illustration-q2-stan-traceplot-lambda.pdf",
       width = 14, height = 5)


##---------------------------------
##  Checking for factor switching  
##---------------------------------

rsp_names <- paste0("Lambda", paste0(rep(paste0("V", 1:p), times = q), "_", rep(1:q, each = p)))
ord <- (1 + p*(0:(q-1))) + rep(0:(p-1), each = q)

idx_Lk <- map(1:K, \(k) grep(paste0("^Lambda\\[", k, ","), x = colnames(out_stan$L)))
out_stan$L_sep <- map(idx_Lk, \(x) out_stan$L[, x] |> `colnames<-`(rsp_names))
out_rsp <- map(out_stan$L_sep, \(x) rsp_exact(lambda_mcmc = x[, ord], rotate = FALSE))

out_stan$dk <- mcmc_stan(fit, c("delta", "kappa"))
idx_dkk <- map(1:K, \(k) grep(paste0("^.*\\[", k, ","), x = colnames(out_stan$dk)))
out_stan$dk_sep <- map(idx_dkk, \(x) out_stan$dk[, x])
out_stan$theta_rsp <- list()
for(k in 1:K){
  out_stan$theta_rsp[[k]] <- array(dim = c(R, 4*q))
  for(r in 1:R){
    out_stan$theta_rsp[[k]][r, ] <- apply(ifa2::sp_fun_theta(matrix(out_stan$dk_sep[[k]][r, ], 2, q, byrow = TRUE),
                                                             out_rsp[[k]]$permute_vectors[r, ],
                                                             out_rsp[[k]]$sign_vectors[r, ]),
                                          2,
                                          \(x) complete_theta_std(x[1], x[2])) |> 
      c()
  }
}
theta_idx <- c(3, 4) + rep(seq(0, (q - 1) * 4, by = 4), each = 2)

out_stan$dk_rsp <- cbind(out_stan$theta_rsp[[1]][, theta_idx],
                         out_stan$theta_rsp[[2]][, theta_idx])[, c(seq(1, 7, by = 2), seq(2, 8, by = 2))] |> 
  `colnames<-`(colnames(out_stan$dk))

mcmc_trace(out_stan$dk_rsp, facet_args = list(nrow = 2))+
  theme_bw()+
  theme(legend.position = "none")
ggsave(filename = "output/illustration-q2-stan-traceplot-delta-kappa-rsp.pdf",
       width = 8, height = 4)

ess <- map(out_stan$theta_rsp, effectiveSize) |> 
  map_dbl(mean)

save(list = c("out_stan"),
     file = "data/illustration-q2-out-stan.RData")


  
(mcmc_trace(out_rsp[[1]]$lambda_reordered_mcmc,
           facet_args = list(nrow = 4))+
  labs(title = "Component 1")) +
(mcmc_trace(out_rsp[[2]]$lambda_reordered_mcmc,
             facet_args = list(nrow = 4))+
   labs(title = "Component 2")) &
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_blank())
ggsave(filename = "output/illustration-q2-stan-traceplot-lambda-rsp.pdf",
       width = 14, height = 5)
