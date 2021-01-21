summary <- function(fit) {
    diag  <- fit$sampler_diagnostics(inc_warmup=FALSE)
    n.leaps <- sum(diag[,,"n_leapfrog__"])
    time <- fit$time()$chains$total

    time.per.leap <- time/colSums(diag[,,"n_leapfrog__"])

    list(num_leaps = n.leaps, elapsed = time, time_per_leap = time.per.leap,
         ess=fit$summary("lp__")[1,c("ess_bulk", "ess_tail")]
}

benchmark.rk45 <- function(model.file, data.file1, data.file2, seed) {
    file <- file.path(model.file)
    mod <- cmdstan_model(file, quiet=FALSE)
    fit.rk45 <- mod$sample(
                        data = data.file1,
                        iter_warmup = 1000,
                        iter_sampling = 1000,
                        seed = seed,
                        chains = 4,
                        parallel_chains = 1
                        )
    fit.erk45 <- mod$sample(
                        data = data.file2,
                        iter_warmup = 1000,
                        iter_sampling = 1000,
                        seed = seed,
                        chains = 4,
                        parallel_chains = 1
                        )
    list(rk45=summary(fit.rk45), erk45=summary(fit.erk45))
}
