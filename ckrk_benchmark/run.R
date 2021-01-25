summary <- function(fit) {
    n.leaps <- sum(fit$sampler_diagnostics(inc_warmup=FALSE)[,,"n_leapfrog__"])
    time <- fit$time()$total
    time.per.leap <- time/n.leaps
    data.frame(num_leaps = n.leaps, time = time, time_per_leap = time.per.leap,
               time.per.ess.bulk = time/fit$summary("lp__")[1,c("ess_bulk")]$ess_bulk,
               time.per.ess.tail = time/fit$summary("lp__")[1,c("ess_tail")]$ess_tail)
}

benchmark.rk45 <- function(model.file, data.file1, data.file2, seed) {
    file <- file.path(model.file)
    mod <- cmdstan_model(file, quiet=FALSE, force_recompile=FALSE)
    fit.ckrk <- mod$sample(
                        data = data.file2,
                        iter_warmup = 1000,
                        iter_sampling = 1000,
                        seed = seed,
                        chains = 4,
                        save_warmup=TRUE,
                        parallel_chains = 1
                        )
    fit.rk45 <- mod$sample(
                        data = data.file1,
                        iter_warmup = 1000,
                        iter_sampling = 1000,
                        seed = seed,
                        chains = 4,
                        save_warmup=TRUE,
                        parallel_chains = 1
                        )
    res <- rbind(summary(fit.rk45), summary(fit.ckrk))
    res$seed <- seed
    row.names(res) <- c("rk45", "ckrk")
    res
}
