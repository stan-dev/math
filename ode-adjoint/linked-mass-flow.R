##install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

library(cmdstanr)
library(posterior)
library(bayesplot)
color_scheme_set("brightblue")

set_cmdstan_path("~/work/cmdstan")

stan_model_file <- "linked-mass-flow.stan"

mod <- cmdstan_model(stan_model_file)

base_data <- list(rel_tol=1e-6,
                  abs_tol=1e-6,
                  max_num_steps=10000,
                  adjoint_integrator=0,
                  system_size=2,
                  num_obs=8,
                  sigma_sim=2.0,
                  sigma_y=0.2)

run_benchmark <- function(..., num_iter=250, adapt_delta=0.8) {
    fit_args <- list(...)
    fit_data <- modifyList(base_data, fit_args)
    print(str(fit_data))
    fit <- mod$sample(
                   data = fit_data,
                   iter_warmup = num_iter,
                   iter_sampling = num_iter,
                   seed = 123,
                   chains = 1,
                   parallel_chains = 1,
                   refresh = round(num_iter/5),
                   save_warmup=TRUE,
                   adapt_delta=adapt_delta
                   ##verbose=TRUE
               )
    fit
}

bdf_fit <- run_benchmark(adjoint_integrator=0)
adjoint_fit <- run_benchmark(adjoint_integrator=1)

adjoint_fit_2 <- run_benchmark(adjoint_integrator=1, adapt_delta=0.65)

adjoint_fit_3 <- run_benchmark(adjoint_integrator=1, adapt_delta=0.95)

adjoint_fit_4 <- run_benchmark(adjoint_integrator=1, abs_tol=1E-10, rel_tol=1E-10)

summarize_benchmark <- function(fit) {
    sdiag_all  <- fit$sampler_diagnostics(inc_warmup=TRUE)
    sdiag_main  <- fit$sampler_diagnostics(inc_warmup=FALSE)

    num_divergent <- sum(sdiag_main[,,"divergent__"])
    accept_stat <- mean(sdiag_main[,,"accept_stat__"])
    treedepth <- mean(sdiag_main[,,"treedepth__"])
    num_leaps <- sum(sdiag_all[,,"n_leapfrog__"])
    fit_time <- fit$time()$total

    ess_speed <- fit$summary("lp__")[1,c("ess_bulk", "ess_tail")] / fit_time

    list(num_leaps=num_leaps, time=fit_time, time_per_leap=1E3 * fit_time / num_leaps,
         step_size=fit$metadata()$step_size_adaptation,
         num_divergent=num_divergent, accept_stat=accept_stat, treedepth=treedepth,
         ess_lp_speed=ess_speed)
}

mcmc_rhat(rhat = rhat(bdf_fit))  + yaxis_text(hjust = 0)

str(summarize_benchmark(bdf_fit))
str(summarize_benchmark(adjoint_fit))


str(summarize_benchmark(adjoint_fit_2))
str(summarize_benchmark(adjoint_fit_3))

bdf_fit$metadata()

bdf_fit$time()


