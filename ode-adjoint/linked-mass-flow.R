
if(FALSE) {
    ## run as needed to get latest cmdstanr
    install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
    ## or even install the latest dev version with profiling support
    Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
    remotes::install_github("stan-dev/cmdstanr")
}

library(cmdstanr)

set_cmdstan_path("~/work/cmdstan")

## model which scales number of states and parameters
stan_model_file <- "linked-mass-flow.stan"
## model which scales only number of states, but recycles parameters
stan_model_fixed_file <- "linked-mass-flow-fixed.stan"

mod <- cmdstan_model(stan_model_file)
mod_fixed <- cmdstan_model(stan_model_fixed_file)

rel_tol_f <- 1E-8
abs_tol_f <- 1E-8
base_data <- list(rel_tol_f=rel_tol_f,
                  abs_tol_f=abs_tol_f,
                  rel_tol_b=rel_tol_f * 100,
                  abs_tol_b=abs_tol_f * 100,
                  rel_tol_q=rel_tol_f * 1000,
                  abs_tol_q=abs_tol_f * 1000,
                  adjoint_integrator=0,
                  max_num_steps=10000,
                  num_checkpoints=150,
                  interpolation_polynomial=1,
                  solver_f=2,
                  solver_b=2,
                  system_size=2,
                  num_obs=8,
                  sigma_sim=2.0,
                  sigma_y=0.2
                  )

run_benchmark <- function(..., model, num_iter=250, adapt_delta=0.8) {
    fit_args <- list(...)
    fit_data <- modifyList(base_data, fit_args)
    print(str(fit_data))
    fit <- model$sample(
                   data = fit_data,
                   iter_warmup = num_iter,
                   iter_sampling = num_iter,
                   seed = 123,
                   chains = 1,
                   parallel_chains = 1,
                   refresh = round(num_iter/5),
                   save_warmup=TRUE,
                   adapt_delta=adapt_delta,
                   init=0.1
                   ##verbose=TRUE
               )
    fit
}


summarize_benchmark <- function(fit) {
    sdiag_all  <- fit$sampler_diagnostics(inc_warmup=TRUE)
    sdiag_main  <- fit$sampler_diagnostics(inc_warmup=FALSE)

    num_divergent <- sum(sdiag_main[,,"divergent__"])
    accept_stat <- mean(sdiag_main[,,"accept_stat__"])
    treedepth <- mean(sdiag_main[,,"treedepth__"])
    num_leaps <- sum(sdiag_all[,,"n_leapfrog__"])
    fit_time <- fit$time()$total

    ess_speed <- as.vector(fit$summary("lp__")[1,c("ess_bulk", "ess_tail")]) / fit_time
    lp_est <- as.vector(fit$summary("lp__"))

    list(num_leaps=num_leaps, time=fit_time, time_per_leap=1E3 * fit_time / num_leaps,
         step_size=fit$metadata()$step_size_adaptation,
         num_divergent=num_divergent, accept_stat=accept_stat, treedepth=treedepth,
         ess_lp_speed=ess_speed, lp_est=lp_est)
}


bdf_fit <- run_benchmark(adjoint_integrator=0, model=mod, num_iter=50)
adjoint_fit <- run_benchmark(adjoint_integrator=1, model=mod, num_iter=50)
## vary solvers, e.g.
##adjoint_fit <- run_benchmark(adjoint_integrator=1, solver_f=2, solver_b=1, num_iter=50)

str(summarize_benchmark(bdf_fit))
str(summarize_benchmark(adjoint_fit))

Lbdf_fit <- run_benchmark(adjoint_integrator=0, system_size=6, model=mod, num_iter=50)
Ladjoint_fit <- run_benchmark(adjoint_integrator=1, solver_f=2, solver_b=2, system_size=6, model=mod, num_iter=50)
Ladjoint_fit_2 <- run_benchmark(adjoint_integrator=1, solver_f=2, solver_b=2, system_size=6, model=mod, num_iter=50)

str(summarize_benchmark(Lbdf_fit))
str(summarize_benchmark(Ladjoint_fit))
str(summarize_benchmark(Ladjoint_fit_2))

## there is no gain by adjoint method if there are not more parameters
## than states
Lbdf_fit_fixed <- run_benchmark(adjoint_integrator=0, system_size=5, model=mod_fixed, num_iter=50)
Ladjoint_fit_fixed <- run_benchmark(adjoint_integrator=1, solver_f=2, solver_b=2, system_size=5, model=mod_fixed, interpolation_polynomial=1, num_iter=50)

str(summarize_benchmark(Lbdf_fit_fixed))
str(summarize_benchmark(Ladjoint_fit_fixed))

