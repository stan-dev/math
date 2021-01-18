library(cmdstanr)
library(posterior)
library(bayesplot)
color_scheme_set("brightblue")

if(!("sbc_dir" %in%  ls())) {
    sbc_dir <- getwd()
}

stan_model_file <- normalizePath(file.path(sbc_dir, "linked-mass-flow.stan"))

stan_model <- cmdstan_model(stan_model_file)

setup_system <- function(data, job, ...) {
    args <- list(...)
    base_stan_data <- list(rel_tol=1e-6,
                           abs_tol=1e-6,
                           max_num_steps=10000,
                           adjoint_integrator=0,
                           system_size=1,
                           num_obs=8,
                           sigma_sim=1.0,
                           sigma_y=0.2)
    modifyList(base_stan_data, args)
}

run_benchmark <- function(data, job, instance, adjoint_integrator, ...) {
    fit_data <- modifyList(instance, list(adjoint_integrator=adjoint_integrator))
    r <- job$repl
    fit <- data$stan_model$sample(
        data = fit_data,
        iter_warmup = data$num_warmup,
        iter_sampling = data$num_sampling,
        seed = 123 + r,
        init=0,
        chains = 1,
        parallel_chains = 1,
        refresh = 100,
        save_warmup=TRUE,
        adapt_delta=data$adapt_delta
        ##verbose=TRUE
        )
    summarize_benchmark(fit)
}


summarize_benchmark <- function(fit) {
    sdiag_all  <- fit$sampler_diagnostics(inc_warmup=TRUE)
    sdiag_main  <- fit$sampler_diagnostics(inc_warmup=FALSE)

    num_divergent <- sum(sdiag_main[,,"divergent__"])
    accept_stat <- mean(sdiag_main[,,"accept_stat__"])
    treedepth <- mean(sdiag_main[,,"treedepth__"])
    num_leaps <- sum(sdiag_all[,,"n_leapfrog__"])
    fit_time <- fit$time()$total

    ess_speed <- as.vector(unlist(fit$summary("lp__")[1,c("ess_bulk", "ess_tail")])) / fit_time
    names(ess_speed) <- c("bulk", "tail")
    lp_est <- as.vector(fit$summary("lp__"))

    rank_vars <- grep("^rank", fit$metadata()$stan_variables, value=TRUE)
    ranks <- fit$summary(rank_vars)
    num_iter <- fit$metadata()$iter_sampling

    count_ranks <- ranks$mean * num_iter
    names(count_ranks) <- ranks$variable

    list(num_leaps=num_leaps, time=fit_time, time_per_leap=1E3 * fit_time / num_leaps,
         step_size=fit$metadata()$step_size_adaptation,
         num_divergent=num_divergent, accept_stat=accept_stat, treedepth=treedepth,
         ess_lp_speed=ess_speed, lp_est=lp_est, ranks=count_ranks)
}


## Submits to batchtools cluster with fault tolerance, i.e.
## resubmitting failed jobs max_num_tries times
auto_submit <- function(jobs, registry, resources=list(), max_num_tries = 10) {
  all_unfinished_jobs <- jobs

  num_unfinished_jobs <- nrow(all_unfinished_jobs)
  num_all_jobs <- num_unfinished_jobs
  remaining_tries <- max_num_tries
  all_jobs_finished <- FALSE
  while (remaining_tries > 0 && !all_jobs_finished) {
    remaining_tries <- remaining_tries - 1

    message("Submitting jobs at ", Sys.time())
    # Once things run fine let's submit this work to the cluster.
    submitJobs(all_unfinished_jobs, resources=resources)
    # Wait for results.
    waitForJobs()
    message("Finished waiting for jobs at ", Sys.time())

    # Check status:
    print(getStatus())

    # Ensure that all jobs are done
    if (nrow(findNotDone()) != 0) {
      not_done_jobs <- findNotDone()
      print(getErrorMessages(not_done_jobs))
      ##browser()
      ##invisible(readline(prompt="Press [enter] to continue"))

      message("Some jobs did not complete. Please check the batchtools registry ", registry$file.dir)
      all_unfinished_jobs <- inner_join(not_done_jobs, all_unfinished_jobs)

      if (num_unfinished_jobs == nrow(all_unfinished_jobs) &&  nrow(all_unfinished_jobs) > 0.25 * num_all_jobs)
      {
        # Unfinished job count did not change -> retrying will probably not help. Abort!
        warning("Error: unfinished job count is not decreasing. Aborting job retries.")
        remaining_tries <- 0
      }

      if (num_unfinished_jobs == nrow(jobs))
      {
        # All jobs errored -> retrying will probably not help. Abort!
        warning("Error: all jobs errored. Aborting job retries.")
        remaining_tries <- 0
      }

      num_unfinished_jobs <- nrow(all_unfinished_jobs)
      message("Trying to resubmit jobs. Remaining tries: ", remaining_tries, " / ", max_num_tries)
    } else {
      all_jobs_finished <- TRUE
    }
  }

  invisible(all_jobs_finished)
}


reduce_run <- function(x) {
    c(as.list(x$rank),
      ess_speed=as.list(x$ess_lp_speed),
      x[c("num_leaps", "time", "time_per_leap", "num_divergent")]
      )
}

reduce_run_bench <- function(x) {
    c(
      ess_speed=as.list(x$ess_lp_speed),
      x[c("num_leaps", "time", "time_per_leap", "num_divergent")]
      )
}

