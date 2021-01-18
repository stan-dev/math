library(cmdstanr)
library(posterior)
library(batchtools)
set.seed(453453)
## load utilities and current dev version of OncoBayes2
source("sbc_tools.R")

#' according to the docs this speeds up the reduce step
options(batchtools.progress = FALSE)

registry_tmp <- Sys.getenv("TMP_NFS", tempdir())

reg <- makeExperimentRegistry(
    file.dir = tempfile("sbc-", registry_tmp),
    ## use the default configured batchtools configuration batchtools.conf
    ## - found in the environemnt variable R_BATCHTOOLS_SEARCH_PATH
    ## - current working directory
    ## - $(HOME)/.batchtools.conf
    ## conf.file = NA,
    seed = 47845854,
    ## our worker functions and package loading
    source="sbc_tools.R")

## resources of each job: Less than 55min, 2000MB RAM and 1 core
job_resources <- list(walltime=55, memory=2000, ncpus=1, max.concurrent.jobs=300)

## choose number of sample as being dividable by 2 minus one. This gives us exactly by 2 divisable number of bins.
base_data  <- list(num_warmup=2^8 - 1, num_sampling=2^8 - 1, adapt_delta=0.9, stan_model=stan_model)

addProblem("linked_mass_flow",
           data = base_data,
           fun = setup_system,
           seed = 2345,
           ## caching speeds up the reduce step
           cache = TRUE
)

addAlgorithm("stan_cvodes", run_benchmark)

pdes <- list(linked_mass_flow = data.frame(system_size=c(1,4)))
ades <- list(stan_cvodes= data.frame(adjoint_integrator=c(0,1)))

##addExperiments(pdes, ades, repls = 1000)
addExperiments(pdes, ades, repls = 2)

summarizeExperiments()

## approximate per job runtime in s
target_runtime <- 100
## for 30 min runtime we then need
chunk_size <- round(30 * 60 / target_runtime)
chunk_size

##chunk_size <- 10

ids <- unwrap(getJobPars())

ids[,chunk:=chunk(job.id, chunk.size=chunk_size)]

auto_submit(ids, reg, job_resources)


if(FALSE) {
job1 <- testJob(4000)

names(job1)


job1$summary("lp__")
job1$draws("lp__")

job1$metadata()$iter_sampling

as.list(job1$rank)

job1_def <- makeJob(4000)

job1_def$instance

job1_def$repl

names(job1_def$problem$data)
}

calibration <- ijoin(
  ## grab job parameters
  unwrap(getJobPars()),
    unwrap(reduceResultsDataTable(fun=function(x) c(as.list(x$rank),
                                      ess_speed=as.list(x$ess_lp_speed),
                                      x[c("num_leaps", "time", "time_per_leap", "num_divergent")]
                                                    )) )
    )



head(calibration)

saveRDS(calibration, file = "calibration.rds")

removeRegistry(0)

sessionInfo()
