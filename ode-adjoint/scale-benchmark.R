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
    file.dir = tempfile("scale-bench-", registry_tmp),
    ## use the default configured batchtools configuration batchtools.conf
    ## - found in the environemnt variable R_BATCHTOOLS_SEARCH_PATH
    ## - current working directory
    ## - $(HOME)/.batchtools.conf
    ## conf.file = NA,
    seed = 47845854,
    ## our worker functions and package loading
    source="sbc_tools.R")

## resources of each job: Less than 4h, 2000MB RAM and 1 core
job_resources <- list(walltime=4*60, memory=2000, ncpus=1, max.concurrent.jobs=300)

## choose number of sample as being dividable by 2 minus one. This gives us exactly by 2 divisable number of bins.
base_data  <- list(num_warmup=20, num_sampling=20, adapt_delta=0.9, stan_model=stan_model)

addProblem("linked_mass_flow",
           data = base_data,
           fun = setup_system,
           seed = 2345,
           ## caching speeds up the reduce step
           cache = TRUE
)

addAlgorithm("stan_cvodes", run_benchmark)

pdes <- list(linked_mass_flow = data.frame(system_size=c(1,2,4,8,16)))
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

##ids[,chunk:=chunk(job.id, chunk.size=chunk_size)]

##auto_submit(ids, reg, job_resources)

submitJobs(ids)

getStatus()

getLog(3)

if(FALSE) {

    job1 <- testJob(3)

    names(job1)

    job1


job1$summary("lp__")
job1$draws("lp__")

job1$metadata()$iter_sampling

as.list(job1$rank)

job1_def <- makeJob(4)

job1_def$instance

job1_def$repl

names(job1_def$problem$data)
}


library(ggplot2)
library(dplyr)

scale_benchmark <- ijoin(
  ## grab job parameters
  unwrap(getJobPars()),
  unwrap(reduceResultsDataTable(fun=reduce_run_bench))
)

head(scale_benchmark)

scale_benchmark <- mutate(scale_benchmark,
                          adjoint_integrator=factor(adjoint_integrator))

saveRDS(scale_benchmark, file = "scale_benchmark.rds")

removeRegistry(0)

sessionInfo()

scale_benchmark <- as.data.frame(scale_benchmark)
scale_benchmark <- mutate(scale_benchmark, repl=factor(rep(1:2, 10)), set=paste(repl,adjoint_integrator,sep="-"))


scale_benchmark_ref <- left_join(scale_benchmark, filter(scale_benchmark, adjoint_integrator==0)[c("repl", "system_size", "time_per_leap")], suffix=c("", "_ref"), by=c("repl", "system_size"))

scale_benchmark_ref

theme_set(theme_bw())
ggplot(scale_benchmark, aes(system_size, time_per_leap, group=set, colour=adjoint_integrator, shape=repl)) +
   geom_point()  + geom_line() +
    scale_x_log10(breaks=c(1,2,4,8,16)) +
    scale_y_log10() +
    ggtitle("Time per leapfrog of adjoint and forward per replication")

ggsave("scale_time_per_leap.png", width=1.6 *4, height=4, dpi=120)

ggplot(filter(scale_benchmark_ref, adjoint_integrator==1), aes(system_size, time_per_leap_ref/time_per_leap, group=set, shape=repl)) +
   geom_point()  + geom_line() +
    scale_x_log10(breaks=c(1,2,4,8,16)) +
    scale_y_log10() +
    geom_hline(yintercept=1.0, linetype=2) +
    ggtitle("Speedup of adjoint vs forward per replication")

ggsave("scale_speedup.png", width=1.6 *4, height=4, dpi=120)
