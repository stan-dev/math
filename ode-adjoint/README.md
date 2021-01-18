# Codes to run SBC & Speed Benchmarks for adjoint ODE

The work is dispatched using batchtools (see
https://mllg.github.io/batchtools/articles/batchtools.html
). Batchtools allows to use a cluster as backend or just your laptop.

The R codes rely on cmdstanr. So please ensure that the environment
variable CMDSTAN points to a properly setup cmdstan checkout which
itself is configured to use a math checkout of the branch
feature/adjoint-odes .

Files:

sbc.R runs SBC
sbc_report.R plots SBC results
scale-benchmark.R Speed benchmark which varries the system size
sbc_tools.R shared R functions
