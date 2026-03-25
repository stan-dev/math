using BridgeStan

stan_code = """
functions {
  real partial_sum(array[] real slice, int start, int end,
                   tuple(real, int) params) {
    real mu = params.1;
    int K = params.2;
    real lp = 0;
    for (i in 1:size(slice)) {
      lp += normal_lpdf(slice[i] | mu, K);
    }
    return lp;
  }
}
data {
  int N;
  int K;
  array[N] real y;
}
parameters {
  real mu;
}
model {
  mu ~ normal(0, 10);
  target += reduce_sum(partial_sum, y, 1, (mu, K));
}
"""

stan_math = ENV["MWE_RUN_DIR"]
label = ENV["MWE_LABEL"]

workdir = mktempdir()
stan_file = joinpath(workdir, "mwe.stan")
write(stan_file, stan_code)

lib = compile_model(stan_file; make_args=["STAN_MATH=$stan_math", "STAN_THREADS=true"])

data = """{"N": 5, "K": 2, "y": [1.0, 2.0, 3.0, 4.0, 5.0]}"""
sm = StanModel(lib, data)

params = [3.0]  # mu (unconstrained)
lp = log_density(sm, params)
lp_grad, grad = log_density_gradient(sm, params)

println("[$label] log_density = $lp")
println("[$label] gradient = $grad")

@assert isfinite(lp) "log_density should be finite"
@assert all(isfinite, grad) "gradient should be finite"
@assert lp ≈ lp_grad "log_density values should match"
@assert length(grad) == 1 "gradient should have 1 element"

println("[$label] All checks passed!")
