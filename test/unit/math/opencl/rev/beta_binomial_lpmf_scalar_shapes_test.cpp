#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <test/unit/math/opencl/rev/beta_binomial_lpmf_reference.hpp>
#include <vector>

TEST(ProbDistributionsBetaBinomial, opencl_large_shapes_reference_scalar) {
  // Scalar shapes: alpha + beta and the shape terms are formed on the host.
  using beta_binomial_opencl_test_internal::expect_reference;
  using beta_binomial_opencl_test_internal::test_values;
  using stan::math::var;
  for (const auto& t : test_values) {
    const std::vector<int> n{t.n};
    const std::vector<int> N{t.N};
    stan::math::matrix_cl<int> n_cl(n);
    stan::math::matrix_cl<int> N_cl(N);
    var alpha = t.alpha;
    var beta = t.beta;
    var lp = stan::math::beta_binomial_lpmf(n_cl, N_cl, alpha, beta);
    lp.grad();
    expect_reference(t, lp, alpha.adj(), beta.adj(), "scalar shapes");
    stan::math::recover_memory();
  }
}

#endif
