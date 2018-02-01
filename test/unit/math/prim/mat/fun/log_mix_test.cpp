#include <stan/math/prim/mat/fun/log_mix.hpp>
#include <stan/math/prim/scal/fun/log_mix.hpp>
#include <stan/math/prim/mat/fun/log_sum_exp.hpp>
#include <stan/math/prim/mat/fun/log.hpp>
#include <stan/math/prim/mat/fun/typedefs.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

TEST(MatrixFunctions, LogMix_Values) {
  stan::math::vector_d prob(5, 1);
  prob << 0.1, 0.3, 0.25, 0.15, 0.2;

  stan::math::vector_d dens(5, 1);
  dens << -5.65, -7.62, -12.63, -55.62, -2.35;

  std::vector<double> std_dens(5);
  std_dens[0] = -5.65;
  std_dens[1] = -7.62;
  std_dens[2] = -12.63;
  std_dens[3] = -55.62;
  std_dens[4] = -2.35;

  std::vector<stan::math::vector_d> std_dens_vec{dens, dens, dens, dens};

  double log_mix_stan_1 = stan::math::log_mix(prob, std_dens);

  double log_mix_stan_2 = stan::math::log_mix(prob, std_dens_vec);

  EXPECT_FLOAT_EQ(log_mix_stan_1 * 4, log_mix_stan_2);
}