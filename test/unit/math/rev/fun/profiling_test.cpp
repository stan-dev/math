#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/fun/util.hpp>
#include <boost/random/mersenne_twister.hpp>

stan::math::profiles profiles;
TEST(Profiling, simple) {
  using stan::math::cholesky_decompose;
  using stan::math::matrix_d;
  using stan::math::matrix_v;
  using stan::math::start_profiling;
  using stan::math::stop_profiling;
  using stan::math::var;
  using stan::math::wishart_rng;

  boost::random::mt19937 rng;
  matrix_d I(1500, 1500);
  I.setZero();
  I.diagonal().setOnes();

  matrix_v Y = wishart_rng(5000, I, rng);
  matrix_v PP = Eigen::VectorXd::Ones(1500);

  start_profiling(0, profiles);
  matrix_v L = cholesky_decompose(Y);
  stop_profiling(0, profiles);

  start_profiling(1, profiles);
  matrix_v P = L * PP;
  stop_profiling(1, profiles);

  P(0, 0).grad();

  stan::math::print_profiling(profiles);
}
