#include <stan/math/rev.hpp>
#include <boost/mp11.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/functor/test_fixture_ode.hpp>
#include <test/unit/math/rev/functor/test_fixture_dae_pph.hpp>
#include <test/unit/math/rev/functor/dae_test_functors.hpp>

/**
 *
 * Use same solver functor type for both w & w/o tolerance control
 */
template <typename solve_type, typename... Ts>
using ode_test_tuple = std::tuple<solve_type, solve_type, Ts...>;

/**
 * Outer product of test types
 */
using pph_test_types = boost::mp11::mp_product<
    ode_test_tuple, ::testing::Types<dae_functor>,
    ::testing::Types<double>,                                 // t
    ::testing::Types<double, stan::math::var_value<double>>,  // yy0
    ::testing::Types<double, stan::math::var_value<double>>,  // yp0
    ::testing::Types<stan::math::var_value<double>>           // theta
    >;

TYPED_TEST_SUITE_P(pph_dae_test);
TYPED_TEST_P(pph_dae_test, param_finite_diff) {
  const double h = 0.01;

  // after theta perturbation we need to make IC consistent
  this->theta += h;
  this->yy0[2] = stan::math::value_of(this->theta);
  this->yp0[1] = 1.0 - stan::math::value_of(this->yy0[2]);
  auto res_ub = this->apply_solver();

  // after theta perturbation we need to make IC consistent
  this->theta -= 2.0 * h;
  this->yy0[2] = stan::math::value_of(this->theta);
  this->yp0[1] = 1.0 - stan::math::value_of(this->yy0[2]);
  auto res_lb = this->apply_solver();

  const size_t nt = this->times().size();
  std::vector<Eigen::Matrix<stan::math::var, -1, 1>> grad_fd(nt);
  for (size_t i = 0; i < nt; ++i) {
    grad_fd[i] = (res_ub[i] - res_lb[i]) / (2.0 * h);
  }

  this->theta += h;
  this->yy0[2] = stan::math::value_of(this->theta);
  this->yp0[1] = 1.0 - stan::math::value_of(this->yy0[2]);
  auto res = this->apply_solver();
  std::vector<double> g;
  std::vector<stan::math::var> vec{this->theta};

  for (size_t i = 0; i < nt; i++) {
    for (size_t j = 0; j < res[0].size(); j++) {
      g.clear();
      res[i][j].grad(vec, g);

      double tol = 1.e-5;
      if (j == 2)
        tol = 2.e-5;
      double g_fd = stan::math::value_of(grad_fd[i][j]);

      EXPECT_NEAR(g[0], g_fd, tol)
          << "Gradient of DAE solver failed with initial positions"
          << " known and parameters unknown at time index " << i
          << ", equation index " << j;
      stan::math::set_zero_all_adjoints();
    }
  }
}

REGISTER_TYPED_TEST_SUITE_P(pph_dae_test, param_finite_diff);
INSTANTIATE_TYPED_TEST_SUITE_P(StanOde, pph_dae_test, pph_test_types);
