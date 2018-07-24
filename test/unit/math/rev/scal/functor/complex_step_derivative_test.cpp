#include <gtest/gtest.h>
#include <stan/math.hpp>
#include <stan/math/rev/scal/functor/complex_step_derivative.hpp>
#include <test/unit/util.hpp>
#include <iostream>
#include <limits>
#include <sstream>
#include <vector>

struct ComplexStepDerivativeScalTest : public ::testing::Test {
  struct Fexp {
    template <typename T>
    inline T operator()(const T &theta, const std::vector<double> &x_r,
                        const std::vector<int> &x_i, std::ostream *msgs) const {
      return exp(theta) / sqrt(theta) - 0.5 * exp(theta) * pow(theta, -1.5);
    }
  };

  void SetUp() {}
  Fexp f;
  const std::vector<double> x_r;
  const std::vector<int> x_i;
  std::ostream *msgs = nullptr;
};

TEST_F(ComplexStepDerivativeScalTest, func_exp_sqrt) {
  using stan::math::complex_step_derivative;
  using stan::math::var;

  /* f near x = 0 has very large derivative */
  var x = 0.01;
  var y = complex_step_derivative(f, x, x_r, x_i, msgs);

  ASSERT_FLOAT_EQ(y.val(), f(x.val(), x_r, x_i, msgs));

  std::vector<stan::math::var> xv{x};
  std::vector<double> g1, g;
  var y1 = f(x, x_r, x_i, msgs);
  stan::math::set_zero_all_adjoints();
  y1.grad(xv, g1);
  stan::math::set_zero_all_adjoints();
  y.grad(xv, g);
  ASSERT_FLOAT_EQ(g[0], g1[0]);
}
