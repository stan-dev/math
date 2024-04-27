#include <stan/math/rev.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(RevFunctor, partialDerivsAllScalar) {
  using stan::math::partial_derivatives;

  auto f = [](const auto& x1, const auto& x2) {
    return x1 * x1 * x2 + 3.0 * x2 * x2;
  };
  double x1 = 5;
  double x2 = 7;
  std::tuple<double, std::tuple<double, double>> res
    = partial_derivatives(f, x1, x2);
  double fx = std::get<0>(res);
  std::tuple<double, double> grads = std::get<1>(res);

  EXPECT_FLOAT_EQ(f(x1, x2), fx);
  EXPECT_FLOAT_EQ(2 * x1 * x2, std::get<0>(grads));
  EXPECT_FLOAT_EQ(x1 * x1 + 3 * 2 * x2, std::get<1>(grads));
}

TEST(RevFunctor, partialDerivsVecScalar) {
  using stan::math::partial_derivatives;
  using stan::math::to_vector;
  using stan::math::square;
  using stan::math::sum;

  auto f = [](const auto& x1, const auto& x2) {
    auto x1_arr = stan::math::as_array_or_scalar(x1);
    return sum(square(x1_arr) * x2 + 3.0 * x2 * x2);
  };

  Eigen::VectorXd x1_vec(3);
  x1_vec << 2, 0.4, 10.8;
  double x2 = 7;

  std::tuple<double, std::tuple<Eigen::VectorXd, double>> res_vec
    = partial_derivatives(f, x1_vec, x2);
  double fx_array = std::get<0>(res_vec);
  std::tuple<Eigen::VectorXd, double> grads_vec = std::get<1>(res_vec);

  EXPECT_FLOAT_EQ(f(x1_vec, x2), fx_array);
  EXPECT_MATRIX_EQ(2 * x1_vec * x2, std::get<0>(grads_vec));
  auto x1_arr = stan::math::as_array_or_scalar(x1_vec);
  EXPECT_FLOAT_EQ(sum(square(x1_arr) + 3 * 2 * x2),
                  std::get<1>(grads_vec));
}
