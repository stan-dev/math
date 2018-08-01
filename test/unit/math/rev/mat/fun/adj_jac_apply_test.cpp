#include <stan/math/rev/core.hpp>
#include <test/unit/math/rev/mat/util.hpp>
#include <gtest/gtest.h>
#include <sstream>

struct SinFunctor {
  int N_;
  double *x_mem_;
  Eigen::VectorXd operator()(const Eigen::VectorXd &x) {
    N_ = x.size();
    Eigen::VectorXd out(N_);
    x_mem_
        = stan::math::ChainableStack::instance().memalloc_.alloc_array<double>(
            N_);

    for (int n = 0; n < N_; ++n) {
      x_mem_[n] = x(n);
      out(n) = sin(x(n));
    }

    return out;
  }

  Eigen::VectorXd multiply_adjoint_jacobian(const Eigen::VectorXd &adj) {
    Eigen::VectorXd out(N_);

    for (int n = 0; n < N_; ++n) {
      out(n) = cos(x_mem_[n]) * adj(n);
    }

    return out;
  }
};

TEST(AgradRev, test_stack) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> x1(1), x2(2), y1(1), y2(2);
  x1 << 1.0;
  x2 << 2.0, 1.0;

  y1 = stan::math::adj_jac_apply<SinFunctor>(x1);
  y2 = stan::math::adj_jac_apply<SinFunctor>(x2);

  test::check_varis_on_stack(y1);
  test::check_varis_on_stack(y2);
}

TEST(AgradRev, test_values) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> x1(1), x2(2), y1(1), y2(2);
  x1 << 1.0;
  x2 << 2.0, 1.0;

  y1 = stan::math::adj_jac_apply<SinFunctor>(x1);
  y2 = stan::math::adj_jac_apply<SinFunctor>(x2);

  EXPECT_NEAR(y1(0).val(), 0.841470984807897, 1e-10);
  EXPECT_NEAR(y2(0).val(), 0.909297426825682, 1e-10);
  EXPECT_NEAR(y2(1).val(), 0.841470984807897, 1e-10);
}

TEST(AgradRev, test_multiple_jac) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> x1(1), x2(2), y1(1), y2(2);
  x1 << 1.0;
  x2 << 2.0, 1.0;

  y1 = stan::math::adj_jac_apply<SinFunctor>(x1);
  y2 = stan::math::adj_jac_apply<SinFunctor>(x2);

  y1(0).grad();
  EXPECT_NEAR(x1(0).adj(), 0.5403023058681398, 1e-10);
  EXPECT_NEAR(x2(0).adj(), 0.0, 1e-10);
  EXPECT_NEAR(x2(1).adj(), 0.0, 1e-10);

  stan::math::set_zero_all_adjoints();

  y2(0).grad();
  EXPECT_NEAR(x1(0).adj(), 0.0, 1e-10);
  EXPECT_NEAR(x2(0).adj(), -0.4161468365471424, 1e-10);
  EXPECT_NEAR(x2(1).adj(), 0.0, 1e-10);

  stan::math::set_zero_all_adjoints();

  y2(1).grad();
  EXPECT_NEAR(x1(0).adj(), 0.0, 1e-10);
  EXPECT_NEAR(x2(0).adj(), 0.0, 1e-10);
  EXPECT_NEAR(x2(1).adj(), 0.5403023058681398, 1e-10);

  stan::math::set_zero_all_adjoints();

  stan::math::var sum_y2 = (1.73 * y2(0) + 1.57 * y2(1));
  sum_y2.grad();
  EXPECT_NEAR(x1(0).adj(), 0.0, 1e-10);
  EXPECT_NEAR(x2(0).adj(), 1.73 * -0.4161468365471424, 1e-10);
  EXPECT_NEAR(x2(1).adj(), 1.57 * 0.5403023058681398, 1e-10);
}
