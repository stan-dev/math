#include <stan/math/rev/functor/solve_itp.hpp>
#include <gtest/gtest.h>
#include <Eigen/Core>
#include <cmath>
#include <ostream>
#include <vector>

namespace {

struct cubic_parameter_functor {
  template <typename T_x, typename T_theta>
  auto operator()(const T_x& x, std::ostream*, const T_theta& theta) const {
    return x * x * x - theta;
  }
};

struct affine_parameter_functor {
  template <typename T_x, typename T_scale, typename T_offset>
  auto operator()(const T_x& x, std::ostream*, const T_scale& scale,
                  const T_offset& offset) const {
    return scale * x - offset;
  }
};

struct vector_parameter_functor {
  template <typename T_x, typename T_theta>
  auto operator()(const T_x& x, std::ostream*, const T_theta& theta) const {
    return theta(0) * x + theta(1) - theta(2);
  }
};

struct singular_parameter_functor {
  template <typename T_x, typename T_theta>
  auto operator()(const T_x& x, std::ostream*, const T_theta& theta) const {
    return x * x - theta;
  }
};

}  // namespace

TEST(MathRevFunctorSolveItp, ImplicitDerivativeScalarParameter) {
  using stan::math::var;
  var theta = 2.0;
  var root = stan::math::solve_itp(cubic_parameter_functor(), 0.0, 2.0, nullptr,
                                   theta);

  std::vector<var> inputs{theta};
  std::vector<double> gradient;
  root.grad(inputs, gradient);

  const double expected_root = std::cbrt(2.0);
  const double expected_derivative
      = 1.0 / (3.0 * expected_root * expected_root);
  EXPECT_NEAR(expected_root, root.val(), 1e-10);
  EXPECT_NEAR(expected_derivative, gradient[0], 1e-9);
  stan::math::recover_memory();
}

TEST(MathRevFunctorSolveItp, ImplicitDerivativeVariadicParameters) {
  using stan::math::var;
  var scale = 3.0;
  var offset = 2.0;
  var root = stan::math::solve_itp(affine_parameter_functor(), -10.0, 10.0,
                                   nullptr, scale, offset);

  std::vector<var> inputs{scale, offset};
  std::vector<double> gradient;
  root.grad(inputs, gradient);

  EXPECT_NEAR(2.0 / 3.0, root.val(), 1e-10);
  EXPECT_NEAR(-2.0 / 9.0, gradient[0], 1e-10);
  EXPECT_NEAR(1.0 / 3.0, gradient[1], 1e-10);
  stan::math::recover_memory();
}

TEST(MathRevFunctorSolveItp, ImplicitDerivativeVectorParameter) {
  using stan::math::var;
  Eigen::Matrix<var, Eigen::Dynamic, 1> theta(3);
  theta << 2.0, 1.0, 5.0;
  var root = stan::math::solve_itp(vector_parameter_functor(), -10.0, 10.0,
                                   nullptr, theta);

  std::vector<var> inputs{theta(0), theta(1), theta(2)};
  std::vector<double> gradient;
  root.grad(inputs, gradient);

  EXPECT_NEAR(2.0, root.val(), 1e-10);
  EXPECT_NEAR(-1.0, gradient[0], 1e-10);
  EXPECT_NEAR(-0.5, gradient[1], 1e-10);
  EXPECT_NEAR(0.5, gradient[2], 1e-10);
  stan::math::recover_memory();
}

TEST(MathRevFunctorSolveItp, RejectsSingularImplicitDerivative) {
  using stan::math::var;
  var theta = 0.0;
  EXPECT_THROW(stan::math::solve_itp(singular_parameter_functor(), -1.0, 0.0,
                                     nullptr, theta),
               std::domain_error);
  stan::math::recover_memory();
}
