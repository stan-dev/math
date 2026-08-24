#include <stan/math/mix/functor/hessian_times_vector_variadic.hpp>
#include <gtest/gtest.h>
#include <type_traits>
#include <utility>

namespace hessian_times_vector_variadic_test {

struct cubic_functor {
  template <typename T0, typename T1, typename T2>
  auto operator()(const T0& theta, std::ostream* /* msgs */, const T1& u,
                  const T2& scale) const {
    auto u_dot_theta = u(0) * theta(0);
    for (Eigen::Index i = 1; i < theta.size(); ++i) {
      u_dot_theta += u(i) * theta(i);
    }
    return scale * u_dot_theta * u_dot_theta * u_dot_theta / 6.0;
  }
};

template <typename ThetaVec, typename VVec, typename = void>
struct is_hessian_times_vector_callable : std::false_type {};

template <typename ThetaVec, typename VVec>
struct is_hessian_times_vector_callable<
    ThetaVec, VVec,
    std::void_t<decltype(stan::math::hessian_times_vector(
        cubic_functor{}, std::declval<const ThetaVec&>(),
        std::declval<const VVec&>(), nullptr,
        std::declval<const Eigen::VectorXd&>(), 1.0))>> : std::true_type {};

using fvar_vector
    = Eigen::Matrix<stan::math::fvar<double>, Eigen::Dynamic, 1>;
static_assert(is_hessian_times_vector_callable<Eigen::Vector3d,
                                               Eigen::Vector3d>::value);
static_assert(!is_hessian_times_vector_callable<fvar_vector,
                                                fvar_vector>::value);

}  // namespace hessian_times_vector_variadic_test

TEST(MathMixFunctor, HessianTimesVectorVariadicArithmetic) {
  Eigen::VectorXd theta(3);
  Eigen::VectorXd v(3);
  Eigen::VectorXd u(3);
  theta << 0.4, -0.7, 1.2;
  v << 1.1, -0.3, 0.6;
  u << 0.8, -0.5, 0.2;
  const double scale = 1.3;

  Eigen::VectorXd hv = stan::math::hessian_times_vector(
      hessian_times_vector_variadic_test::cubic_functor{}, theta, v, nullptr,
      u, scale);
  Eigen::VectorXd expected = scale * u.dot(theta) * u.dot(v) * u;
  EXPECT_TRUE(hv.isApprox(expected, 1e-12));
}

TEST(MathMixFunctor, HessianTimesVectorVariadicReverse) {
  using stan::math::var;
  Eigen::VectorXd theta(3);
  Eigen::VectorXd v(3);
  Eigen::VectorXd u(3);
  Eigen::VectorXd w(3);
  theta << 0.4, -0.7, 1.2;
  v << 1.1, -0.3, 0.6;
  u << 0.8, -0.5, 0.2;
  w << -0.5, 0.8, 0.2;
  const double scale = 1.3;

  Eigen::Matrix<var, Eigen::Dynamic, 1> theta_var = theta;
  Eigen::Matrix<var, Eigen::Dynamic, 1> v_var = v;
  Eigen::Matrix<var, Eigen::Dynamic, 1> u_var = u;
  var scale_var = scale;
  auto hv = stan::math::hessian_times_vector(
      hessian_times_vector_variadic_test::cubic_functor{}, theta_var, v_var,
      nullptr, u_var, scale_var);
  var target = w.dot(hv);
  target.grad();

  const double ut = u.dot(theta);
  const double uv = u.dot(v);
  const double uw = u.dot(w);
  Eigen::VectorXd expected = scale * ut * uv * u;
  Eigen::VectorXd theta_adj = scale * uv * uw * u;
  Eigen::VectorXd v_adj = scale * ut * uw * u;
  Eigen::VectorXd u_adj
      = scale * (uv * uw * theta + ut * uw * v + ut * uv * w);
  EXPECT_TRUE(hv.val().isApprox(expected, 1e-12));
  EXPECT_TRUE(theta_var.adj().isApprox(theta_adj, 1e-12));
  EXPECT_TRUE(v_var.adj().isApprox(v_adj, 1e-12));
  EXPECT_TRUE(u_var.adj().isApprox(u_adj, 1e-12));
  EXPECT_NEAR(scale_var.adj(), ut * uv * uw, 1e-12);
  stan::math::recover_memory();
}

TEST(MathMixFunctor, HessianTimesVectorVariadicReverseVarMatrix) {
  using stan::math::var;
  using stan::math::var_value;
  Eigen::VectorXd theta(3);
  Eigen::VectorXd v(3);
  Eigen::VectorXd u(3);
  Eigen::VectorXd w(3);
  theta << 0.4, -0.7, 1.2;
  v << 1.1, -0.3, 0.6;
  u << 0.8, -0.5, 0.2;
  w << -0.5, 0.8, 0.2;
  const double scale = 1.3;

  var_value<Eigen::VectorXd> theta_var = theta;
  var_value<Eigen::VectorXd> v_var = v;
  var_value<Eigen::VectorXd> u_var = u;
  auto hv = stan::math::hessian_times_vector(
      hessian_times_vector_variadic_test::cubic_functor{}, theta_var, v_var,
      nullptr, u_var, scale);
  var target = w.dot(hv);
  target.grad();

  const double ut = u.dot(theta);
  const double uv = u.dot(v);
  const double uw = u.dot(w);
  Eigen::VectorXd expected = scale * ut * uv * u;
  Eigen::VectorXd theta_adj = scale * uv * uw * u;
  Eigen::VectorXd v_adj = scale * ut * uw * u;
  Eigen::VectorXd u_adj
      = scale * (uv * uw * theta + ut * uw * v + ut * uv * w);
  EXPECT_TRUE(hv.val().isApprox(expected, 1e-12));
  EXPECT_TRUE(theta_var.adj().isApprox(theta_adj, 1e-12));
  EXPECT_TRUE(v_var.adj().isApprox(v_adj, 1e-12));
  EXPECT_TRUE(u_var.adj().isApprox(u_adj, 1e-12));
  stan::math::recover_memory();
}

TEST(MathMixFunctor, HessianTimesVectorVariadicVarArgument) {
  using stan::math::var;
  Eigen::VectorXd theta(3);
  Eigen::VectorXd v(3);
  Eigen::VectorXd u(3);
  Eigen::VectorXd w(3);
  theta << 0.4, -0.7, 1.2;
  v << 1.1, -0.3, 0.6;
  u << 0.8, -0.5, 0.2;
  w << -0.5, 0.8, 0.2;
  var scale = 1.3;

  auto hv = stan::math::hessian_times_vector(
      hessian_times_vector_variadic_test::cubic_functor{}, theta, v, nullptr,
      u, scale);
  w.dot(hv).grad();
  EXPECT_NEAR(scale.adj(), u.dot(theta) * u.dot(v) * u.dot(w), 1e-12);
  stan::math::recover_memory();
}

TEST(MathMixFunctor, HessianTimesVectorVariadicChecks) {
  Eigen::VectorXd theta(2);
  Eigen::VectorXd v(3);
  Eigen::VectorXd u(2);
  theta << 0.4, -0.7;
  v << 1.1, -0.3, 0.6;
  u << 0.8, -0.5;
  EXPECT_THROW(stan::math::hessian_times_vector(
                   hessian_times_vector_variadic_test::cubic_functor{}, theta,
                   v, nullptr, u, 1.3),
               std::invalid_argument);

  Eigen::VectorXd empty(0);
  EXPECT_EQ(stan::math::hessian_times_vector(
                hessian_times_vector_variadic_test::cubic_functor{}, empty,
                empty, nullptr, empty, 1.3)
                .size(),
            0);
}
