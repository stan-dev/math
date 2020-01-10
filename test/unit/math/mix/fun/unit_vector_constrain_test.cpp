#include <test/unit/math/test_ad.hpp>

template <typename T>
typename Eigen::Matrix<typename stan::scalar_type<T>::type, -1, -1> g1(
    const T& x) {
  return stan::math::unit_vector_constrain(x);
}
template <typename T>
typename Eigen::Matrix<typename stan::scalar_type<T>::type, -1, -1> g2(
    const T& x) {
  typename stan::scalar_type<T>::type lp = 0;
  auto a = stan::math::unit_vector_constrain(x, lp);
  return a;
}
template <typename T>
typename stan::scalar_type<T>::type g3(const T& x) {
  typename stan::scalar_type<T>::type lp = 0;
  stan::math::unit_vector_constrain(x, lp);
  return lp;
}

template <typename T>
void expect_unit_vector_constrain(const T& x) {
  stan::test::ad_tolerances tols;
  tols.gradient_fvar_grad_ = 1e0;
  tols.hessian_grad_ = 1e0;
  tols.hessian_hessian_ = 1e0;
  tols.hessian_fvar_grad_ = 1e0;
  tols.hessian_fvar_hessian_ = 1e0;
  tols.grad_hessian_hessian_ = 1e0;
  tols.grad_hessian_grad_hessian_ = 1e1;

  auto f1 = [](const auto& x) { return g1(x); };
  auto f2 = [](const auto& x) { return g2(x); };
  auto f3 = [](const auto& x) { return g3(x); };
  stan::test::expect_ad(tols, f1, x);
  stan::test::expect_ad(tols, f2, x);
  stan::test::expect_ad(tols, f3, x);
}

TEST(MathMixMatFun, unitVectorConstrain) {
  Eigen::VectorXd v0;
  expect_unit_vector_constrain(v0);

  Eigen::VectorXd v1(1);
  v1 << 0.7;
  expect_unit_vector_constrain(v1);

  Eigen::VectorXd v3(3);
  v3 << 2, 3, -1;
  expect_unit_vector_constrain(v3);

  Eigen::VectorXd v3b(3);
  v3b << 0.6, 3, -1;
  expect_unit_vector_constrain(v3b);

  Eigen::VectorXd v6(6);
  v6 << 1, 2, -3, 1.5, 0.2, 2;
  expect_unit_vector_constrain(v6);
}
