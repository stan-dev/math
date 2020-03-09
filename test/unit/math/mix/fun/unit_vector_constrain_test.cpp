#include <test/unit/math/test_ad.hpp>
#include <vector>

template <typename T>
auto g1(const T& x) {
  auto x_cons = stan::math::unit_vector_constrain(x);
  auto x_free = stan::math::unit_vector_free(x_cons);
  return x_free;
}
template <typename T>
auto g2(const T& x) {
  stan::scalar_type_t<T> lp = 0;
  auto x_cons = stan::math::unit_vector_constrain(x, lp);
  auto x_free = stan::math::unit_vector_free(x_cons);
  return x_free;
}
template <typename T>
auto g3(const T& x) {
  stan::scalar_type_t<T> lp = 0;
  auto x_cons = stan::math::unit_vector_constrain(x, lp);
  auto x_free = stan::math::unit_vector_free(x_cons);
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

TEST(MathMixMatFun, unitVectorConstrainEig) {
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

TEST(MathMixMatFun, unitVectorConstrainVec) {
  std::vector<double> v0(0);
  expect_unit_vector_constrain(v0);

  std::vector<double> v1({1});
  expect_unit_vector_constrain(v1);

  std::vector<double> v2({3, -1});
  expect_unit_vector_constrain(v2);

  std::vector<double> v3({-12, 3, -1.9});
  expect_unit_vector_constrain(v3);

  std::vector<double> v4({-1, 0, -1.1, 0.5});
  expect_unit_vector_constrain(v4);

  std::vector<double> v5({1, -3, 2, 0, -1});
  expect_unit_vector_constrain(v5);
}
