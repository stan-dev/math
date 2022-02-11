#include <test/unit/math/test_ad.hpp>
#include <stan/math.hpp>

namespace unit_vector_constrain_test {
template <typename T>
auto g1(const T& x) {
  stan::scalar_type_t<T> lp = 0;
  return stan::math::unit_vector_constrain<false>(x, lp);
}
template <typename T>
auto g2(const T& x) {
  stan::scalar_type_t<T> lp = 0;
  auto a = stan::math::unit_vector_constrain<true>(x, lp);
  return a;
}
template <typename T>
auto g3(const T& x) {
  stan::scalar_type_t<T> lp = 0;
  stan::math::unit_vector_constrain<true>(x, lp);
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
  stan::test::expect_ad_matvar(tols, f1, x);
  stan::test::expect_ad(tols, f2, x);
  stan::test::expect_ad_matvar(tols, f2, x);
  stan::test::expect_ad(tols, f3, x);
  stan::test::expect_ad_matvar(tols, f3, x);
}
}  // namespace unit_vector_constrain_test

TEST(MathMixMatFun, unitVectorConstrain) {
  Eigen::VectorXd v0;
  unit_vector_constrain_test::expect_unit_vector_constrain(v0);

  Eigen::VectorXd v1(1);
  v1 << 0.7;
  unit_vector_constrain_test::expect_unit_vector_constrain(v1);

  Eigen::VectorXd v3(3);
  v3 << 2, 3, -1;
  unit_vector_constrain_test::expect_unit_vector_constrain(v3);

  Eigen::VectorXd v3b(3);
  v3b << 0.6, 3, -1;
  unit_vector_constrain_test::expect_unit_vector_constrain(v3b);

  Eigen::VectorXd v6(6);
  v6 << 1, 2, -3, 1.5, 0.2, 2;
  unit_vector_constrain_test::expect_unit_vector_constrain(v6);
}
