#include <test/unit/math/test_ad.hpp>

namespace stochastic_column_constrain_test {
template <typename T>
T g1(const T& x) {
  stan::scalar_type_t<T> lp = 0;
  return stan::math::stochastic_column_constrain<false>(x, lp);
}
template <typename T>
T g2(const T& x) {
  stan::scalar_type_t<T> lp = 0;
  return stan::math::stochastic_column_constrain<true>(x, lp);
}
template <typename T>
typename stan::scalar_type<T>::type g3(const T& x) {
  stan::scalar_type_t<T> lp = 0;
  stan::math::stochastic_column_constrain<true>(x, lp);
  return lp;
}

template <typename T>
void expect_simplex_transform(const T& x) {
  auto f1 = [](const auto& x) { return g1(x); };
  auto f2 = [](const auto& x) { return g2(x); };
  auto f3 = [](const auto& x) { return g3(x); };
  stan::test::expect_ad(f1, x);
  stan::test::expect_ad_matvar(f1, x);
  stan::test::expect_ad(f2, x);
  stan::test::expect_ad_matvar(f2, x);
  stan::test::expect_ad(f3, x);
  stan::test::expect_ad_matvar(f3, x);
}
}  // namespace stochastic_column_constrain_test

TEST(MathMixMatFun, simplexColumnTransform) {
  Eigen::MatrixXd v0(0, 0);
  stochastic_column_constrain_test::expect_simplex_transform(v0);

  Eigen::MatrixXd v1(1, 3);
  v1 << 1, 2, 3;
  stochastic_column_constrain_test::expect_simplex_transform(v1);

  Eigen::MatrixXd v2(2, 3);
  v2 << 3, -1, 3, -1, 3, -1;
  stochastic_column_constrain_test::expect_simplex_transform(v2);

  Eigen::MatrixXd v3(3, 3);
  v3 << 2, 3, -1, 2, 3, -1, 2, 3, -1;
  stochastic_column_constrain_test::expect_simplex_transform(v3);

  Eigen::MatrixXd v4(4, 3);
  v4 << 2, -1, 0, -1.1, 2, -1, 0, -1.1, 2, -1, 0, -1.1;
  stochastic_column_constrain_test::expect_simplex_transform(v4);

  Eigen::MatrixXd v5(5, 3);
  v5 << 1, -3, 2, 0, -1, 1, -3, 2, 0, -1, 1, -3, 2, 0, -1;
  stochastic_column_constrain_test::expect_simplex_transform(v5);
}
