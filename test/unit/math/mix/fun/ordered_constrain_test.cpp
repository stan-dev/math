#include <test/unit/math/test_ad.hpp>

namespace ordered_constrain_test {
template <typename T>
T g1(const T& x) {
  return stan::math::ordered_constrain(x);
}
template <typename T>
T g2(const T& x) {
  typename stan::scalar_type<T>::type lp = 0;
  return stan::math::ordered_constrain(x, lp);
}
template <typename T>
typename stan::scalar_type<T>::type g3(const T& x) {
  typename stan::scalar_type<T>::type lp = 0;
  stan::math::ordered_constrain(x, lp);
  return lp;
}

template <typename T>
void expect_ordered_transform(const T& x) {
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
}  // namespace ordered_constrain_test

TEST(MathMixMatFun, orderedTransform) {
  Eigen::VectorXd v0(0);
  ordered_constrain_test::expect_ordered_transform(v0);

  Eigen::VectorXd v1(1);
  v1 << -1;
  ordered_constrain_test::expect_ordered_transform(v1);

  Eigen::VectorXd v2(2);
  v2 << 3, -1;
  ordered_constrain_test::expect_ordered_transform(v2);

  Eigen::VectorXd v3(3);
  v3 << -12, 3, -1.9;
  ordered_constrain_test::expect_ordered_transform(v3);

  Eigen::VectorXd v4(4);
  v4 << -1, 0, -1.1, 0.5;
  ordered_constrain_test::expect_ordered_transform(v4);

  Eigen::VectorXd v5(5);
  v5 << 1, -3, 2, 0, -1;
  ordered_constrain_test::expect_ordered_transform(v5);
}
