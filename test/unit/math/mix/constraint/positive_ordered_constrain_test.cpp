#include <test/unit/math/test_ad.hpp>

namespace positive_ordered_constrain_test {
template <typename T>
T g1(const T& x) {
  stan::scalar_type_t<T> lp = 0;
  return stan::math::positive_ordered_constrain<false>(x, lp);
}
template <typename T>
T g2(const T& x) {
  stan::scalar_type_t<T> lp = 0;
  return stan::math::positive_ordered_constrain<true>(x, lp);
}
template <typename T>
stan::scalar_type_t<T> g3(const T& x) {
  typename stan::scalar_type<T>::type lp = 0;
  stan::math::positive_ordered_constrain<true>(x, lp);
  return lp;
}

template <typename T>
void expect_positive_ordered_transform(const T& x) {
  auto f1 = [](const auto& x) { return g1(x); };
  auto f2 = [](const auto& x) { return g2(x); };
  auto f3 = [](const auto& x) { return g3(x); };
  stan::test::expect_ad(f1, x);
  stan::test::expect_ad(f2, x);
  stan::test::expect_ad(f3, x);
  stan::test::expect_ad_matvar(f1, x);
  stan::test::expect_ad_matvar(f2, x);
  stan::test::expect_ad_matvar(f3, x);
}
}  // namespace positive_ordered_constrain_test

TEST(MathMixMatFun, positiveOrderedTransform) {
  Eigen::VectorXd v0(0);
  positive_ordered_constrain_test::expect_positive_ordered_transform(v0);

  Eigen::VectorXd v1(1);
  v1 << -1;
  positive_ordered_constrain_test::expect_positive_ordered_transform(v1);

  Eigen::VectorXd v2(2);
  v2 << 3, -1;
  positive_ordered_constrain_test::expect_positive_ordered_transform(v2);

  Eigen::VectorXd v3(3);
  v3 << -12, 3, -1.9;
  positive_ordered_constrain_test::expect_positive_ordered_transform(v3);

  Eigen::VectorXd v4(4);
  v4 << -1, 0, -1.1, 0.5;
  positive_ordered_constrain_test::expect_positive_ordered_transform(v4);

  Eigen::VectorXd v5(5);
  v5 << 1, -3, 2, 0, -1;
  positive_ordered_constrain_test::expect_positive_ordered_transform(v5);
}
