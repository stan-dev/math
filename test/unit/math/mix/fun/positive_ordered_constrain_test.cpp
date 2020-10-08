#include <test/unit/math/test_ad.hpp>

namespace positive_ordered_constrain_test {
template <typename T>
T g1(const T& x) {
  auto x_con = stan::math::positive_ordered_constrain(x);
  auto x_free = stan::math::positive_ordered_free(x_con);
  return x_free;
}
template <typename T>
T g2(const T& x) {
  typename stan::scalar_type<T>::type lp = 0;
  auto x_con = stan::math::positive_ordered_constrain(x, lp);
  auto x_free = stan::math::positive_ordered_free(x_con);
  return x_free;
}
template <typename T>
typename stan::scalar_type<T>::type g3(const T& x) {
  typename stan::scalar_type<T>::type lp = 0;
  auto x_con = stan::math::positive_ordered_constrain(x, lp);
  auto x_free = stan::math::positive_ordered_free(x_con);
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
