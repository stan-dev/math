#include <test/unit/math/test_ad.hpp>

namespace sum_to_zero_constrain_test {
template <typename T>
inline T g1(const T& x) {
  stan::scalar_type_t<T> lp = 0;
  return stan::math::sum_to_zero_constrain<false>(x, lp);
}
template <typename T>
inline T g2(const T& x) {
  stan::scalar_type_t<T> lp = 0;
  return stan::math::sum_to_zero_constrain<true>(x, lp);
}
template <typename T>
typename stan::scalar_type<T>::type g3(const T& x) {
  stan::scalar_type_t<T> lp = 0;
  stan::math::sum_to_zero_constrain<true>(x, lp);
  return lp;
}

template <typename T>
inline void expect_sum_to_zero_transform(const T& x) {
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
}  // namespace sum_to_zero_constrain_test

TEST(MathMixMatFun, sum_to_zeroTransform) {
  Eigen::VectorXd v0(0);
  sum_to_zero_constrain_test::expect_sum_to_zero_transform(v0);

  Eigen::VectorXd v1(1);
  v1 << 1;
  sum_to_zero_constrain_test::expect_sum_to_zero_transform(v1);

  Eigen::VectorXd v2(2);
  v2 << 3, -1;
  sum_to_zero_constrain_test::expect_sum_to_zero_transform(v2);

  Eigen::VectorXd v3(3);
  v3 << 2, 3, -1;
  sum_to_zero_constrain_test::expect_sum_to_zero_transform(v3);

  Eigen::VectorXd v4(4);
  v4 << 2, -1, 0, -1.1;
  sum_to_zero_constrain_test::expect_sum_to_zero_transform(v4);

  Eigen::VectorXd v5(5);
  v5 << 1, -3, 2, 0, -1;
  sum_to_zero_constrain_test::expect_sum_to_zero_transform(v5);
}

TEST(MathMixMatFun, sum_to_zero_matrixTransform) {
  Eigen::MatrixXd m0_0(0, 0);
  sum_to_zero_constrain_test::expect_sum_to_zero_transform(m0_0);

  Eigen::MatrixXd m1_1(1, 1);
  m1_1 << 1;
  sum_to_zero_constrain_test::expect_sum_to_zero_transform(m1_1);

  Eigen::MatrixXd m2_2(2, 2);
  m2_2 << 1, 2, -3, 4;
  sum_to_zero_constrain_test::expect_sum_to_zero_transform(m2_2);

  Eigen::MatrixXd m3_4(3, 4);
  m3_4 << 1, 2, -3, 4, 5, 6, -7, 8, 9, -10, 11, -12;

  sum_to_zero_constrain_test::expect_sum_to_zero_transform(m3_4);

  Eigen::MatrixXd m4_3 = m3_4.transpose();
  sum_to_zero_constrain_test::expect_sum_to_zero_transform(m4_3);
}
