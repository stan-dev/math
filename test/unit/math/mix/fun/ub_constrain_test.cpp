#include <test/unit/math/test_ad.hpp>
#include <limits>

template <typename T, typename U>
auto g1(const T& x, const U& ub) {
  auto x_cons = stan::math::ub_constrain(x, ub);
  auto x_free = stan::math::ub_free(x_cons, ub);
  return x_free;
}
template <typename T, typename U>
auto g2(const T& x, const U& ub) {
  stan::return_type_t<T, U> lp = 0;
  auto x_cons = stan::math::ub_constrain(x, ub, lp);
  auto x_free = stan::math::ub_free(x_cons, ub);
  return x_free;
}
template <typename T, typename U>
auto g3(const T& x, const U& ub) {
  stan::return_type_t<T, U> lp = 0;
  auto x_cons = stan::math::ub_constrain(x, ub, lp);
  auto x_free = stan::math::ub_free(x_cons, ub);
  return lp;
}

template <typename T>
void expect_ub_constrain(const T& x, double ub) {
  auto f1 = [](const auto& x, const auto& ub) { return g1(x, ub); };
  auto f2 = [](const auto& x, const auto& ub) { return g2(x, ub); };
  auto f3 = [](const auto& x, const auto& ub) { return g3(x, ub); };
  stan::test::expect_ad(f1, x, ub);
  stan::test::expect_ad(f2, x, ub);
  stan::test::expect_ad(f3, x, ub);
}

TEST(mathMixScalFun, ub_constrain_scalar) {
  expect_ub_constrain(-1, 2);
  expect_ub_constrain(2, 4);
}

TEST(mathMixScalFun, ub_constrain_mat) {
  Eigen::Matrix<double, -1, -1> mat(2, 2);
  mat << -1, 1, 2, 3;
  expect_ub_constrain(mat, 2);
  expect_ub_constrain(mat, 4);
}

TEST(mathMixScalFun, ub_constrain_vec) {
  std::vector<double> vec({-1, 1, 2, 3, 4});
  expect_ub_constrain(vec, 5);
}
