#include <test/unit/math/test_ad.hpp>
#include <limits>
#include <vector>

template <typename T, typename U>
auto g1(const T& x, const U& lb) {
  auto x_cons = stan::math::lb_constrain(x, lb);
  auto x_free = stan::math::lb_free(x_cons, lb);
  return x_free;
}
template <typename T, typename U>
auto g2(const T& x, const U& lb) {
  stan::return_type_t<T, U> lp = 0;
  auto x_cons = stan::math::lb_constrain(x, lb, lp);
  auto x_free = stan::math::lb_free(x_cons, lb);
  return x_free;
}

template <typename T, typename U>
auto g3(const T& x, const U& lb) {
  stan::return_type_t<T, U> lp = 0;
  auto x_cons = stan::math::lb_constrain(x, lb, lp);
  auto x_free = stan::math::lb_free(x_cons, lb);
  return lp;
}

template <typename T>
void expect_lb_constrain(const T& x, double lb) {
  auto f1 = [](const auto& x, const auto& lb) { return g1(x, lb); };
  auto f2 = [](const auto& x, const auto& lb) { return g2(x, lb); };
  auto f3 = [](const auto& x, const auto& lb) { return g3(x, lb); };
  stan::test::expect_ad(f1, x, lb);
  stan::test::expect_ad(f2, x, lb);
  stan::test::expect_ad(f3, x, lb);
}

TEST(mathMixScalFun, lbConstrain) {
  expect_lb_constrain(-1, 2);
  expect_lb_constrain(2, 4);
}

TEST(mathMixScalFun, lb_constrain_mat) {
  Eigen::Matrix<double, -1, -1> mat(2, 2);
  mat << -1, 1, 2, 3;
  expect_lb_constrain(mat, 2);
  expect_lb_constrain(mat, 4);
}

TEST(mathMixScalFun, lb_constrain_vec) {
  std::vector<double> vec({-1, 1, 2, 3, 4});
  expect_lb_constrain(vec, 5);
}
