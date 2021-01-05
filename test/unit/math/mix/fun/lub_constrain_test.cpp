#include <test/unit/math/test_ad.hpp>

namespace lub_constrain_test {
template <typename T, typename L, typename U>
auto g1(const T& x, const L& lb, const U& ub) {
  return stan::math::lub_constrain(x, lb, ub);
}

template <typename T, typename L, typename U>
auto g2(const T& x, const L& lb, const U& ub) {
  stan::return_type_t<T, L, U> lp = 0;
  return stan::math::lub_constrain(x, lb, ub, lp);
}

template <typename T, typename L, typename U>
auto g3(const T& x, const L& lb, const U& ub) {
  stan::return_type_t<T, L, U> lp = 0;
  stan::math::lub_constrain(x, lb, ub, lp);
  return lp;
}
}  // namespace lub_constrain_test

TEST(mathMixMatFun, lub_constrain_scalar) {
  auto f1 = [](const auto& x, const auto& lb, const auto& ub) {
    return lub_constrain_test::g1(x, lb, ub);
  };
  auto f2 = [](const auto& x, const auto& lb, const auto& ub) {
    return lub_constrain_test::g2(x, lb, ub);
  };
  auto f3 = [](const auto& x, const auto& lb, const auto& ub) {
    return lub_constrain_test::g3(x, lb, ub);
  };

  double x1 = 0.7;
  double x2 = -1.1;
  double lb = -2.0;
  double lbi = -stan::math::INFTY;
  double ub = 1.9;
  double ubi = stan::math::INFTY;

  stan::test::expect_ad(f1, x1, lb, ub);
  stan::test::expect_ad(f1, x1, lb, ubi);
  stan::test::expect_ad(f1, x1, lbi, ub);
  stan::test::expect_ad(f1, x1, lbi, ubi);
  stan::test::expect_ad(f1, x2, lb, ub);
  stan::test::expect_ad(f1, x2, lb, ubi);
  stan::test::expect_ad(f1, x2, lbi, ub);
  stan::test::expect_ad(f1, x2, lbi, ubi);

  stan::test::expect_ad(f2, x1, lb, ub);
  stan::test::expect_ad(f2, x1, lb, ubi);
  stan::test::expect_ad(f2, x1, lbi, ub);
  stan::test::expect_ad(f2, x1, lbi, ubi);
  stan::test::expect_ad(f2, x2, lb, ub);
  stan::test::expect_ad(f2, x2, lb, ubi);
  stan::test::expect_ad(f2, x2, lbi, ub);
  stan::test::expect_ad(f2, x2, lbi, ubi);

  stan::test::expect_ad(f3, x1, lb, ub);
  stan::test::expect_ad(f3, x1, lb, ubi);
  stan::test::expect_ad(f3, x1, lbi, ub);
  stan::test::expect_ad(f3, x1, lbi, ubi);
  stan::test::expect_ad(f3, x2, lb, ub);
  stan::test::expect_ad(f3, x2, lb, ubi);
  stan::test::expect_ad(f3, x2, lbi, ub);
  stan::test::expect_ad(f3, x2, lbi, ubi);
}

TEST(mathMixMatFun, lub_mat_constrain) {
  auto f1 = [](const auto& x, const auto& lb, const auto& ub) {
    return lub_constrain_test::g1(x, lb, ub);
  };
  auto f2 = [](const auto& x, const auto& lb, const auto& ub) {
    return lub_constrain_test::g2(x, lb, ub);
  };
  auto f3 = [](const auto& x, const auto& lb, const auto& ub) {
    return lub_constrain_test::g3(x, lb, ub);
  };

  Eigen::MatrixXd x1(2, 2);
  x1 << 0.7, -1.0, 5.7, -3.8;
  Eigen::MatrixXd x2(2, 2);
  x2 << 1.1, 2.0, -0.3, 1.1;
  double lb = 2.7;
  double lbi = -stan::math::INFTY;
  double ub = -1.3;
  double ubi = stan::math::INFTY;

  stan::test::expect_ad(f1, x1, lb, ub);
  stan::test::expect_ad(f1, x1, lb, ubi);
  stan::test::expect_ad(f1, x1, lbi, ub);
  stan::test::expect_ad(f1, x1, lbi, ubi);
  /*stan::test::expect_ad(f1, x2, lb, ub);
  stan::test::expect_ad(f1, x2, lb, ubi);
  stan::test::expect_ad(f1, x2, lbi, ub);
  stan::test::expect_ad(f1, x2, lbi, ubi);*/

  /*stan::test::expect_ad(f2, x1, lb, ub);
  stan::test::expect_ad(f2, x1, lb, ubi);
  stan::test::expect_ad(f2, x1, lbi, ub);
  stan::test::expect_ad(f2, x1, lbi, ubi);
  stan::test::expect_ad(f2, x2, lb, ub);
  stan::test::expect_ad(f2, x2, lb, ubi);
  stan::test::expect_ad(f2, x2, lbi, ub);
  stan::test::expect_ad(f2, x2, lbi, ubi);

  stan::test::expect_ad(f3, x1, lb, ub);
  stan::test::expect_ad(f3, x1, lb, ubi);
  stan::test::expect_ad(f3, x1, lbi, ub);
  stan::test::expect_ad(f3, x1, lbi, ubi);
  stan::test::expect_ad(f3, x2, lb, ub);
  stan::test::expect_ad(f3, x2, lb, ubi);
  stan::test::expect_ad(f3, x2, lbi, ub);
  stan::test::expect_ad(f3, x2, lbi, ubi);*/
}
