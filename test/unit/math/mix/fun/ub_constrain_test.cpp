#include <test/unit/math/test_ad.hpp>
#include <limits>

namespace ub_constrain_test {
template <typename T, typename U>
auto g1(const T& x, const U& ub) {
  return stan::math::ub_constrain(x, ub);
}
template <typename T, typename U>
auto g2(const T& x, const U& ub) {
  stan::return_type_t<T, U> lp = 0;
  return stan::math::ub_constrain(x, ub, lp);
}
template <typename T, typename U>
auto g3(const T& x, const U& ub) {
  stan::return_type_t<T, U> lp = 0;
  stan::math::ub_constrain(x, ub, lp);
  return lp;
}

void expect_ub_constrain(double x, double ub) {
  auto f1 = [](const auto& x, const auto& ub) { return g1(x, ub); };
  auto f2 = [](const auto& x, const auto& ub) { return g2(x, ub); };
  auto f3 = [](const auto& x, const auto& ub) { return g3(x, ub); };
  stan::test::expect_ad(f1, x, ub);
  stan::test::expect_ad(f2, x, ub);
  stan::test::expect_ad(f3, x, ub);
}
}  // namespace ub_constrain_test

TEST(mathMixScalFun, ub_constrain) {
  ub_constrain_test::expect_ub_constrain(-1, 2);
  ub_constrain_test::expect_ub_constrain(2, 4);
}

TEST(mathMixMatFun, ub_mat_constrain) {
  using stan::scalar_type_t;
  using stan::math::ub_constrain;
  using stan::math::promote_scalar_t;

  auto f1 = [](const auto& x, const auto& ub) {
    return ub_constrain_test::g1(x, ub);
  };

  auto f2 = [](const auto& x, const auto& ub) {
    return ub_constrain_test::g2(x, ub);
  };

  auto f3 = [](const auto& x, const auto& ub) {
    return ub_constrain_test::g3(x, ub);
  };
  
  Eigen::MatrixXd A(Eigen::MatrixXd::Random(2, 2));
  Eigen::MatrixXd ubm(2, 2);
  ubm << 1.0, -5.0, stan::math::INFTY, 100.0;

  double ubd1 = -5.0;
  double ubd2 = stan::math::INFTY;

  stan::test::expect_ad(f1, A, ubm);
  stan::test::expect_ad(f1, A, ubd1);
  stan::test::expect_ad(f1, A, ubd2);
  stan::test::expect_ad_matvar(f1, A, ubm);
  stan::test::expect_ad_matvar(f1, A, ubd1);
  stan::test::expect_ad_matvar(f1, A, ubd2);
  stan::test::expect_ad(f2, A, ubm);
  stan::test::expect_ad(f2, A, ubd1);
  stan::test::expect_ad(f2, A, ubd2);
  stan::test::expect_ad_matvar(f2, A, ubm);
  stan::test::expect_ad_matvar(f2, A, ubd1);
  stan::test::expect_ad_matvar(f2, A, ubd2);
  stan::test::expect_ad(f3, A, ubm);
  stan::test::expect_ad(f3, A, ubd1);
  stan::test::expect_ad(f3, A, ubd2);
  stan::test::expect_ad_matvar(f3, A, ubm);
  stan::test::expect_ad_matvar(f3, A, ubd1);
  stan::test::expect_ad_matvar(f3, A, ubd2);
}
