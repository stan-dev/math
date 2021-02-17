#include <test/unit/math/test_ad.hpp>
#include <limits>

namespace ub_constrain_test {
template <typename T1, typename T2>
void expect_ub_constrain_matvar(const T1& x, const T2& ub) {
  auto f1 = [](const auto& x, const auto& ub) {
    return stan::math::ub_constrain(x, ub);
  };
  auto f2 = [](const auto& x, const auto& ub) {
    stan::return_type_t<decltype(x), decltype(ub)> lp = 0;
    return stan::math::ub_constrain(x, ub, lp);
  };
  auto f3 = [](const auto& x, const auto& ub) {
    stan::return_type_t<decltype(x), decltype(ub)> lp = 0;
    stan::math::ub_constrain(x, ub, lp);
    return lp;
  };
  stan::test::expect_ad_matvar(f1, x, ub);
  stan::test::expect_ad_matvar(f2, x, ub);
  stan::test::expect_ad_matvar(f3, x, ub);
}
}  // namespace ub_constrain_test

TEST(mathMixMatFun, ub_matvar_constrain) {
  using stan::scalar_type_t;
  using stan::math::promote_scalar_t;
  using stan::math::ub_constrain;

  Eigen::MatrixXd A(2, 3);
  A << -1.1, 0.0, 1.0, 2.0, 0.0, 0.0005;
  Eigen::MatrixXd ubm(2, 3);
  ubm << 1.0, 2.0, 2.0, 33.7, 0.0, 0.00005;

  double ubd1 = 10.1;
  ub_constrain_test::expect_ub_constrain_matvar(A, ubm);
  ub_constrain_test::expect_ub_constrain_matvar(A, ubd1);
}

TEST(mathMixMatFun, ub_matvar_constrain_inf) {

  Eigen::MatrixXd A(2, 2);
  A << -1.1, 0.0, 1.0, 2.0;
  Eigen::MatrixXd ubm(2, 2);
  ubm << 1.0, 2.0, stan::math::INFTY, 33.7;

  double ubd1 = stan::math::INFTY;

  ub_constrain_test::expect_ub_constrain_matvar(A, ubm);
  ub_constrain_test::expect_ub_constrain_matvar(A, ubd1);
}
