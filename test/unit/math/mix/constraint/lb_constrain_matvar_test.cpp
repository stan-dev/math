#include <test/unit/math/test_ad.hpp>
#include <limits>

namespace lb_constrain_test {
template <typename T1, typename T2>
void expect_matvar(const T1& x, const T2& lb) {
  auto f1 = [](const auto& x, const auto& lb) {
    stan::return_type_t<decltype(x), decltype(lb)> lp = 0;
    return stan::math::lb_constrain<false>(x, lb, lp);
  };
  auto f2 = [](const auto& x, const auto& lb) {
    stan::return_type_t<decltype(x), decltype(lb)> lp = 0;
    return stan::math::lb_constrain<true>(x, lb, lp);
  };
  auto f3 = [](const auto& x, const auto& lb) {
    stan::return_type_t<decltype(x), decltype(lb)> lp = 0;
    stan::math::lb_constrain<true>(x, lb, lp);
    return lp;
  };
  auto f4 = [](const auto& x, const auto& lb) {
    stan::return_type_t<decltype(x), decltype(lb)> lp = 0;
    auto xx = stan::math::lb_constrain<true>(x, lb, lp);
    return stan::math::add(lp, stan::math::sum(xx));
  };

  stan::test::expect_ad_matvar(f1, x, lb);
  stan::test::expect_ad_matvar(f2, x, lb);
  stan::test::expect_ad_matvar(f3, x, lb);
  stan::test::expect_ad_matvar(f4, x, lb);
}

template <typename T1, typename T2>
void expect_vec_matvar(const T1& x, const T2& lb) {
  auto f1 = [](const auto& x, const auto& lb) {
    stan::return_type_t<decltype(x), decltype(lb)> lp = 0;
    return stan::math::lb_constrain<false>(x, lb, lp);
  };
  auto f2 = [](const auto& x, const auto& lb) {
    stan::return_type_t<decltype(x), decltype(lb)> lp = 0;
    return stan::math::lb_constrain<true>(x, lb, lp);
  };
  auto f3 = [](const auto& x, const auto& lb) {
    stan::return_type_t<decltype(x), decltype(lb)> lp = 0;
    stan::math::lb_constrain<true>(x, lb, lp);
    return lp;
  };
  auto f4 = [](const auto& x, const auto& lb) {
    stan::return_type_t<decltype(x), decltype(lb)> lp = 0;
    auto xx = stan::math::eval(stan::math::lb_constrain<true>(x, lb, lp));
    stan::return_type_t<decltype(x), decltype(lb)> xx_acc = 0;
    for (size_t i = 0; i < xx.size(); ++i) {
      xx_acc += stan::math::sum(xx[i]);
    }
    return stan::math::add(lp, xx_acc);
  };
  stan::test::expect_ad_matvar(f1, x, lb);
  stan::test::expect_ad_matvar(f2, x, lb);
  stan::test::expect_ad_matvar(f3, x, lb);
  stan::test::expect_ad_matvar(f4, x, lb);
}
}  // namespace lb_constrain_test

TEST(mathMixMatFun, lb_matvar_constrain) {
  using stan::scalar_type_t;
  using stan::math::lb_constrain;
  using stan::math::promote_scalar_t;
  Eigen::MatrixXd A(2, 2);
  A << 5.0, 2.0, 0.0, 0.005;
  Eigen::MatrixXd lbm(2, 2);
  lbm << 7.0, 5.0, 0.0, 0.0005;
  Eigen::MatrixXd lbm_bad(2, 1);
  lbm_bad << 7.0, 5.0;
  lb_constrain_test::expect_matvar(A, lbm);
  lb_constrain_test::expect_matvar(A, lbm_bad);
  double lbd = 6.0;
  lb_constrain_test::expect_matvar(A, lbd);
}

TEST(mathMixMatFun, lb_matvar_constrain_neg_inf) {
  Eigen::MatrixXd A(2, 2);
  A << 5.0, 2.0, 4.0, -2.0;
  Eigen::MatrixXd lbm(2, 2);
  lbm << 7.0, 5.0, stan::math::NEGATIVE_INFTY, 100.0;
  lb_constrain_test::expect_matvar(A, lbm);
  lb_constrain_test::expect_matvar(A, stan::math::NEGATIVE_INFTY);
}

// matrix[], matrix
// matrix[], real
TEST(mathMixMatFun, lb_stdvec_mat_mat_constrain_matvar) {
  Eigen::MatrixXd A_inner(2, 3);
  A_inner << 5.0, 2.0, 4.0, -2.0, 0.0, 0.005;
  Eigen::MatrixXd lbm_inner(2, 3);
  lbm_inner << 7.0, 5.0, 6.0, 100.0, 0.0, 0.0005;
  Eigen::MatrixXd lbm_inner_bad(2, 2);
  lbm_inner_bad << 7.0, 5.0, 6.0, 100.0;
  Eigen::MatrixXd A_inner2 = 2.0 * A_inner;
  std::vector<Eigen::MatrixXd> A;
  A.push_back(A_inner);
  A.push_back(A_inner2);
  lb_constrain_test::expect_vec_matvar(A, lbm_inner);
  lb_constrain_test::expect_vec_matvar(A, lbm_inner_bad);
  double lbd = 6.0;
  lb_constrain_test::expect_vec_matvar(A, lbd);
}

TEST(mathMixMatFun, lb_stdvec_mat_mat_constrain_matvar_neg_inf) {
  Eigen::MatrixXd A_inner(2, 3);
  A_inner << 5.0, 2.0, 4.0, -2.0, 0.0, 0.005;
  Eigen::MatrixXd lbm_inner(2, 3);
  lbm_inner << 7.0, 5.0, 6.0, stan::math::NEGATIVE_INFTY, 0.0, 0.0005;
  Eigen::MatrixXd A_inner2 = 2.0 * A_inner;
  std::vector<Eigen::MatrixXd> A;
  A.push_back(A_inner);
  A.push_back(A_inner2);
  lb_constrain_test::expect_vec_matvar(A, lbm_inner);
  double lbi = stan::math::NEGATIVE_INFTY;
  lb_constrain_test::expect_vec_matvar(A, lbi);
}

// matrix[], matrix[]
TEST(mathMixMatFun, lb_stdvec_mat_constrain_matvar) {
  Eigen::MatrixXd A_inner(2, 3);
  A_inner << 5.0, 2.0, 4.0, -2.0, 0.0, 0.005;
  Eigen::MatrixXd lbm_inner(2, 3);
  lbm_inner << 7.0, 5.0, 6.0, 100.0, 0.0, 0.0005;
  Eigen::MatrixXd A_inner2 = 2.0 * A_inner;
  Eigen::MatrixXd lbm_inner2 = 3.0 * lbm_inner;
  std::vector<Eigen::MatrixXd> A;
  A.push_back(A_inner);
  A.push_back(A_inner2);
  std::vector<Eigen::MatrixXd> lbm;
  lbm.push_back(lbm_inner);
  lbm.push_back(lbm_inner2);
  std::vector<Eigen::MatrixXd> lbm_bad1;
  lbm_bad1.push_back(lbm_inner);
  Eigen::MatrixXd lbm_inner_bad(2, 2);
  lbm_inner_bad << 7.0, 5.0, 6.0, 100.0;
  std::vector<Eigen::MatrixXd> lbm_bad2;
  lbm_bad2.push_back(lbm_inner_bad);
  lb_constrain_test::expect_vec_matvar(A, lbm);
  lb_constrain_test::expect_vec_matvar(A, lbm_bad1);
  lb_constrain_test::expect_vec_matvar(A, lbm_bad2);
}

TEST(mathMixMatFun, lb_stdvec_mat_constrain_matvar_neg_inf) {
  Eigen::MatrixXd A_inner(2, 2);
  A_inner << 5.0, 2.0, 4.0, -2.0;
  Eigen::MatrixXd lbm_inner(2, 2);
  lbm_inner << 7.0, 5.0, stan::math::NEGATIVE_INFTY, 100.0;
  Eigen::MatrixXd A_inner2 = 2 * A_inner;
  Eigen::MatrixXd lbm_inner2(2, 2);
  lbm_inner << 7.0, stan::math::NEGATIVE_INFTY, 5.0, 100.0;
  std::vector<Eigen::MatrixXd> A{A_inner, A_inner2};
  std::vector<Eigen::MatrixXd> lbm{lbm_inner, lbm_inner2};
  lb_constrain_test::expect_vec_matvar(A, lbm);
  lb_constrain_test::expect_vec_matvar(A, A_inner);
  lb_constrain_test::expect_vec_matvar(A, stan::math::NEGATIVE_INFTY);
}
