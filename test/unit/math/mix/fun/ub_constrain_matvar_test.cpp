#include <test/unit/math/test_ad.hpp>
#include <limits>

namespace ub_constrain_test {
template <typename T1, typename T2>
void expect_matvar(const T1& x, const T2& ub) {
  auto f1 = [](const auto& x, const auto& ub) {
    stan::return_type_t<decltype(x), decltype(ub)> lp = 0;
    return stan::math::ub_constrain<false>(x, ub, lp);
  };
  auto f2 = [](const auto& x, const auto& ub) {
    stan::return_type_t<decltype(x), decltype(ub)> lp = 0;
    return stan::math::ub_constrain<true>(x, ub, lp);
  };
  auto f3 = [](const auto& x, const auto& ub) {
    stan::return_type_t<decltype(x), decltype(ub)> lp = 0;
    stan::math::ub_constrain<true>(x, ub, lp);
    return lp;
  };
  auto f4 = [](const auto& x, const auto& ub) {
    stan::return_type_t<decltype(x), decltype(ub)> lp = 0;
    auto xx = stan::math::ub_constrain<true>(x, ub, lp);
    return stan::math::add(lp, stan::math::sum(xx));
  };

  stan::test::expect_ad_matvar(f1, x, ub);
  stan::test::expect_ad_matvar(f2, x, ub);
  stan::test::expect_ad_matvar(f3, x, ub);
  stan::test::expect_ad_matvar(f4, x, ub);
}

template <typename T1, typename T2>
void expect_vec_matvar(const T1& x, const T2& ub) {
  auto f1 = [](const auto& x, const auto& ub) {
    stan::return_type_t<decltype(x), decltype(ub)> lp = 0;
    return stan::math::ub_constrain<false>(x, ub, lp);
  };
  auto f2 = [](const auto& x, const auto& ub) {
    stan::return_type_t<decltype(x), decltype(ub)> lp = 0;
    return stan::math::ub_constrain<true>(x, ub, lp);
  };
  auto f3 = [](const auto& x, const auto& ub) {
    stan::return_type_t<decltype(x), decltype(ub)> lp = 0;
    stan::math::ub_constrain<true>(x, ub, lp);
    return lp;
  };
  auto f4 = [](const auto& x, const auto& ub) {
    stan::return_type_t<decltype(x), decltype(ub)> lp = 0;
    auto xx = stan::math::eval(stan::math::ub_constrain<true>(x, ub, lp));
    stan::return_type_t<decltype(x), decltype(ub)> xx_acc = 0;
    for (size_t i = 0; i < xx.size(); ++i) {
      xx_acc += stan::math::sum(xx[i]);
    }
    return stan::math::add(lp, xx_acc);
  };
  stan::test::expect_ad_matvar(f1, x, ub);
  stan::test::expect_ad_matvar(f2, x, ub);
  stan::test::expect_ad_matvar(f3, x, ub);
  stan::test::expect_ad_matvar(f4, x, ub);
}
}  // namespace ub_constrain_test

TEST(mathMixMatFun, ub_matvar_constrain) {
  using stan::scalar_type_t;
  using stan::math::promote_scalar_t;
  using stan::math::ub_constrain;
  Eigen::MatrixXd A(2, 2);
  A << 5.0, 2.0, 0.0, 0.005;
  Eigen::MatrixXd ubm(2, 2);
  ubm << 7.0, 5.0, 0.0, 0.0005;
  Eigen::MatrixXd ubm_bad(2, 1);
  ubm_bad << 7.0, 5.0;
  ub_constrain_test::expect_matvar(A, ubm);
  ub_constrain_test::expect_matvar(A, ubm_bad);
  double ubd = 6.0;
  ub_constrain_test::expect_matvar(A, ubd);
}

TEST(mathMixMatFun, ub_matvar_constrain_neg_inf) {
  Eigen::MatrixXd A(2, 2);
  A << 5.0, 2.0, 4.0, -2.0;
  Eigen::MatrixXd ubm(2, 2);
  ubm << 7.0, 5.0, stan::math::INFTY, 100.0;
  ub_constrain_test::expect_matvar(A, ubm);
  ub_constrain_test::expect_matvar(A, stan::math::INFTY);
}

// matrix[], matrix
// matrix[], real
TEST(mathMixMatFun, ub_stdvec_mat_mat_constrain_matvar) {
  Eigen::MatrixXd A_inner(2, 3);
  A_inner << 5.0, 2.0, 4.0, -2.0, 0.0, 0.005;
  Eigen::MatrixXd ubm_inner(2, 3);
  ubm_inner << 7.0, 5.0, 6.0, 100.0, 0.0, 0.0005;
  Eigen::MatrixXd ubm_inner_bad(2, 2);
  ubm_inner_bad << 7.0, 5.0, 6.0, 100.0;
  Eigen::MatrixXd A_inner2 = 2.0 * A_inner;
  std::vector<Eigen::MatrixXd> A;
  A.push_back(A_inner);
  A.push_back(A_inner2);
  ub_constrain_test::expect_vec_matvar(A, ubm_inner);
  ub_constrain_test::expect_vec_matvar(A, ubm_inner_bad);
  double ubd = 6.0;
  ub_constrain_test::expect_vec_matvar(A, ubd);
}

TEST(mathMixMatFun, ub_stdvec_mat_mat_constrain_matvar_neg_inf) {
  Eigen::MatrixXd A_inner(2, 3);
  A_inner << 5.0, 2.0, 4.0, -2.0, 0.0, 0.005;
  Eigen::MatrixXd ubm_inner(2, 3);
  ubm_inner << 7.0, 5.0, 6.0, stan::math::INFTY, 0.0, 0.0005;
  Eigen::MatrixXd A_inner2 = 2.0 * A_inner;
  std::vector<Eigen::MatrixXd> A;
  A.push_back(A_inner);
  A.push_back(A_inner2);
  ub_constrain_test::expect_vec_matvar(A, ubm_inner);
  double ubi = stan::math::INFTY;
  ub_constrain_test::expect_vec_matvar(A, ubi);
}

// matrix[], matrix[]
TEST(mathMixMatFun, ub_stdvec_mat_constrain_matvar) {
  Eigen::MatrixXd A_inner(2, 3);
  A_inner << 5.0, 2.0, 4.0, -2.0, 0.0, 0.005;
  Eigen::MatrixXd ubm_inner(2, 3);
  ubm_inner << 7.0, 5.0, 6.0, 100.0, 0.0, 0.0005;
  Eigen::MatrixXd A_inner2 = 2.0 * A_inner;
  Eigen::MatrixXd ubm_inner2 = 3.0 * ubm_inner;
  std::vector<Eigen::MatrixXd> A;
  A.push_back(A_inner);
  A.push_back(A_inner2);
  std::vector<Eigen::MatrixXd> ubm;
  ubm.push_back(ubm_inner);
  ubm.push_back(ubm_inner2);
  std::vector<Eigen::MatrixXd> ubm_bad1;
  ubm_bad1.push_back(ubm_inner);
  Eigen::MatrixXd ubm_inner_bad(2, 2);
  ubm_inner_bad << 7.0, 5.0, 6.0, 100.0;
  std::vector<Eigen::MatrixXd> ubm_bad2;
  ubm_bad2.push_back(ubm_inner_bad);
  ub_constrain_test::expect_vec_matvar(A, ubm);
  ub_constrain_test::expect_vec_matvar(A, ubm_bad1);
  ub_constrain_test::expect_vec_matvar(A, ubm_bad2);
}

TEST(mathMixMatFun, ub_stdvec_mat_constrain_matvar_neg_inf) {
  Eigen::MatrixXd A_inner(2, 2);
  A_inner << 5.0, 2.0, 4.0, -2.0;
  Eigen::MatrixXd ubm_inner(2, 2);
  ubm_inner << 7.0, 5.0, stan::math::INFTY, 100.0;
  Eigen::MatrixXd A_inner2 = 2 * A_inner;
  Eigen::MatrixXd ubm_inner2(2, 2);
  ubm_inner << 7.0, stan::math::INFTY, 5.0, 100.0;
  std::vector<Eigen::MatrixXd> A{A_inner, A_inner2};
  std::vector<Eigen::MatrixXd> ubm{ubm_inner, ubm_inner2};
  ub_constrain_test::expect_vec_matvar(A, ubm);
  ub_constrain_test::expect_vec_matvar(A, A_inner);
  ub_constrain_test::expect_vec_matvar(A, stan::math::INFTY);
}
