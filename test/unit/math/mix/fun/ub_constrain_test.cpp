#include <test/unit/math/test_ad.hpp>
#include <limits>

namespace ub_constrain_test {
template <typename T1, typename T2>
void expect(const T1& x, const T2& ub) {
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

  auto f4 = [](const auto& x, const auto& ub) {
    stan::return_type_t<decltype(x), decltype(ub)> lp = 0;
    auto xx = stan::math::ub_constrain(x, ub, lp);
    return stan::math::add(lp, stan::math::sum(xx));
  };
  stan::test::expect_ad(f1, x, ub);
  stan::test::expect_ad(f2, x, ub);
  stan::test::expect_ad(f3, x, ub);
  stan::test::expect_ad(f4, x, ub);
}

template <typename T1, typename T2>
void expect_vec(const T1& x, const T2& ub) {
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
  auto f4 = [](const auto& x, const auto& ub) {
    stan::return_type_t<decltype(x), decltype(ub)> lp = 0;
    auto xx = stan::math::eval(stan::math::ub_constrain(x, ub, lp));
    stan::return_type_t<decltype(x), decltype(ub)> xx_acc = 0;
    for (size_t i = 0; i < xx.size(); ++i) {
      xx_acc += stan::math::sum(xx[i]);
    }
    return stan::math::add(lp, xx_acc);
  };
  stan::test::expect_ad(f1, x, ub);
  stan::test::expect_ad(f2, x, ub);
  stan::test::expect_ad(f3, x, ub);
  stan::test::expect_ad(f4, x, ub);
}
}  // namespace ub_constrain_test

TEST(mathMixScalFun, ub_constrain) {
  ub_constrain_test::expect(-1, 2);
  ub_constrain_test::expect(2, 4);
}

TEST(mathMixScalFun, ub_constrain_inf) {
  ub_constrain_test::expect(-1, stan::math::INFTY);
  ub_constrain_test::expect(2, stan::math::INFTY);
}

TEST(mathMixMatFun, ub_mat_constrain) {
  Eigen::MatrixXd A(2, 3);
  A << -1.1, 0.0, 1.0, 2.0, 0.0, 0.0005;
  Eigen::MatrixXd ubm(2, 3);
  ubm << 1.0, 2.0, 2.0, 33.7, 0.0, 0.00005;

  double ubd1 = 10.1;
  ub_constrain_test::expect(A, ubm);
  ub_constrain_test::expect(A, ubd1);
}

TEST(mathMixMatFun, ub_mat_constrain_inf) {
  Eigen::MatrixXd A(2, 2);
  A << -1.1, 0.0, 1.0, 2.0;
  Eigen::MatrixXd ubm(2, 2);
  ubm << 1.0, 2.0, stan::math::INFTY, 33.7;

  double ubd1 = stan::math::INFTY;

  ub_constrain_test::expect(A, ubm);
  ub_constrain_test::expect(A, ubd1);
}

TEST(mathMixMatFun, ub_stdvec_constrain) {
  std::vector<double> A{5.0, 2.0, 4.0, -2.0};
  std::vector<double> ubm{7.0, 5.0, 6.0, 100.0};
  ub_constrain_test::expect(A, ubm);
  double ubd = 6.0;
  ub_constrain_test::expect(A, ubd);
}

TEST(mathMixMatFun, ub_stdvec_constrain_neg_inf) {
  std::vector<double> A{5.0, 2.0, 4.0, -2.0};
  std::vector<double> ubm{7.0, 5.0, stan::math::NEGATIVE_INFTY, 100.0};
  ub_constrain_test::expect(A, ubm);
  ub_constrain_test::expect(A, stan::math::NEGATIVE_INFTY);
}

TEST(mathMixMatFun, ub_stdvec_mat_mat_constrain) {
  Eigen::MatrixXd A_inner(2, 3);
  A_inner << 5.0, 2.0, 4.0, -2.0, 0.0, 0.005;
  Eigen::MatrixXd ubm_inner(2, 3);
  ubm_inner << 7.0, 5.0, 6.0, 100.0, 0.0, 0.0005;
  std::vector<Eigen::MatrixXd> A;
  A.push_back(A_inner);
  A.push_back(A_inner);
  // A.push_back(A_inner);
  ub_constrain_test::expect_vec(A, ubm_inner);
  double ubd = 6.0;
  ub_constrain_test::expect_vec(A, ubd);
}

TEST(mathMixMatFun, ub_stdvec_mat_constrain) {
  Eigen::MatrixXd A_inner(2, 3);
  A_inner << 5.0, 2.0, 4.0, -2.0, 0.0, 0.005;
  Eigen::MatrixXd ubm_inner(2, 3);
  ubm_inner << 7.0, 5.0, 6.0, 100.0, 0.0, 0.0005;
  std::vector<Eigen::MatrixXd> A;
  A.push_back(A_inner);
  A.push_back(A_inner);
  A.push_back(A_inner);
  std::vector<Eigen::MatrixXd> ubm;
  ubm.push_back(ubm_inner);
  ubm.push_back(ubm_inner);
  ubm.push_back(ubm_inner);
  ub_constrain_test::expect_vec(A, ubm);
  double ubd = 6.0;
  ub_constrain_test::expect_vec(A, ubd);
}

TEST(mathMixMatFun, ub_stdvec_mat_constrain_neg_inf) {
  Eigen::MatrixXd A_inner(2, 2);
  A_inner << 5.0, 2.0, 4.0, -2.0;
  Eigen::MatrixXd ubm_inner(2, 2);
  ubm_inner << 7.0, 5.0, stan::math::NEGATIVE_INFTY, 100.0;
  std::vector<Eigen::MatrixXd> A{A_inner, A_inner, A_inner};
  std::vector<Eigen::MatrixXd> ubm{ubm_inner, ubm_inner, ubm_inner};
  ub_constrain_test::expect_vec(A, ubm);
  ub_constrain_test::expect_vec(A, stan::math::NEGATIVE_INFTY);
}
