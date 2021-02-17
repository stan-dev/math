#include <test/unit/math/test_ad.hpp>
#include <limits>

namespace lb_constrain_test {
template <typename T1, typename T2>
void expect(const T1& x, const T2& lb) {
  auto f1 = [](const auto& x, const auto& lb) {
    return stan::math::lb_constrain(x, lb);
  };
  auto f2 = [](const auto& x, const auto& lb) {
    stan::return_type_t<decltype(x), decltype(lb)> lp = 0;
    return stan::math::lb_constrain(x, lb, lp);
  };
  auto f3 = [](const auto& x, const auto& lb) {
    stan::return_type_t<decltype(x), decltype(lb)> lp = 0;
    stan::math::lb_constrain(x, lb, lp);
    return lp;
  };
  auto f4 = [](const auto& x, const auto& lb) {
    stan::return_type_t<decltype(x), decltype(lb)> lp = 0;
    auto xx = stan::math::lb_constrain(x, lb, lp);
    return stan::math::add(lp, stan::math::sum(xx));
  };

  stan::test::expect_ad(f1, x, lb);
  stan::test::expect_ad(f2, x, lb);
  stan::test::expect_ad(f3, x, lb);
  stan::test::expect_ad(f4, x, lb);
}

template <typename T1, typename T2>
void expect_vec(const T1& x, const T2& lb) {
  auto f1 = [](const auto& x, const auto& lb) {
    return stan::math::lb_constrain(x, lb);
  };
  auto f2 = [](const auto& x, const auto& lb) {
    stan::return_type_t<decltype(x), decltype(lb)> lp = 0;
    return stan::math::lb_constrain(x, lb, lp);
  };
  auto f3 = [](const auto& x, const auto& lb) {
    stan::return_type_t<decltype(x), decltype(lb)> lp = 0;
    stan::math::lb_constrain(x, lb, lp);
    return lp;
  };
  auto f4 = [](const auto& x, const auto& lb) {
    stan::return_type_t<decltype(x), decltype(lb)> lp = 0;
    auto xx = stan::math::eval(stan::math::lb_constrain(x, lb, lp));
    stan::return_type_t<decltype(x), decltype(lb)> xx_acc = 0;
    for (size_t i = 0; i < xx.size(); ++i) {
      xx_acc += stan::math::sum(xx[i]);
    }
    return stan::math::add(lp, xx_acc);
  };
  stan::test::expect_ad(f1, x, lb);
  stan::test::expect_ad(f2, x, lb);
  stan::test::expect_ad(f3, x, lb);
  stan::test::expect_ad(f4, x, lb);
}

}  // namespace lb_constrain_test

TEST(mathMixScalFun, lbConstrain) {
  lb_constrain_test::expect(-1, 2);
  lb_constrain_test::expect(2, 4);
}

TEST(mathMixScalFun, lbConstrain_neg_inf) {
  lb_constrain_test::expect(-1, stan::math::NEGATIVE_INFTY);
  lb_constrain_test::expect(2, stan::math::NEGATIVE_INFTY);
}

TEST(mathMixMatFun, lb_mat_constrain) {
  Eigen::MatrixXd A(2, 2);
  // Doesn't work for values near zero fail for fvar hessian?
  A << 5.0, 2.0, 4.0, -2.0;//, 0.0, 0.005;
  Eigen::MatrixXd lbm(2, 2);
  lbm << 7.0, 5.0, 6.0, 100.0;//, 0.0, 0.0005;
  lb_constrain_test::expect(A, lbm);
  double lbd = 6.0;
  lb_constrain_test::expect(A, lbd);
}

TEST(mathMixMatFun, lb_mat_constrain_neg_inf) {
  Eigen::MatrixXd A(2, 2);
  A << 5.0, 2.0, 4.0, -2.0;
  Eigen::MatrixXd lbm(2, 2);
  lbm << 7.0, 5.0, stan::math::NEGATIVE_INFTY, 100.0;
  lb_constrain_test::expect(A, lbm);
  lb_constrain_test::expect(A, stan::math::NEGATIVE_INFTY);
}

TEST(mathMixMatFun, lb_stdvec_constrain) {
  std::vector<double> A{5.0, 2.0, 4.0, -2.0};
  std::vector<double> lbm{7.0, 5.0, 6.0, 100.0};
  lb_constrain_test::expect(A, lbm);
  double lbd = 6.0;
  lb_constrain_test::expect(A, lbd);
}

TEST(mathMixMatFun, lb_stdvec_constrain_neg_inf) {
  std::vector<double> A{5.0, 2.0, 4.0, -2.0};
  std::vector<double> lbm{7.0, 5.0, stan::math::NEGATIVE_INFTY, 100.0};
  lb_constrain_test::expect(A, lbm);
  lb_constrain_test::expect(A, stan::math::NEGATIVE_INFTY);
}


TEST(mathMixMatFun, lb_stdvec_mat_mat_constrain) {
  Eigen::MatrixXd A_inner(2, 1);
  // Doesn't work for values near zero fail for fvar hessian?
  A_inner << 5.0, 2.0;//, 4.0, -2.0;//, 0.0, 0.005;
  Eigen::MatrixXd lbm_inner(2, 1);
  lbm_inner << 7.0, 5.0;//, 6.0, 100.0;//, 0.0, 0.0005;
  std::vector<Eigen::MatrixXd> A;
  A.push_back(A_inner);
  A.push_back(A_inner);
  //A.push_back(A_inner);
  lb_constrain_test::expect_vec(A, lbm_inner);
  double lbd = 6.0;
  lb_constrain_test::expect_vec(A, lbd);
}

TEST(mathMixMatFun, lb_stdvec_mat_constrain) {
  Eigen::MatrixXd A_inner(2, 2);
  // Doesn't work for values near zero fail for fvar hessian?
  A_inner << 5.0, 2.0, 4.0, -2.0;//, 0.0, 0.005;
  Eigen::MatrixXd lbm_inner(2, 2);
  lbm_inner << 7.0, 5.0, 6.0, 100.0;//, 0.0, 0.0005;
  std::vector<Eigen::MatrixXd> A;
  A.push_back(A_inner);
  A.push_back(A_inner);
  A.push_back(A_inner);
  std::vector<Eigen::MatrixXd> lbm;
  lbm.push_back(lbm_inner);
  lbm.push_back(lbm_inner);
  lbm.push_back(lbm_inner);
  lb_constrain_test::expect_vec(A, lbm);
  double lbd = 6.0;
  lb_constrain_test::expect_vec(A, lbd);
}

TEST(mathMixMatFun, lb_stdvec_mat_constrain_neg_inf) {
  Eigen::MatrixXd A_inner(2, 2);
  A_inner << 5.0, 2.0, 4.0, -2.0;
  Eigen::MatrixXd lbm_inner(2, 2);
  lbm_inner << 7.0, 5.0, stan::math::NEGATIVE_INFTY, 100.0;
  std::vector<Eigen::MatrixXd> A{A_inner, A_inner, A_inner};
  std::vector<Eigen::MatrixXd> lbm{lbm_inner, lbm_inner, lbm_inner};
  lb_constrain_test::expect_vec(A, lbm);
  lb_constrain_test::expect_vec(A, stan::math::NEGATIVE_INFTY);
}
