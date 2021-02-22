#include <test/unit/math/test_ad.hpp>

namespace lub_constrain_tests {
template <typename T1, typename T2, typename T3>
void expect(const T1& x, const T2& lb, const T3& ub) {
  auto f1 = [](const auto& x, const auto& lb, const auto& ub) {
    return stan::math::lub_constrain(x, lb, ub);
  };
  auto f2 = [](const auto& x, const auto& lb, const auto& ub) {
    stan::return_type_t<decltype(x), decltype(lb), decltype(ub)> lp = 0;
    return stan::math::lub_constrain(x, lb, ub, lp);
  };
  auto f3 = [](const auto& x, const auto& lb, const auto& ub) {
    stan::return_type_t<decltype(x), decltype(lb), decltype(ub)> lp = 0;
    stan::math::lub_constrain(x, lb, ub, lp);
    return lp;
  };
  auto f4 = [](const auto& x, const auto& lb, const auto& ub) {
    stan::return_type_t<decltype(x), decltype(lb), decltype(ub)> lp = 0;
    auto xx = stan::math::lub_constrain(x, lb, ub, lp);
    return stan::math::add(lp, stan::math::sum(xx));
  };

  stan::test::expect_ad(f1, x, lb, ub);
  stan::test::expect_ad(f2, x, lb, ub);
  stan::test::expect_ad(f3, x, lb, ub);
  stan::test::expect_ad(f4, x, lb, ub);
}
template <typename T1, typename T2, typename T3>
void expect_vec(const T1& x, const T2& lb, const T3& ub) {
  auto f1 = [](const auto& x, const auto& lb, const auto& ub) {
    return stan::math::lub_constrain(x, lb, ub);
  };
  auto f2 = [](const auto& x, const auto& lb, const auto& ub) {
    stan::return_type_t<decltype(x), decltype(lb), decltype(ub)> lp = 0;
    return stan::math::lub_constrain(x, lb, ub, lp);
  };
  auto f3 = [](const auto& x, const auto& lb, const auto& ub) {
    stan::return_type_t<decltype(x), decltype(lb), decltype(ub)> lp = 0;
    stan::math::lub_constrain(x, lb, ub, lp);
    return lp;
  };
  auto f4 = [](const auto& x, const auto& lb, const auto& ub) {
    stan::return_type_t<decltype(x), decltype(lb), decltype(ub)> lp = 0;
    auto xx = stan::math::lub_constrain(x, lb, ub, lp);
    stan::return_type_t<decltype(x), decltype(lb), decltype(ub)> xx_acc = 0;
    for (size_t i = 0; i < xx.size(); ++i) {
      xx_acc += stan::math::sum(xx[i]);
    }
    return stan::math::add(lp, xx_acc);
  };

  stan::test::expect_ad(f1, x, lb, ub);
  stan::test::expect_ad(f2, x, lb, ub);
  stan::test::expect_ad(f3, x, lb, ub);
  stan::test::expect_ad(f4, x, lb, ub);
}

}

TEST(mathMixMatFun, lub_constrain_scalars) {

  double x1 = 0.7;
  double x2 = -38.1;
  double lb = -2.0;
  double ub = 3.5;
  lub_constrain_tests::expect(x1, lb, ub);
  lub_constrain_tests::expect(x2, lb, ub);
  lub_constrain_tests::expect(x1, lb, lb);
  lub_constrain_tests::expect(x2, lb, lb);
  // ub inf
  auto ub_inf = stan::math::INFTY;
  lub_constrain_tests::expect(x1, lb, ub_inf);
  lub_constrain_tests::expect(x2, lb, ub_inf);

  // lb inf
  auto lb_inf = stan::math::NEGATIVE_INFTY;
  lub_constrain_tests::expect(x1, lb_inf, ub);
  lub_constrain_tests::expect(x2, lb_inf, ub);

  // both inf
  lub_constrain_tests::expect(x1, lb_inf, ub_inf);
  lub_constrain_tests::expect(x2, lb_inf, ub_inf);

}

TEST(mathMixMatFun, lub_constrain_vector_scalar_scalar) {

  Eigen::MatrixXd x1(1, 1);
  x1 << 0.7;
  Eigen::MatrixXd x2(1, 1);
  x2 << -1.1;
  double lb = -2.0;
  double ub = 3.5;
  lub_constrain_tests::expect(x1, lb, ub);
  lub_constrain_tests::expect(x2, lb, ub);
  lub_constrain_tests::expect(x1, lb, lb);
  lub_constrain_tests::expect(x2, lb, lb);


  // ub inf
  auto ub_inf = stan::math::INFTY;
  lub_constrain_tests::expect(x1, lb, ub_inf);
  lub_constrain_tests::expect(x2, lb, ub_inf);

  // lb inf
  auto lb_inf = stan::math::NEGATIVE_INFTY;
  lub_constrain_tests::expect(x1, lb_inf, ub);
  lub_constrain_tests::expect(x2, lb_inf, ub);

  // both inf
  lub_constrain_tests::expect(x1, lb_inf, ub_inf);
  lub_constrain_tests::expect(x2, lb_inf, ub_inf);

}

TEST(mathMixMatFun, lub_constrain_vector_vector_scalar) {

  Eigen::MatrixXd x1(2, 2);
  x1 << 5.0, 2.0, 4.0, 5.0;
  Eigen::MatrixXd x2(2, 2);
  x2 << -1.1, 0.005, 1.0, 3.0;
  Eigen::MatrixXd lb(2, 2);
  lb << -3.0, 5.0, -6.0, 6.0;
  double ub = 13.5;
  lub_constrain_tests::expect(x1, lb, ub);
  lub_constrain_tests::expect(x2, lb, ub);

  // ub inf
  auto ub_inf = stan::math::INFTY;
  lub_constrain_tests::expect(x1, lb, ub_inf);
  lub_constrain_tests::expect(x2, lb, ub_inf);


  // lb inf
  lb(1, 1) =  stan::math::NEGATIVE_INFTY;
  lub_constrain_tests::expect(x1, lb, ub);
  lub_constrain_tests::expect(x2, lb, ub);

  // both inf
  lub_constrain_tests::expect(x1, lb, ub_inf);
  lub_constrain_tests::expect(x2, lb, ub_inf);

}

TEST(mathMixMatFun, lub_constrain_vector_scalar_vector) {
  Eigen::MatrixXd x1(2, 2);
  x1 << 5.0, 2.0, 4.0, 5.0;
  Eigen::MatrixXd x2(2, 2);
  x2 << -1.1, 0.005, 1.0, 3.0;
  double lb = -2.0;
  Eigen::MatrixXd ub(2, 2);
  ub << -1.0, 5.0, 0.0, 38.0;

  lub_constrain_tests::expect(x1, lb, ub);
  lub_constrain_tests::expect(x2, lb, ub);

  // ub inf
  ub(1, 1) = stan::math::INFTY;
  lub_constrain_tests::expect(x1, lb, ub);
  lub_constrain_tests::expect(x2, lb, ub);


  // lb inf
  ub(1, 1) = 38.0;
  auto lb_inf = stan::math::NEGATIVE_INFTY;
  lub_constrain_tests::expect(x1, lb_inf, ub);
  lub_constrain_tests::expect(x2, lb_inf, ub);

  // both inf
  ub(1, 1) = stan::math::INFTY;
  lub_constrain_tests::expect(x1, lb_inf, ub);
  lub_constrain_tests::expect(x2, lb_inf, ub);

}

TEST(mathMixMatFun, lub_constrain_vector_vector_vector) {

  Eigen::MatrixXd x1(2, 2);
  x1 << 5.0, 2.0, 4.0, 5.0;
  Eigen::MatrixXd x2(2, 2);
  x2 << -1.1, 0.005, 1.0, 3.0;
  Eigen::MatrixXd lb(2, 2);
  lb << -3.0, 3.0, -6.0, 6.0;
  Eigen::MatrixXd ub(2, 2);
  ub << -1.0, 5.0, 0.0, 38.0;
  lub_constrain_tests::expect(x1, lb, ub);
  lub_constrain_tests::expect(x2, lb, ub);

  // ub inf
  ub(1, 1) = stan::math::INFTY;
  lub_constrain_tests::expect(x1, lb, ub);
  lub_constrain_tests::expect(x2, lb, ub);


  // lb inf
  ub(1, 1) = 38.0;
  lb(1, 1) = stan::math::NEGATIVE_INFTY;
  lub_constrain_tests::expect(x1, lb, ub);
  lub_constrain_tests::expect(x2, lb, ub);

  // both inf
  ub(1, 1) = stan::math::INFTY;
  lub_constrain_tests::expect(x1, lb, ub);
  lub_constrain_tests::expect(x2, lb, ub);

}


// real[], real[], real[]
// real[], real, real[]
// real[], real[], real
TEST(mathMixMatFun, lub_stdvec_constrain) {
  std::vector<double> A{5.0, 2.0, 4.0, -2.0};
  std::vector<double> lbm{-3.0, 3.0, -6.0, 6.0};
  std::vector<double> ubm{-1.0, 5.0, 0.0, 38.0};
  lub_constrain_tests::expect_vec(A, lbm, ubm);
  double lbd = -6.0;
  lub_constrain_tests::expect_vec(A, lbd, ubm);
  double ubd = 8.0;
  lub_constrain_tests::expect_vec(A, lbd, ubd);
  lub_constrain_tests::expect_vec(A, lbm, ubd);
}

TEST(mathMixMatFun, lub_stdvec_constrain_neg_inf) {
  std::vector<double> A{5.0, 2.0, 4.0, -2.0};
  std::vector<double> lbm{stan::math::NEGATIVE_INFTY, 3.0, stan::math::NEGATIVE_INFTY, 6.0};
  std::vector<double> ubm{-1.0, stan::math::INFTY, stan::math::INFTY, 38.0};
  lub_constrain_tests::expect_vec(A, lbm, ubm);
  double lbd = stan::math::NEGATIVE_INFTY;
  lub_constrain_tests::expect_vec(A, lbd, ubm);
  double ubd = stan::math::INFTY;
  lub_constrain_tests::expect_vec(A, lbd, ubd);
  lub_constrain_tests::expect_vec(A, lbm, ubd);
}

// array matrix[], array matrix[], array matrix[]
// array matrix[], array matrix[], matrix[]
// array matrix[], matrix[], array matrix[]
// array matrix[], matrix[], matrix[]
// array matrix[], array matrix[], real
// array matrix[], real, array matrix[]
// array matrix[], matrix[], real
// array matrix[], real, matrix[]
// array matrix[], real, real
TEST(mathMixMatFun, lub_stdvec_mat_mat_constrain) {
  Eigen::MatrixXd A_inner(2, 3);
  // swapping 0.05 for 0 causes a failure for the hessian?
  A_inner << 5.0, 2.0, 4.0, -2.0, 0.05, 0.1;
  Eigen::MatrixXd lb_inner(2, 3);
  lb_inner << -1.0, 1.0, -6.0, 1.0, 0.0, 0.01;
  Eigen::MatrixXd ub_inner(2, 3);
  ub_inner << 6.0, 3.0, 12.0, 38.0, 0.1, 0.15;

  std::vector<Eigen::MatrixXd> A{A_inner, A_inner};
  std::vector<Eigen::MatrixXd> lb_vec{lb_inner, lb_inner};
  std::vector<Eigen::MatrixXd> ub_vec{ub_inner, ub_inner};
  lub_constrain_tests::expect_vec(A, lb_vec, ub_vec);
  lub_constrain_tests::expect_vec(A, lb_vec, ub_inner);
  lub_constrain_tests::expect_vec(A, lb_inner, ub_vec);
  lub_constrain_tests::expect_vec(A, lb_inner, ub_inner);
  double lb_scal = -1.0;
  double ub_scal = 7.0;
  lub_constrain_tests::expect_vec(A, lb_vec, ub_scal);
  lub_constrain_tests::expect_vec(A, lb_scal, ub_vec);
  lub_constrain_tests::expect_vec(A, lb_scal, ub_inner);
  lub_constrain_tests::expect_vec(A, lb_inner, ub_scal);
  lub_constrain_tests::expect_vec(A, lb_scal, ub_scal);
}

TEST(mathMixMatFun, lub_stdvec_mat_mat_constrain_infty) {
  Eigen::MatrixXd A_inner(2, 3);
  A_inner << 5.0, 2.0, 4.0, -2.0, 0.05, 0.1;
  Eigen::MatrixXd lb_inner(2, 3);
  lb_inner << stan::math::NEGATIVE_INFTY, 1.0, stan::math::NEGATIVE_INFTY, 1.0, 0.0, 0.01;
  Eigen::MatrixXd ub_inner(2, 3);
  ub_inner << 6.0, stan::math::INFTY, stan::math::INFTY, 38.0, 0.1, 0.15;

  std::vector<Eigen::MatrixXd> A{A_inner, A_inner};
  std::vector<Eigen::MatrixXd> lb_vec{lb_inner, lb_inner};
  std::vector<Eigen::MatrixXd> ub_vec{ub_inner, ub_inner};
  lub_constrain_tests::expect_vec(A, lb_vec, ub_vec);
  lub_constrain_tests::expect_vec(A, lb_vec, ub_inner);
  lub_constrain_tests::expect_vec(A, lb_inner, ub_vec);
  lub_constrain_tests::expect_vec(A, lb_inner, ub_inner);
  double lb_scal = stan::math::NEGATIVE_INFTY;
  double ub_scal = stan::math::INFTY;
  lub_constrain_tests::expect_vec(A, lb_vec, ub_scal);
  lub_constrain_tests::expect_vec(A, lb_scal, ub_vec);
  lub_constrain_tests::expect_vec(A, lb_scal, ub_inner);
  lub_constrain_tests::expect_vec(A, lb_inner, ub_scal);
  lub_constrain_tests::expect_vec(A, lb_scal, ub_scal);
}
