#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/constraint/lub_constrain_helpers.hpp>

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
  std::vector<double> lbm{stan::math::NEGATIVE_INFTY, 3.0,
                          stan::math::NEGATIVE_INFTY, 6.0};
  std::vector<double> ubm{-1.0, stan::math::INFTY, stan::math::INFTY, 38.0};
  lub_constrain_tests::expect_vec(A, lbm, ubm);
  double lbd = stan::math::NEGATIVE_INFTY;
  lub_constrain_tests::expect_vec(A, lbd, ubm);
  double ubd = stan::math::INFTY;
  lub_constrain_tests::expect_vec(A, lbd, ubd);
  lub_constrain_tests::expect_vec(A, lbm, ubd);
}
