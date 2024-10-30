#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/constraint/offset_multiplier_constrain_matvar_helpers.hpp>

// real[], real[], real[]
// real[], real, real[]
// real[], real[], real
TEST(mathMixMatFun, offset_multiplier_constrain_matvar_stdvec_constrain) {
  std::vector<double> A{5.0, 2.0, 4.0, -2.0};
  std::vector<double> mum{-3.0, 3.0, -6.0, 6.0};
  std::vector<double> sigmam{-1.0, 5.0, 0.0, 38.0};
  offset_multiplier_constrain_tests::expect_vec_matvar(A, mum, sigmam);
  double mud = -6.0;
  offset_multiplier_constrain_tests::expect_vec_matvar(A, mud, sigmam);
  double sigmad = 8.0;
  offset_multiplier_constrain_tests::expect_vec_matvar(A, mud, sigmad);
  offset_multiplier_constrain_tests::expect_vec_matvar(A, mum, sigmad);
}
