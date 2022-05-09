#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, invPhiLog) {
  auto f = [](const auto& x1) { return stan::math::inv_Phi_log(x1); };
  //  stan::test::expect_common_unary_vectorized(f);
  stan::test::expect_unary_vectorized(f, log(0.02425),
                                      log(0.97575));  // breakpoints
  stan::test::expect_unary_vectorized(f, -100.25, -2, 0.01, 0.1, 0.98, 0.5,
                                      2.0);
}

TEST(mathMixScalLogFun, invPhiLog) {
  auto f = [](const auto& x1) { return stan::math::inv_Phi_log(x1); };
  // stan::test::expect_common_nonzero_binary(f);
  stan::test::expect_ad(f, -1.1);
}

TEST(mathMixMatFunLog, invPhiLog) {
  auto f = [](const auto& x1) { return stan::math::inv_Phi_log(x1); };
  // stan::test::expect_common_unary_vectorized(f);
  stan::test::expect_unary_vectorized(f, log(0.02425),
                                      log(0.97575));  // breakpoints
  stan::test::expect_unary_vectorized(f, -100.25, -2, 0.01, 0.1, 0.98, 0.5,
                                      2.0);
}

// TEST(mathMixMatFunLog, invPhiLog_varmat) {
//   using stan::math::vec_concat;
//   using stan::test::expect_ad_vector_matvar;
//   using stan::test::internal::common_args;
//   auto f = [](const auto& x1) {
//     return stan::math::inv_Phi_log(x1);
//   };
//   std::vector<double> com_args = common_args();
//   std::vector<double> args{log(0.02425), log(0.97575), 2, log(0.1), 0.5,
//   log(2.0)}; auto all_args = vec_concat(com_args, args); Eigen::VectorXd
//   A(all_args.size()); for (int i = 0; i < all_args.size(); ++i) {
//     A(i) = all_args[i];
//   }
//   expect_ad_vector_matvar(f, A);
// }
