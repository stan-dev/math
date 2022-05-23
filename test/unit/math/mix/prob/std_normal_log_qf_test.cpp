#include <stan/math/rev.hpp>
#include <stan/math/prim.hpp>
#include <test/unit/math/rev/util.hpp>
#include <test/unit/math/test_ad.hpp>
#include <stan/math/fwd/prob/std_normal_log_qf.hpp>

TEST(mathMixLogFun, stdNormalLogQf) {
  auto f = [](const auto& x1) { return stan::math::std_normal_log_qf(x1); };
  stan::test::expect_ad(f, -100.25);
  stan::test::expect_unary_vectorized(
      f, -100.25, -2, 0.01, 0.1, 0.98, 0.5, 2.0, -1.3, 0.49, 0.99, 1.01,
      stan::math::not_a_number(), stan::math::positive_infinity(),
      stan::math::negative_infinity());
  stan::test::expect_unary_vectorized(f, log(0.02425), log(0.97575));
}

TEST(mathMixScalLogFun, stdNormalLogQfInt) {
  auto f = [](const auto& x1) { return stan::math::std_normal_log_qf(x1); };
  int y = 1;
  stan::test::expect_ad(f, y);
  y = -1;
  stan::test::expect_ad(f, y);
}

TEST(mathMixZeroLogFun, stdNormalLogQfZero) {
  auto f = [](const auto& x1) { return stan::math::std_normal_log_qf(x1); };
  int y_int = 0;
  stan::test::expect_ad(f, y_int);

  double y = 0;
  stan::test::expect_ad(f, y);
}

TEST(mathMixMatFunLog, stdNormalLogQfVarmat) {
  using stan::math::vec_concat;
  using stan::test::expect_ad_vector_matvar;
  using stan::test::internal::common_args;
  auto f = [](const auto& x1) { return stan::math::std_normal_log_qf(x1); };
  std::vector<double> com_args = common_args();
  std::vector<double> args{log(0.02425), log(0.97575), 2,
                           log(0.1),     0.5,          log(2.0)};
  auto all_args = vec_concat(com_args, args);
  Eigen::VectorXd A(all_args.size());
  for (int i = 0; i < all_args.size(); ++i) {
    A(i) = all_args[i];
  }
  expect_ad_vector_matvar(f, A);
}
