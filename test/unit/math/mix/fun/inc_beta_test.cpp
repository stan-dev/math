#include <test/unit/math/test_ad.hpp>
#include <limits>

TEST(mathMixScalFun, inc_beta) {
  auto f = [](const auto& x1, const auto& x2, const auto& x3) {
    return stan::math::inc_beta(x1, x2, x3);
  };

  // replicate pre-existing test cases
  stan::test::expect_ad(f, 0.6, 0.3, 0.5);
  stan::test::expect_ad(f, 0.9, 0.9, 0.9);
  stan::test::expect_ad(f, 1.0, 1.0, 0.4);

  // integer instantiations
  stan::test::expect_ad(f, 3, 2, 0.2);
  stan::test::expect_ad(f, 3, 2.0, 0.2);
  stan::test::expect_ad(f, 3.0, 2, 0.2);
  stan::test::expect_ad(f, 3.0, 2.0, 0.2);

  // test all nan instantiations
  double nan = std::numeric_limits<double>::quiet_NaN();
  stan::test::expect_ad(f, 0.6, 0.3, nan);
  stan::test::expect_ad(f, 0.6, nan, 0.5);
  stan::test::expect_ad(f, 0.6, nan, nan);
  stan::test::expect_ad(f, nan, 0.3, 0.5);
  stan::test::expect_ad(f, nan, 0.3, nan);
  stan::test::expect_ad(f, nan, nan, 0.5);
  stan::test::expect_ad(f, nan, nan, nan);
}

TEST(mathMixScalFun, inc_beta_upper_bound_infinite_density) {
  // At z = 1 with b < 1 the density is infinite; the forward and reverse
  // sweeps of fvar<var> must not throw there.
  using stan::math::fvar;
  using stan::math::var;
  var a_v = 1.05;
  var b_v = 0.35;
  fvar<var> a(a_v, 1.0);
  fvar<var> b(b_v, 0.0);
  fvar<var> z(1.0, 0.0);
  fvar<var> f;
  EXPECT_NO_THROW(f = stan::math::inc_beta(a, b, z));
  EXPECT_EQ(1.0, f.val_.val());
  EXPECT_NO_THROW(f.val_.grad());
  stan::math::set_zero_all_adjoints();
  EXPECT_NO_THROW(f.d_.grad());
  stan::math::recover_memory();
}

TEST(mathMixScalFun, inc_beta_vec) {
  auto f = [](const auto& x1, const auto& x2, const auto& x3) {
    return stan::math::inc_beta(x1, x2, x3);
  };

  Eigen::VectorXd in1(2);
  in1 << 3, 1;
  Eigen::VectorXd in2(2);
  in2 << 10.2, 3.4;
  Eigen::VectorXd in3(2);
  in3 << 0.1, 0.4;
  stan::test::expect_ad_vectorized_ternary(f, in1, in2, in3);
}
