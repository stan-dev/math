#include <test/unit/math/test_ad.hpp>
#include <cmath>
#include <limits>
#include <vector>

TEST(mathMixMatFun, log1pExp) {
  auto f = [](const auto& x1) { return stan::math::log1p_exp(x1); };
  stan::test::expect_common_nonzero_unary_vectorized<
      stan::test::ScalarSupport::Real>(f);
  stan::test::expect_unary_vectorized(f, -2.6, -2, -1, -0.5, -0.2, 0.5, 1.0,
                                      1.3, 2, 3);

  std::vector<double> com_args = stan::test::internal::common_nonzero_args();
  std::vector<double> args{0.1, -2.5, 5.5};

  stan::test::expect_ad_vector_matvar(f, stan::math::to_vector(com_args));
  stan::test::expect_ad_vector_matvar(f, stan::math::to_vector(args));
}

TEST(mathMixMatFun, log1pExpHigherDerivativesLargeArgs) {
  using stan::math::fvar;
  using stan::math::log1p_exp;
  using stan::math::var;
  // {x, log1p_exp''(x), log1p_exp'''(x)}: mpmath, 50 digits; 0 where the
  // value underflows (e^-800)
  std::vector<std::vector<double>> points{
      {-800, 0, 0},
      {-40, 4.2483542552915889592e-18, 4.2483542552915889231e-18},
      {-20, 2.0611536139418493437e-9, 2.0611536054451408856e-9},
      {-1, 1.9661193324148185254e-1, 9.0857747672948409442e-2},
      {0, 0.25, 0},
      {1, 1.9661193324148185254e-1, -9.0857747672948409442e-2},
      {20, 2.0611536139418493437e-9, -2.0611536054451408856e-9},
      {30, 9.3576229688384233028e-14, -9.3576229688366720006e-14},
      {37, 8.533047625744064338e-17, -8.5330476257440628818e-17},
      {40, 4.2483542552915889592e-18, -4.2483542552915889231e-18},
      {100, 3.720075976020835963e-44, -3.720075976020835963e-44},
      {700, 9.8596765437597708567e-305, -9.8596765437597708567e-305},
      {800, 0, 0}};
  for (const auto& p : points) {
    double x = p[0];
    double d2 = p[1];
    double d3 = p[2];
    double tol2 = 1e-13 * std::fmax(d2, std::numeric_limits<double>::min());
    double tol3
        = 1e-13 * std::fmax(std::fabs(d3), std::numeric_limits<double>::min());

    fvar<fvar<double>> ff(fvar<double>(x, 1), fvar<double>(1, 0));
    EXPECT_NEAR(log1p_exp(ff).d_.d_, d2, tol2) << "x = " << x;

    fvar<var> fv(x, 1);
    log1p_exp(fv).d_.grad();
    EXPECT_NEAR(fv.val_.adj(), d2, tol2) << "x = " << x;
    stan::math::recover_memory();

    fvar<fvar<var>> ffv(fvar<var>(x, 1), fvar<var>(1, 0));
    log1p_exp(ffv).d_.d_.grad();
    EXPECT_NEAR(ffv.val_.val_.adj(), d3, tol3) << "x = " << x;
    stan::math::recover_memory();
  }
}
