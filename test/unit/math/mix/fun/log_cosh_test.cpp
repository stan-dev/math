#include <test/unit/math/test_ad.hpp>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

TEST(mathMixMatFun, logCosh) {
  auto f = [](const auto& x1) { return stan::math::log_cosh(x1); };
  stan::test::expect_common_unary_vectorized<stan::test::ScalarSupport::Real>(
      f);
  stan::test::expect_unary_vectorized<stan::test::ScalarSupport::Real>(
      f, -30.0, -2.6, -1.0, -0.5, 0.0, 0.5, 0.999, 1.0, 1.001, 2.6, 30.0);
}

TEST(mathMixMatFun, logCoshVarMat) {
  using stan::math::vec_concat;
  using stan::test::expect_ad_vector_matvar;
  using stan::test::internal::common_nonzero_args;
  auto f = [](const auto& x1) { return stan::math::log_cosh(x1); };
  std::vector<double> com_args = common_nonzero_args();
  std::vector<double> args{-30.0, -2.6, -0.5, 0.999, 1.001, 2.6, 30.0};
  auto all_args = vec_concat(com_args, args);
  Eigen::VectorXd A(all_args.size());
  for (int i = 0; i < all_args.size(); ++i) {
    A(i) = all_args[i];
  }
  expect_ad_vector_matvar(f, A);
}

// Values and derivatives against mpmath, on both sides of the branch at
// |x| = 1, at a tiny argument (where log(cosh(x)) = x^2 / 2), and where
// cosh(x) overflows.
TEST(mathMixMatFun, logCoshReference) {
  using stan::math::fvar;
  using stan::math::var;
  // x, log(cosh(x)), tanh(x)
  std::vector<std::array<double, 3>> cases
      = {{1e-10, 5.0e-21, 1.0e-10},
         {0.5, 0.12011450695827752, 0.46211715726000976},
         {0.999, 0.43301944642088479, 0.76117386166056995},
         {1.0, 0.43378083048302719, 0.76159415595576489},
         {-3.0, 2.3093285045777851, -0.99505475368673045},
         {30.0, 29.306852819440055, 1.0},
         {800.0, 799.30685281944005, 1.0}};
  for (const auto& c : cases) {
    const double x = c[0];
    const double tol_v = 4 * std::numeric_limits<double>::epsilon() * c[1];
    const double tol_d
        = 4 * std::numeric_limits<double>::epsilon() * std::fabs(c[2]);
    EXPECT_NEAR(c[1], stan::math::log_cosh(x), tol_v) << "x = " << x;

    var x_v = x;
    var y_v = stan::math::log_cosh(x_v);
    y_v.grad();
    EXPECT_NEAR(c[1], y_v.val(), tol_v) << "x = " << x;
    EXPECT_NEAR(c[2], x_v.adj(), tol_d) << "x = " << x;
    stan::math::recover_memory();

    fvar<double> y_f = stan::math::log_cosh(fvar<double>(x, 1.0));
    EXPECT_NEAR(c[1], y_f.val_, tol_v) << "x = " << x;
    EXPECT_NEAR(c[2], y_f.d_, tol_d) << "x = " << x;
  }
}
