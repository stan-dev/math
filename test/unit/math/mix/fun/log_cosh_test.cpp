#include <test/unit/math/test_ad.hpp>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

// Every autodiff mode of log_cosh must agree with finite differences on
// both sides of the |x| = 1 branch.
TEST(mathMixMatFun, logCosh) {
  auto f = [](auto&& x1) { return stan::math::log_cosh(x1); };
  stan::test::expect_common_unary_vectorized<stan::test::ScalarSupport::Real>(
      f);
  stan::test::expect_unary_vectorized<stan::test::ScalarSupport::Real>(
      f, -30.0, -2.6, -1.0, -0.5, 0.0, 0.5, 0.999, 1.0, 1.001, 2.6, 30.0);
}

// The var_value overload has its own reverse pass; it must match the
// Eigen<var> path.
TEST(mathMixMatFun, logCoshVarMat) {
  using stan::math::vec_concat;
  using stan::test::expect_ad_vector_matvar;
  using stan::test::internal::common_nonzero_args;
  auto f = [](auto&& x1) { return stan::math::log_cosh(x1); };
  std::vector<double> com_args = common_nonzero_args();
  std::vector<double> args{-30.0, -2.6, -0.5, 0.999, 1.001, 2.6, 30.0};
  auto all_args = vec_concat(com_args, args);
  Eigen::VectorXd A(all_args.size());
  for (int i = 0; i < all_args.size(); ++i) {
    A(i) = all_args[i];
  }
  expect_ad_vector_matvar(f, A);
}

// log(cosh(x)) rounds to 0 near x = 0 and overflows at 800; log_cosh must
// match mpmath to a few eps there and across the |x| = 1 branch.
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
  for (auto&& c : cases) {
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

// The var_value Cholesky factor takes tanh and sech from this helper; all
// three outputs must stay accurate through the underflow of sech at 800.
TEST(mathMixMatFun, logCoshTanhSech) {
  const double eps = std::numeric_limits<double>::epsilon();
  // x, log(cosh(x)), tanh(x), sech(x)
  std::vector<std::array<double, 4>> cases
      = {{0.0, 0.0, 0.0, 1.0},
         {1e-10, 5.0e-21, 1.0e-10, 1.0},
         {0.5, 0.12011450695827752, 0.46211715726000976, 0.88681888397007391},
         {0.999, 0.43301944642088479, 0.76117386166056995, 0.64854787974731328},
         {1.0, 0.43378083048302719, 0.76159415595576489, 0.6480542736638854},
         {1.001, 0.43454263451956291, 0.76201381055100676, 0.64756077130222642},
         {-3.0, 2.3093285045777851, -0.99505475368673045, 0.099327927419433208},
         {30.0, 29.306852819440055, 1.0, 1.8715245937680349e-13},
         {700.0, 699.30685281944005, 1.0, 1.9719353087519542e-304},
         {800.0, 799.30685281944005, 1.0, 0.0}};
  for (auto&& c : cases) {
    const double x = c[0];
    const stan::math::internal::log_cosh_parts parts
        = stan::math::internal::log_cosh_tanh_sech(x);
    EXPECT_NEAR(c[1], parts.log_cosh, 4 * eps * std::fabs(c[1])) << "x = " << x;
    EXPECT_NEAR(c[2], parts.tanh, 4 * eps * std::fabs(c[2])) << "x = " << x;
    EXPECT_NEAR(c[3], parts.sech, 4 * eps * std::fabs(c[3])) << "x = " << x;
  }
}
