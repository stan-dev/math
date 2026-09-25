#include <stan/math/mix.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/fun/util.hpp>

TEST(ProbInternalMath, gradRegIncGamma_typical) {
  double a = 0.5;
  double b = 1.0;
  double g = 1.77245;
  double dig = -1.96351;

  // d/da Q(0.5, 1.0) from mpmath at 80 digits. The g and dig arguments are
  // accepted for signature compatibility and do not affect the result.
  EXPECT_FLOAT_EQ(0.38983726432851057,
                  stan::math::grad_reg_inc_gamma(a, b, g, dig));
}

TEST(ProbInternalMath, gradRegIncGamma_infLoopInVersion2_0_1) {
  double a = 8.01006;
  double b = 2.47579e+215;
  double g = 5143.28;
  double dig = 2.01698;

  EXPECT_FLOAT_EQ(0, stan::math::grad_reg_inc_gamma(a, b, g, dig));
}

TEST(ProbInternalMath, gradRegIncGamma_largeZ) {
  double a = 3;
  double z = 48;
  double g = 2.0;
  double dig = 0.9227843;

  EXPECT_FLOAT_EQ(5.08294581508e-18,
                  stan::math::grad_reg_inc_gamma(a, z, g, dig));
}

TEST(ProbInternalMath, gradRegIncGamma_fd) {
  using stan::math::fvar;

  fvar<double> a = 0.5;
  fvar<double> b = 1.0;
  fvar<double> g = 1.77245;
  fvar<double> dig = -1.96351;

  EXPECT_FLOAT_EQ(0.38983726432851057,
                  stan::math::grad_reg_inc_gamma(a, b, g, dig).val());
}
TEST(ProbInternalMath, gradRegIncGamma_ffd) {
  using stan::math::fvar;

  fvar<fvar<double> > a = 0.5;
  fvar<fvar<double> > b = 1.0;
  fvar<fvar<double> > g = 1.77245;
  fvar<fvar<double> > dig = -1.96351;

  EXPECT_FLOAT_EQ(0.38983726432851057,
                  stan::math::grad_reg_inc_gamma(a, b, g, dig).val_.val_);
}

TEST(ProbInternalMath, gradRegIncGamma_fv) {
  using stan::math::digamma;
  using stan::math::fvar;
  using stan::math::var;

  fvar<var> a = 0.5;
  fvar<var> b = 1.0;
  fvar<var> g = 1.77245;
  fvar<var> dig = digamma(a);

  EXPECT_FLOAT_EQ(0.38983726432851057,
                  stan::math::grad_reg_inc_gamma(a, b, g, dig).val_.val());
}

TEST(ProbInternalMath, gradRegIncGamma_fv_1stderiv) {
  using stan::math::digamma;
  using stan::math::fvar;
  using stan::math::tgamma;
  using stan::math::var;

  fvar<var> a = 0.5;
  a.d_ = 1.0;
  fvar<var> b = 1.0;
  fvar<var> g = tgamma(a);
  fvar<var> dig = digamma(a);
  a.d_ = 1.0;

  fvar<var> z = stan::math::grad_reg_inc_gamma(a, b, g, dig);

  std::vector<stan::math::var> y1{a.val_};
  std::vector<double> grad1;
  z.val_.grad(y1, grad1);
  EXPECT_NEAR(0.2134999674954450667, grad1[0], 1e-6);
}

TEST(ProbInternalMath, gradRegIncGamma_fv_2ndderiv) {
  using stan::math::digamma;
  using stan::math::fvar;
  using stan::math::tgamma;
  using stan::math::var;

  fvar<var> a = 0.5;
  a.d_ = 1.0;
  fvar<var> b = 1.0;
  fvar<var> g = tgamma(a);
  fvar<var> dig = digamma(a);
  a.d_ = 1.0;

  fvar<var> z = stan::math::grad_reg_inc_gamma(a, b, g, dig);

  std::vector<stan::math::var> y1{a.val_};
  std::vector<double> grad1;
  z.d_.grad(y1, grad1);
  EXPECT_NEAR(-0.546236927878295422, grad1[0], 1e-6);
}

namespace {
struct gamma_cdf_shape_functor {
  double y_;
  explicit gamma_cdf_shape_functor(double y) : y_(y) {}
  template <typename T>
  T operator()(const Eigen::Matrix<T, Eigen::Dynamic, 1>& x) const {
    return stan::math::gamma_cdf(y_, x(0), 1.0);
  }
};
}  // namespace

/**
 * The first derivative must not depend on how it is obtained. Reverse mode
 * reaches the gradient root with `double` partials and `hessian` reaches
 * it with `fvar<var>`; the two must agree. This needs no external
 * reference value.
 */
TEST(ProbInternalMath, gradRegIncGamma_gradient_matches_hessian) {
  using stan::math::var;

  for (double alpha : {5.0, 20.0, 100.0, 171.0, 180.0, 300.0}) {
    const double y = alpha;

    double rev_grad;
    {
      stan::math::nested_rev_autodiff nested;
      var a = alpha;
      var p = stan::math::gamma_cdf(y, a, 1.0);
      p.grad();
      rev_grad = a.adj();
    }

    Eigen::Matrix<double, Eigen::Dynamic, 1> x(1);
    x(0) = alpha;
    double fx;
    Eigen::Matrix<double, Eigen::Dynamic, 1> grad(1);
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> hess(1, 1);
    stan::math::hessian(gamma_cdf_shape_functor(y), x, fx, grad, hess);

    ASSERT_TRUE(std::isfinite(rev_grad)) << "reverse mode, alpha = " << alpha;
    ASSERT_TRUE(std::isfinite(grad(0))) << "hessian(), alpha = " << alpha;
    EXPECT_NEAR(rev_grad, grad(0), 1e-12 * std::fabs(rev_grad))
        << "alpha = " << alpha;
  }
}
