#include <stan/math/mix.hpp>
#include <test/unit/math/rev/util.hpp>
#include <gtest/gtest.h>
#include <array>
#include <cmath>
#include <vector>

// First and second derivatives of logistic_lpdf with respect to y in nested
// forward mode and in fvar<var>, against mpmath:
// d/dy = -tanh(z / 2) / sigma and d2/dy2 = -sech(z / 2)^2 / (2 sigma^2)
// with z = (y - mu) / sigma. The second derivative is the derivative of
// tanh(z / 2), which must not round to 0 for large |z|.
TEST_F(AgradRev, logistic_lpdf_second_derivative_y) {
  using stan::math::fvar;
  using stan::math::var;
  // y, mu, sigma, d/dy, d2/dy2
  std::vector<std::array<double, 5>> cases
      = {{1.0, 0.0, 1.0, -0.46211715726000976, -0.39322386648296371},
         {-20.0, 0.0, 1.0, 0.99999999587769276, -4.1223072278836987e-9},
         {40.0, 0.0, 1.0, -0.99999999999999999, -8.4967085105831779e-18},
         {100.0, 0.0, 1.0, -1.0, -7.4401519520416719e-44},
         {80.5, 0.5, 2.0, -0.5, -2.1241771276457945e-18}};
  for (const auto& c : cases) {
    const double y = c[0];
    const double mu = c[1];
    const double sigma = c[2];

    fvar<fvar<double>> y_ffd(fvar<double>(y, 1.0), fvar<double>(1.0, 0.0));
    fvar<fvar<double>> lp_ffd = stan::math::logistic_lpdf(y_ffd, mu, sigma);
    EXPECT_NEAR(c[3], lp_ffd.d_.val_, 1e-14 * std::fabs(c[3])) << "y = " << y;
    EXPECT_NEAR(c[4], lp_ffd.d_.d_, 1e-14 * std::fabs(c[4])) << "y = " << y;

    var y_v = y;
    fvar<var> y_fv(y_v, 1.0);
    fvar<var> lp_fv = stan::math::logistic_lpdf(y_fv, mu, sigma);
    lp_fv.d_.grad();
    EXPECT_NEAR(c[3], lp_fv.d_.val(), 1e-14 * std::fabs(c[3])) << "y = " << y;
    EXPECT_NEAR(c[4], y_v.adj(), 1e-14 * std::fabs(c[4])) << "y = " << y;
    stan::math::recover_memory();
  }
}

// The same through the vectorized tanh of an Eigen argument. With the
// direction (1, 1) the second directional derivative is the sum of the two
// second derivatives.
TEST_F(AgradRev, logistic_lpdf_second_derivative_y_vector) {
  using stan::math::fvar;
  Eigen::Matrix<fvar<fvar<double>>, Eigen::Dynamic, 1> y(2);
  y << fvar<fvar<double>>(fvar<double>(40.0, 1.0), fvar<double>(1.0, 0.0)),
      fvar<fvar<double>>(fvar<double>(100.0, 1.0), fvar<double>(1.0, 0.0));
  fvar<fvar<double>> lp = stan::math::logistic_lpdf(y, 0.0, 1.0);
  const double d2 = -8.4967085105831779e-18 + -7.4401519520416719e-44;
  EXPECT_NEAR(d2, lp.d_.d_, 1e-14 * std::fabs(d2));
}
