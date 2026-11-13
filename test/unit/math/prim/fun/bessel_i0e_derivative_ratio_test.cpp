#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathFunctions, bessel_i0e_derivative_ratio) {
  using stan::math::bessel_i0e_derivative_ratio;
  // n, z and rho_n(z), computed with mpmath (60 digits); rho_0 = 1 - I1 / I0.
  // The points include the switches of its evaluation, at 7.75 and 25.
  const double ref[][3] = {{0, 0, 1.0},
                           {0, 1e-8, 0.999999995},
                           {0, 0.5, 0.75750038741919805},
                           {0, 2, 0.30222534203599202},
                           {0, 7.7, 0.067397060996946867},
                           {0, 7.75, 0.066943433087775686},
                           {0, 15, 0.033930436013491875},
                           {0, 21, 0.024107586113781809},
                           {0, 24.5, 0.020625503352569103},
                           {0, 25, 0.020208546509484069},
                           {0, 30, 0.016810444634663907},
                           {0, 100, 0.0050126269948312344},
                           {0, 1e4, 5.0001250125019535e-5},
                           {0, 1e12, 5.00000000000125e-13},
                           {0, 1e300, 5.0e-301},
                           {1, 0, 1.5},
                           {1, 0.5, 1.3597373239451467},
                           {1, 15, 0.10186236411783493},
                           {1, 1e4, 0.00015000375056261956},
                           {1, 1e12, 1.500000000000375e-12},
                           {3, 0, 1.75},
                           {3, 0.5, 1.7003232014160718},
                           {3, 15, 0.23805719735624314},
                           {3, 1e4, 0.00035000875218819479},
                           {3, 1e12, 3.500000000000875e-12}};
  for (const auto& r : ref) {
    // below 25, rho_0 comes from approximations of I0 and I1
    const double tol = r[0] == 0 && r[1] < 25 ? 1e-13 : 1e-14;
    EXPECT_NEAR(r[2], bessel_i0e_derivative_ratio(static_cast<int>(r[0]), r[1]),
                tol * std::fabs(r[2]));
  }
  EXPECT_NEAR(2.941176470588235e-309, bessel_i0e_derivative_ratio(0, 1.7e308),
              1e-322);
  EXPECT_EQ(0.0, bessel_i0e_derivative_ratio(
                     0, std::numeric_limits<double>::infinity()));
  // integer argument
  EXPECT_NEAR(0.30222534203599202, bessel_i0e_derivative_ratio(0, 2), 1e-14);
  // the largest order; reference from a backward recurrence in mpmath
  EXPECT_NEAR(1.9836288172292427, bessel_i0e_derivative_ratio(1000, 480.0),
              2e-13);

  Eigen::VectorXd z(2);
  z << 0.5, 1e12;
  Eigen::VectorXd v = bessel_i0e_derivative_ratio(0, z);
  EXPECT_NEAR(0.75750038741919805, v(0), 1e-14);
  EXPECT_NEAR(5.00000000000125e-13, v(1), 1e-26);

  EXPECT_THROW(bessel_i0e_derivative_ratio(0, -1.0), std::domain_error);
  EXPECT_THROW(
      bessel_i0e_derivative_ratio(0, std::numeric_limits<double>::quiet_NaN()),
      std::domain_error);
  EXPECT_THROW(bessel_i0e_derivative_ratio(-1, 1.0), std::domain_error);
  EXPECT_THROW(bessel_i0e_derivative_ratio(1001, 1.0), std::domain_error);
}
