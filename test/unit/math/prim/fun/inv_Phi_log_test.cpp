#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(MathFunctionsLog, inv_Phi_log) {
  using stan::math::inv_Phi_log;
  using stan::math::Phi;
  EXPECT_FLOAT_EQ(0.0, inv_Phi_log(log(0.5)));
  double log_p = log(0.123456789);
  EXPECT_FLOAT_EQ(0.123456789, Phi(inv_Phi_log(log_p)));
  log_p = log(8e-311);
  EXPECT_FLOAT_EQ(8e-311, Phi(inv_Phi_log(log_p)));
  log_p = log(0.99);
  EXPECT_FLOAT_EQ(0.99, Phi(inv_Phi_log(log_p)));

  // breakpoints
  log_p = log(0.02425);
  EXPECT_FLOAT_EQ(0.02425, Phi(inv_Phi_log(log_p)));
  log_p = log(0.97575);
  EXPECT_FLOAT_EQ(0.97575, Phi(inv_Phi_log(log_p)));
}

TEST(MathFunctionsLog, Equal) {
  using stan::math::inv_Phi_log;
  // test output generated with WolframAlpha

  double log_p[]
      = {-20.72326583694641115, -16.11809565095831978, -11.51292546497022842,
         -6.907755278982137052, -2.995732273553990993, -1.897119984885881302,
         -1.386294361119890618, -1.049822124498677688, -0.798507696217771610,
         -0.597837000755620449, -0.430782916092454257, -0.287682072451780927,
         -0.162518929497774913, -0.051293294387550533, -0.001000500333583533,
         -0.000010000050000333, -1.000000050000003e-7, -1.000000000500000e-9};

  double exact[]
      = {-5.997807015007686871, -5.199337582192816931, -4.264890793922824628,
         -3.090232306167813541, -1.644853626951472714, -1.036433389493789579,
         -0.674489750196081743, -0.385320466407567623, -0.125661346855074034,
         0.1256613468550740342, 0.3853204664075676238, 0.6744897501960817432,
         1.0364333894937895797, 1.6448536269514727148, 3.0902323061678135415,
         4.2648907939228246284, 5.1993375821928169315, 5.9978070150076868715};

  int numValues = sizeof(log_p) / sizeof(double);

  for (int i = 0; i < numValues; ++i) {
    EXPECT_NEAR(exact[i], inv_Phi_log(log_p[i]), 9.5e-14);
  }
}

TEST(MathFunctionsLog, inv_Phi_log_inf) {
  using stan::math::inv_Phi_log;
  long double log_p = std::numeric_limits<long double>::min();
  const double inf = std::numeric_limits<double>::infinity();
  EXPECT_EQ(inv_Phi_log(-inf), -inf);
  log_p = log(1.);
  EXPECT_EQ(inv_Phi_log(log_p), inf);
}

TEST(MathFunctionsLog, inv_Phi_log_nan) {
  using stan::math::inv_Phi_log;
  double nan = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(inv_Phi_log(nan), std::domain_error);
  EXPECT_THROW(inv_Phi_log(2.1), std::domain_error);
}
