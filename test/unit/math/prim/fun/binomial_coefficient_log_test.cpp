#include <test/unit/math/expect_near_rel.hpp>
#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

template <typename T_N, typename T_n>
void test_binom_coefficient(const T_N& N, const T_n& n) {
  using stan::math::binomial_coefficient_log;
  EXPECT_FLOAT_EQ(lgamma(N + 1) - lgamma(n + 1) - lgamma(N - n + 1),
                  binomial_coefficient_log(N, n))
      << "N = " << N << ", n = " << n;
}

TEST(MathFunctions, binomial_coefficient_log) {
  using stan::math::binomial_coefficient_log;
  EXPECT_FLOAT_EQ(1.0, exp(binomial_coefficient_log(2.0, 2.0)));
  EXPECT_FLOAT_EQ(2.0, exp(binomial_coefficient_log(2.0, 1.0)));
  EXPECT_FLOAT_EQ(3.0, exp(binomial_coefficient_log(3.0, 1.0)));
  EXPECT_NEAR(3.0, exp(binomial_coefficient_log(3.0, 2.0)), 0.0001);

  EXPECT_FLOAT_EQ(29979.16, binomial_coefficient_log(100000, 91116));

  EXPECT_EQ(binomial_coefficient_log(-1, 0), 0); // Needed for neg_binomial_2
  EXPECT_EQ(binomial_coefficient_log(50, 0), 0);
  EXPECT_EQ(binomial_coefficient_log(10000, 0), 0);

  EXPECT_EQ(binomial_coefficient_log(10, 11), -std::numeric_limits<double>::infinity()); 
  EXPECT_EQ(binomial_coefficient_log(10, -1), -std::numeric_limits<double>::infinity()); 

  for (int n = 0; n < 1010; ++n) {
    test_binom_coefficient(1010, n);
    test_binom_coefficient(1010.0, n);
    test_binom_coefficient(1010, static_cast<double>(n));
    test_binom_coefficient(1010.0, static_cast<double>(n));
  }

  test_binom_coefficient(1e9, 1e5);
  test_binom_coefficient(1e50, 1e45);
  test_binom_coefficient(1e20, 1e15);
}

TEST(MathFunctions, binomial_coefficient_log_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_TRUE(std::isnan(stan::math::binomial_coefficient_log(2.0, nan)));
  EXPECT_TRUE(std::isnan(stan::math::binomial_coefficient_log(nan, 2.0)));
  EXPECT_TRUE(std::isnan(stan::math::binomial_coefficient_log(nan, nan)));
}

TEST(MathFunctions, binomial_coefficient_log_errors) {
  using stan::math::binomial_coefficient_log;
  EXPECT_NO_THROW(binomial_coefficient_log(10, 11));
  EXPECT_THROW(binomial_coefficient_log(10, 11.01), std::domain_error);
  EXPECT_THROW(binomial_coefficient_log(-1, 0.3), std::domain_error);
  EXPECT_THROW(binomial_coefficient_log(10, -1.1), std::domain_error);
  EXPECT_NO_THROW(binomial_coefficient_log(10, -0.9));
}

namespace binomial_coefficient_test_internal {
struct TestValue {
  double n;
  double k;
  double val;
};

std::vector<TestValue> testValues = {  {0.00003,-0.9,-2.21375637737528044964112},  {0.00003,3.e-15,1.48040820535563428588846e-19},  {0.00003,3.e-10,1.48039340142161568666212e-14},  {0.00003,3.e-7,1.46560412344434242311062e-11},  {0.00003,6.e-6,2.36865312869367134774661e-10},  {0.00003,0.000015,3.70102051348524044535716e-10},  {0.00003,0.0000291,4.30798787797857754693695e-11},  {0.00003,0.0000299999982,8.88244870007509649977929e-17},  {0.00003,0.00002999999991,4.44122460318735146597007e-18},  {0.00003,0.90003,-2.21375637737528044964112},  {0.002,-0.9,-2.21559326412971099686943},  {0.002,2.e-13,6.57013709556564677684856e-16},  {0.002,2.e-8,6.57007139476544583161173e-11},  {0.002,0.00002,6.50443564072189575550994e-8},  {0.002,0.0004,1.05122171458287350859763e-6},  {0.002,0.001,1.64253373496215313253469e-6},  {0.002,0.00194,1.91190982195918976356429e-7},  {0.002,0.00199999988,3.94208202120835082737684e-13},  {0.002,0.001999999994,1.97104112295366725515452e-14},  {0.002,0.902,-2.21559326412971099686943},  {1,-0.9,-2.85558226198351740582195},  {1,1.e-10,9.9999999988550659331851e-11},  {1,0.00001,9.99988550692664559909352e-6},  {1,0.01,0.00988583703486131052627978},  {1,0.2,0.156457962917688016707705},  {1,0.5,0.241564475270490444691037},  {1,0.97,0.028978328236256312960776},  {1,0.99999994,5.99999958782374313463811e-8},  {1,0.999999997,2.99999998969559340736596e-9},  {1,1.9,-2.85558226198351740582195},  {8,-0.9,-4.22528965320883461943031},  {8,8.e-10,2.17428571372173153982474e-9},  {8,0.00008,0.000217422931805073420417006},  {8,0.08,0.211982267378255838975509},  {8,1.6,2.90678606291134283426263},  {8,4,4.24849524204935898912334},  {8,7.76,0.606274586245453651115361},  {8,7.99999952,1.30457122553768403613331e-6},  {8,7.999999976,6.52285709209869625945566e-8},  {8,8.9,-4.22528965320883461943031},  {1325,-0.9,-8.72360867216657209762532},  {1325,1.325e-7,1.02909578020477960435539e-6},  {1325,0.01325,0.102766042691370430370992},  {1325,13.25,71.9898321274090629975055},  {1325,265,659.435649329029419323398},  {1325,662.5,914.599450340845275100724},  {1325,1285.25,175.786260651191862665015},  {1325,1324.9999205,0.000617452276410452190170437},  {1325,1324.999996025,0.0000308728608380968862741097},  {1325,1325.9,-8.72360867216657209762532},  {845000,-0.9,-14.5350963792733464918229},  {845000,0.0000845,0.00120194816738712136581358},  {845000,8.45,103.738303827736743600251},  {845000,8450,47315.8616457576200611209},  {845000,169000,422833.221695496506553128},  {845000,422500,585702.318235552114086514},  {845000,819650,113851.158132678562120685},  {845000,844999.9493,0.719108776819481762797449},  {845000,844999.997465,0.036053342347290003917417},  {845000,845000.9,-14.5350963792733464918229},  {3000000000000000.0,3,105.120406581508328854183},  {3000000000000000.0,100,3199.99949280231435502243},  {3000000000000000.0,12895,350387.5243605883687667},  {100000000000000000000.0,3,136.363346110414686040237},  {100000000000000000000.0,100,4241.4308104325278778424},  {100000000000000000000.0,12895,484680.092769031900296878},};
}  // namespace binomial_coefficient_test_internal

TEST(MathFunctions, binomial_coefficient_log_precomputed) {
  using binomial_coefficient_test_internal::TestValue;
  using binomial_coefficient_test_internal::testValues;
  using stan::test::expect_near_rel;

  for (TestValue t : testValues) {

    std::ostringstream msg;
    msg << std::setprecision(22) << "n = " << t.n << ", k = " << t.k;

    double val = stan::math::binomial_coefficient_log(t.n, t.k);
    expect_near_rel(msg.str(), val, t.val);

    if(t.k > t.n / 10) {
      // Testing the mirrored binomial coefficient is not performed for
      // small k, as we may loose so much precision computing k2
      // that the test becomes invalid
      std::ostringstream msg2;
      double k2 = t.n - t.k;
      msg2 << std::setprecision(22) << "n = " << t.n << ", k = " << k2;
      double val2 = stan::math::binomial_coefficient_log(t.n, k2);
      expect_near_rel(msg2.str(), val2, t.val);
    }
  }
}
