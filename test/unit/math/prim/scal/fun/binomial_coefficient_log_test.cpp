#include <test/unit/math/expect_near_rel.hpp>
#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

template <typename T_N, typename T_n>
void test_binom_coefficient(const T_N& N, const T_n& n) {
  using stan::math::binomial_coefficient_log;
  EXPECT_FLOAT_EQ(lgamma(N + 1) - lgamma(n + 1) - lgamma(N - n + 1),
                  binomial_coefficient_log(N, n));
}

TEST(MathFunctions, binomial_coefficient_log) {
  using stan::math::binomial_coefficient_log;
  EXPECT_FLOAT_EQ(1.0, exp(binomial_coefficient_log(2.0, 2.0)));
  EXPECT_FLOAT_EQ(2.0, exp(binomial_coefficient_log(2.0, 1.0)));
  EXPECT_FLOAT_EQ(3.0, exp(binomial_coefficient_log(3.0, 1.0)));
  EXPECT_NEAR(3.0, exp(binomial_coefficient_log(3.0, 2.0)), 0.0001);

  EXPECT_FLOAT_EQ(29979.16, binomial_coefficient_log(100000, 91116));

  EXPECT_EQ(binomial_coefficient_log(50, 0), 0);
  EXPECT_EQ(binomial_coefficient_log(10000, 0), 0);

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

namespace binomial_coefficient_test_internal {
  struct TestValue {
    double n;
    double k;
    double val;
  };  

std::vector<TestValue> testValues = {  {0.00003,3.e-20,1.4804082055036736261e-24},  {0.00003,3.e-15,1.4804082053556342859e-19},  {0.00003,3.e-10,1.4803934014216156867e-14},  {0.00003,3.e-7,1.4656041234443424231e-11},  {0.00003,6.e-6,2.3686531286936713477e-10},  {0.00003,0.000015,3.7010205134852404454e-10},  {0.00003,0.0000291,4.3079878779785775469e-11},  {0.00003,0.0000299999982,8.8824487000750964998e-17},  {0.00003,0.00003,1.4804082055036751063e-28},  {0.002,2.e-18,6.5701370962226547789e-21},  {0.002,2.e-13,6.5701370955656467768e-16},  {0.002,2.e-8,6.5700713947654458316e-11},  {0.002,0.00002,6.5044356407218957555e-8},  {0.002,0.0004,1.0512217145828735086e-6},  {0.002,0.001,1.6425337349621531325e-6},  {0.002,0.00194,1.9119098219591897636e-7},  {0.002,0.00199999988,3.9420820212083508274e-13},  {0.002,0.002,6.5701370962226613484e-25},  {1,1.e-15,9.999999999999988551e-16},  {1,1.e-10,9.999999998855065933e-11},  {1,0.00001,9.9998855069266455991e-6},  {1,0.01,0.0098858370348613105263},  {1,0.2,0.15645796291768801671},  {1,0.5,0.24156447527049044469},  {1,0.97,0.028978328236256312961},  {1,0.99999994,5.9999995878237431346e-8},  {1,1.,9.999999999999999999e-20},  {8,8.e-15,2.1742857142857086459e-14},  {8,8.e-10,2.1742857137217315398e-9},  {8,0.00008,0.00021742293180507342042},  {8,0.08,0.21198226737825583898},  {8,1.6,2.9067860629113428343},  {8,4,4.2484952420493589891},  {8,7.76,0.60627458624545365112},  {8,7.99999952,1.3045712255376840361e-6},  {8,8.,2.1742857142857142852e-18},  {1325,1.325e-12,1.0290957946506935679e-11},  {1325,1.325e-7,1.0290957802047796044e-6},  {1325,0.01325,0.10276604269137043037},  {1325,13.25,71.989832127409062998},  {1325,265,659.43564932902941932},  {1325,662.5,914.5994503408452751},  {1325,1285.25,175.78626065119186267},  {1325,1324.9999205,0.00061745227641045219017},  {1325,1325.,1.029095794650838014e-15},  {845000,8.45e-10,1.2019540397111151912e-8},  {845000,0.0000845,0.0012019481673871213658},  {845000,8.45,103.7383038277367436},  {845000,8450,47315.861645757620061},  {845000,169000,422833.22169549650655},  {845000,422500,585702.31823555211409},  {845000,819650,113851.15813267856212},  {845000,844999.9493,0.7191087768194817628},  {845000,845000.,1.2019540397698355631e-12},  {3000000000000000.0,3,105.12040658150832885},  {3000000000000000.0,100,3199.999492802314355},  {3000000000000000.0,12895,350387.52436058836877},  {100000000000000000000.0,3,136.36334611041468604},  {100000000000000000000.0,100,4241.4308104325278778},  {100000000000000000000.0,12895,484680.0927690319003},};

}

TEST(MathFunctions, binomial_coefficient_log_precomputed) {
  using binomial_coefficient_test_internal::TestValue;
  using binomial_coefficient_test_internal::testValues;
  using stan::test::expect_near_rel;

  for (TestValue t : testValues) {
    double val = stan::math::binomial_coefficient_log(t.n, t.k);

    std::ostringstream msg;
    msg << "n = " << t.n << ", k = " << t.k;
    expect_near_rel(msg.str(), val, t.val);
  }
}
