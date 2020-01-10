#include <stan/math/prim/scal.hpp>
#include <stan/math/prim/scal/fun/lgamma_stirling_diff.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <test/unit/math/expect_near_rel.hpp>

TEST(MathFunctions, lgamma_stirling_diff_errors) {
  // TODO[martinmodrak] nan, negative values, ...
}

TEST(MathFunctions, lgamma_stirling_diff_accuracy) {
  using stan::math::internal::lgamma_stirling_diff_big;
  using stan::math::lgamma_stirling_diff;
  using stan::test::expect_near_rel;

  long double start = std::nextafter(10, 11);
  for (long double x = start; x < 2 * lgamma_stirling_diff_big; x *= 1.5) {
    long double stirling
        = 0.5 * std::log(2 * static_cast<long double>(stan::math::pi()))
          + (x - 0.5) * std::log(x) - x;
    long double lgamma_res = std::lgamma(x);
    double diff_actual = static_cast<double>(lgamma_res - stirling);
    double diff = lgamma_stirling_diff(x);

    std::ostringstream msg;
    msg << "x = " << x << "; lgamma = " << lgamma_res
        << "; stirling = " << stirling;
    expect_near_rel(msg.str(), diff, diff_actual, 1e-4);
  }

  double before_big = std::nextafter(lgamma_stirling_diff_big, 0);
  double after_big = std::nextafter(lgamma_stirling_diff_big,
                                    stan::math::positive_infinity());
  expect_near_rel("big cutoff", lgamma_stirling_diff(before_big),
                  lgamma_stirling_diff(after_big));
}

namespace lgamma_stirling_diff_test_internal {
struct TestValue {
  double x;
  double val;
};

std::vector<TestValue> testValues = {
    {10.049787068367863943, 0.0082893206359322200849},
    {10.135335283236612692, 0.0082193992370778811855},
    {10.367879441171442322, 0.0080351590240775284273},
    {11., 0.007573675487951840795},
    {12.718281828459045235, 0.0065508998676477812047},
    {17.389056098930650227, 0.0047917583976502091168},
    {30.085536923187667741, 0.0027697782366142911547},
    {64.598150033144239078, 0.0012900163188678122747},
    {158.41315910257660342, 0.00052604987562270910548},
    {413.42879349273512261, 0.0002015663117649479817},
    {1106.6331584284585993, 0.000075303482848309432847},
    {2990.9579870417282747, 0.000027861753118520084151},
    {8113.0839275753840077, 0.000010271474329002342817},
    {22036.465794806716517, 3.7816106313768371052e-6},
    {59884.141715197818455, 1.3915759823173656792e-6},
    {162764.79141900392081, 5.1198623858832122975e-7},
    {442423.39200892050333, 1.8835652643709875774e-7},
    {1.2026142841647767777e6, 6.9293483730078043652e-8},
    {3.2690273724721106393e6, 2.5491782061865871668e-8},
    {8.8861205205078726368e6, 9.3779206731455076929e-9},
    {2.4154962753575298215e7, 3.4499466707312036902e-9},
    {6.5659979137330511139e7, 1.2691647854325126115e-9},
    {1.7848231096318726084e8, 4.6689967696866715977e-10},
    {4.8516520540979027797e8, 1.717627983295846792e-10},
    {1.3188157444832146972e9, 6.318800308680570434e-11},
    {3.5849128561315915617e9, 2.3245567375731604902e-11},
    {9.7448034562489026e9, 8.5515663509760607976e-12},
    {2.6489122139843472294e10, 3.1459454523782780579e-12},
    {7.2004899347385872524e10, 1.1573286552529392289e-12},
};
}  // namespace lgamma_stirling_diff_test_internal

TEST(MathFunctions, lgamma_stirling_diff_precomputed) {
  using stan::math::lgamma_stirling_diff;
  using stan::test::expect_near_rel;
  using namespace lgamma_stirling_diff_test_internal;

  for (TestValue t : testValues) {
    std::ostringstream msg;
    msg << "x = " << t.x;
    expect_near_rel(msg.str(), lgamma_stirling_diff(t.x), t.val);
  }
}