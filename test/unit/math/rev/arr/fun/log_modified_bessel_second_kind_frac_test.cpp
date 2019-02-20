#include <stan/math/rev/arr/fun/log_modified_bessel_second_kind_frac.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/scal/fun/nan_util.hpp>
#include <test/unit/math/rev/scal/util.hpp>

struct TestValue {
  double v;
  double z;
  double value;
  double grad_v;
  double grad_z;

  TestValue(double _v, double _z, double _value, double _grad_v, double _grad_z)
      : v(_v), z(_z), value(_value), grad_v(_grad_v), grad_z(_grad_z) {}
};

// Test values generated by Mathematica with the following code:
//
// Needs["NumericalCalculus`"]
// vs = {1.11, 2.3,5.68, -5.42,14.1, 31.4};
// zs=  {0.08, 1.35, 5.42, 15.61 };
// Print["std::array<TestValue, ", Length[vs]*Length[zs], "> testValues = {"];
// Block[{$MaxPrecision = 80, $MinPrecision = 40}, {
//   For[i = 1, i <= Length[vs], i++, {
//     For[j = 1, j <= Length[zs], j++, {
//       v = vs[[i]];
//       z = zs[[j]];
//       val = N[Log[BesselK[v,z]]];
//       derivV1 = ND[Log[BesselK[x,z]], x,v, WorkingPrecision -> 40,
//          Terms -> 10];
//       derivV2 =ND[Log[BesselK[x,z]], x, v, Mehod -> NIntegrate,
//          WorkingPrecision -> 40,Terms -> 10];
//       If[derivV1 != derivV2,Throw["Derivatives are not equal"], True];
//
//       derivZAnalytic = N[-(BesselK[v-1,z] / BesselK[v,z]) -v/z];
//
//       Print["  TestValue(",v,",",z,",",NumberForm[val,20],",",
//          NumberForm[derivV1, 20],",",NumberForm[derivZAnalytic,20],"),"]
//     }]
//   }]
// }];
// Print["};"]
std::array<TestValue, 24> testValues
    = {TestValue(1.11, 0.08, 2.818523457573867, 2.829911601437736,
                 -14.03658220085253),

       TestValue(1.11, 1.35, -0.998097434891474, 0.6152521730077149,
                 -1.530923296686704),

       TestValue(1.11, 5.42, -5.955948396465054, 0.1879584909843152,
                 -1.106440428824266),

       TestValue(1.11, 15.61, -16.72766451823992, 0.06892414600018968,
                 -1.03392628816099),

       TestValue(2.3, 0.08, 6.863227971882178, 3.819854901094816,
                 -28.7806690845209),

       TestValue(2.3, 1.35, 0.06987248240261751, 1.154660739860992,
                 -2.093596088937962),

       TestValue(2.3, 5.42, -5.61499495879317, 0.3836097936394976,
                 -1.163530511255285),

       TestValue(2.3, 15.61, -16.6018291360965, 0.1424766644498164,
                 -1.041724687228245),

       TestValue(5.68, 0.08, 21.84072789337136, 4.86529683432042,
                 -71.00854621471304),

       TestValue(5.68, 1.35, 5.694211980821007, 2.059362720336962,
                 -4.348069083335933),

       TestValue(5.68, 5.42, -3.461732931197324, 0.871361663968448,
                 -1.494542640840436),

       TestValue(5.68, 15.61, -15.7728086430974, 0.3465321506980504,
                 -1.092248498203896),

       TestValue(-5.42, 0.08, 20.58236563679835, -4.813975363003773,
                 -67.75904881616597),

       TestValue(-5.42, 1.35, 5.165103157433492, -2.010328485614283,
                 -4.163267880991594),

       TestValue(-5.42, 5.42, -3.68396561124711, -0.838004123015669,
                 -1.462346691964963),

       TestValue(-5.42, 15.61, -15.8609147084551, -0.331194933440333,
                 -1.086946750262522),

       TestValue(14.1, 0.08, 67.50570476690205, 5.829180002337242,
                 -176.2530534042938),

       TestValue(14.1, 1.35, 27.62684751052749, 3.005980765250504,
                 -10.49582397274031),

       TestValue(14.1, 5.42, 7.514275718863735, 1.653392699685558,
                 -2.799615923544776),

       TestValue(14.1, 15.61, -10.894914609875, 0.7957534956863927,
                 -1.36543716043966),

       TestValue(31.4, 0.08, 176.4074573511386, 6.649677370800394,
                 -392.5013157871182),

       TestValue(31.4, 1.35, 87.6613625015311, 3.824335037153727,
                 -23.28145189933024),

       TestValue(31.4, 5.42, 43.79002118380178, 2.441700751731519,
                 -5.881782356655926),

       TestValue(31.4, 15.61, 8.8739384452097, 1.436044516907569,
                 -2.252872461451006)};

double allowed_rel_error = 1e-6;

::testing::AssertionResult check_relative_error(double expected,
                                                double actual) {
  ::testing::AssertionResult base_result = ::testing::AssertionSuccess();
  double difference = fabs(expected - actual);
  double threshold = fabs(allowed_rel_error * expected);
  if (difference > threshold) {
    base_result = ::testing::AssertionFailure();
  }

  return base_result << "Expected: " << expected << " actual: " << actual
                     << " difference: " << difference
                     << " threshold: " << threshold;
}

#define EXPECT_REL_ERROR(a, b) EXPECT_TRUE(check_relative_error(a, b))

TEST(AgradRev, log_modified_bessel_second_kind_frac_double_double) {
  std::ostringstream err_out;

  std::for_each(
      testValues.begin(), testValues.end(), [&err_out](TestValue test) {
        double f1
            = stan::math::log_modified_bessel_second_kind_frac(test.v, test.z);
        EXPECT_REL_ERROR(test.value, f1);
      });
}

TEST(AgradRev, log_modified_bessel_second_kind_frac_double_var) {
  std::ostringstream err_out;

  std::for_each(
      testValues.begin(), testValues.end(), [&err_out](TestValue test) {
        AVAR z(test.z);
        try {
          AVAR f = stan::math::log_modified_bessel_second_kind_frac(test.v, z);
          EXPECT_REL_ERROR(test.value, f.val());

          AVEC x = createAVEC(z);
          VEC g;
          f.grad(x, g);
          EXPECT_REL_ERROR(test.grad_z, g[0]);
        } catch (const std::domain_error& err) {
          std::cout << "\nAt v = " << test.v << ", z = " << test.z << ":\n";
          throw;
        }
      });
}

TEST(AgradRev, log_modified_bessel_second_kind_frac_var_double) {
  std::ostringstream err_out;

  std::for_each(
      testValues.begin(), testValues.end(), [&err_out](TestValue test) {
        try {
          AVAR v(test.v);
          AVAR f = stan::math::log_modified_bessel_second_kind_frac(v, test.z);
          EXPECT_REL_ERROR(test.value, f.val());

          AVEC x = createAVEC(v);
          VEC g;
          f.grad(x, g);
          EXPECT_REL_ERROR(test.grad_v, g[0]);
        } catch (const std::domain_error& err) {
          std::cout << "\nAt v = " << test.v << ", z = " << test.z << ":\n";
          throw;
        }
      });
}

TEST(AgradRev, log_modified_bessel_second_kind_frac_var_var) {
  std::ostringstream err_out;

  std::for_each(
      testValues.begin(), testValues.end(), [&err_out](TestValue test) {
        try {
          AVAR v(test.v);
          AVAR z(test.z);
          AVAR f = stan::math::log_modified_bessel_second_kind_frac(v, z);
          EXPECT_REL_ERROR(test.value, f.val());

          AVEC x = createAVEC(v, z);
          VEC g;
          f.grad(x, g);
          EXPECT_REL_ERROR(test.grad_v, g[0]);
          EXPECT_REL_ERROR(test.grad_z, g[1]);
        } catch (const std::domain_error& err) {
          std::cout << "\nAt v = " << test.v << ", z = " << test.z << ":\n";
          throw;
        }
      });
}
