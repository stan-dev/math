#include <iostream> //TODO remove
#include <stan/math/rev/arr/fun/log_modified_bessel_second_kind_frac.hpp>
#include <stan/math/rev/core/set_zero_all_adjoints.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/scal/fun/nan_util.hpp>
#include <test/unit/math/rev/scal/util.hpp>

//TODO add values close to other boundaries, v = 0
std::array<double, 21> v_to_test = { 3.15e-7, 2.62e-6, 9.2e-5, 0.0026, 0.0843, 0.17345,
  1, 1.63, 7.42, 42.42424, 148.7565, 513.6, 712.456, 714.456, 15330.75,
  stan::math::log_modified_bessel_second_kind_frac_large_v_bound1 - 1,
  stan::math::log_modified_bessel_second_kind_frac_large_v_bound1,
  stan::math::log_modified_bessel_second_kind_frac_large_v_bound1 + 1,
  stan::math::log_modified_bessel_second_kind_frac_large_v_bound2 - 1,
  stan::math::log_modified_bessel_second_kind_frac_large_v_bound2,
  stan::math::log_modified_bessel_second_kind_frac_large_v_bound2 + 1
  };

//TODO: not working for very small z
// std::array<double, 16> z_to_test = {1.48e-7, 3.6e-6, 7.248e-5, 4.32e-4,
//   8.7e-3, 0.04523, 0.17532, 1, 3, 11.32465, 105.6, 1038.4, 105321.6,
//   stan::math::log_modified_bessel_second_kind_frac_medium_z_bound - 1,
//   stan::math::log_modified_bessel_second_kind_frac_medium_z_bound,
//   stan::math::log_modified_bessel_second_kind_frac_medium_z_bound + 1}; 
std::array<double, 12> z_to_test = {
  8.7e-3, 0.04523, 0.17532, 1, 3, 11.32465, 105.6, 1038.4, 15321.6,
  stan::math::log_modified_bessel_second_kind_frac_medium_z_bound - 1,
  stan::math::log_modified_bessel_second_kind_frac_medium_z_bound,
  stan::math::log_modified_bessel_second_kind_frac_medium_z_bound + 1
  }; 

double allowed_recurrence_error = 1e-7;

// TEST(AgradRev, log_modified_bessel_second_kind_frac_can_compute_double) {  
//   std::for_each(v_to_test.begin(), v_to_test.end(), [](double v) {
//     std::for_each(z_to_test.begin(), z_to_test.end(), [v](double z) {
//       try {
//           stan::math::log_modified_bessel_second_kind_frac(v, z);
//       } catch (...) {
//         std::cout << "\nAt v = " << v << ", z = " << z << ":\n";
//         throw;
//       }
//     });
//   });
// }

// TEST(AgradRev, log_modified_bessel_second_kind_frac_can_compute_var) {  
//   std::for_each(v_to_test.begin(), v_to_test.end(), [](double v) {
//     std::for_each(z_to_test.begin(), z_to_test.end(), [v](double z) {
//       try {
//         AVAR v_var(v);
//         AVAR z_var(z);
//         stan::math::log_modified_bessel_second_kind_frac(v_var, z_var);
//       } catch (...) {
//         std::cout << "\nAt v = " << v << ", z = " << z << ":\n";
//         throw;
//       }
//     });
//   });
// }



// Using the recurrence relation (adapted to log)
// http://functions.wolfram.com/Bessel-TypeFunctions/BesselK/17/01/01/

TEST(AgradRev, log_modified_bessel_second_kind_frac_recurrence) {  
  std::for_each(v_to_test.begin(), v_to_test.end(), [](double v) {
    std::for_each(z_to_test.begin(), z_to_test.end(), [v](double z) {
      AVAR v_var(v);
      AVAR z_var(z);

      try {
        AVAR left_hand = 
          stan::math::log_modified_bessel_second_kind_frac(v_var, z_var);
        AVAR right_hand;

        AVAR log_K_vm1 = stan::math::log_modified_bessel_second_kind_frac(v_var - 1, z_var);
        AVAR log_K_vm2 = stan::math::log_modified_bessel_second_kind_frac(v_var - 2, z_var);
        if(fabs(v - 1) < 1e-4) {
          AVAR log_K_vp1 = stan::math::log_modified_bessel_second_kind_frac(v_var + 1, z_var);
          AVAR log_K_vp2 = stan::math::log_modified_bessel_second_kind_frac(v_var + 2, z_var);
          right_hand = stan::math::log_diff_exp(
            log_K_vp2,  
            stan::math::LOG_2 + stan::math::log(v_var + 1) - stan::math::log(z_var)
            + log_K_vp1);
        } else if(v > 1) {
          right_hand = 
            stan::math::log_sum_exp(
              log_K_vm2,
              stan::math::LOG_2 + stan::math::log(v_var - 1) - stan::math::log(z_var)
              + log_K_vm1
          );
        } else {
          right_hand = 
            stan::math::log_diff_exp(
              log_K_vm2,
              stan::math::LOG_2 + stan::math::log(stan::math::fabs(v_var - 1)) - stan::math::log(z_var)
              + log_K_vm1
            );
        } 

        AVAR ratio = left_hand / right_hand;

      std::cout << "v = " << v << "; z = " << z << "; K_v = " << left_hand << "; K_vm1 = " << log_K_vm1 << "; K_vm2 = " << log_K_vm2 <<  std::endl;

        EXPECT_NEAR(ratio.val(), 1.0, allowed_recurrence_error);

        AVEC x = createAVEC(v_var, z_var, left_hand, right_hand);
        VEC g;
        ratio.grad(x, g);
        EXPECT_NEAR(1 + g[0], 1.0, allowed_recurrence_error);
        EXPECT_NEAR(1 + g[1], 1.0, allowed_recurrence_error);
      } catch (...) {
        std::cout << "\nAt v = " << v << ", z = " << z << ":\n";
        throw;
      }
    });
  });
}

struct fun {
  template <typename T_v, typename T_z>
  inline typename boost::math::tools::promote_args<T_v, T_z>::type 
  operator()(const T_v& arg1, const T_z& arg2) const {
    return stan::math::log_modified_bessel_second_kind_frac(arg1, arg2);
  }
};


TEST(AgradRev, log_modified_bessel_second_kind_frac_input) {  
  fun f;
  test_nan(f, 1, 1, true, false);
  EXPECT_THROW(stan::math::log_modified_bessel_second_kind_frac(1.0,-1.0), std::domain_error);
  EXPECT_THROW(stan::math::log_modified_bessel_second_kind_frac(std::numeric_limits<double>::infinity(),1), std::domain_error);
  EXPECT_THROW(stan::math::log_modified_bessel_second_kind_frac(1.0,std::numeric_limits<double>::infinity()), std::domain_error);

  std::for_each(v_to_test.begin(), v_to_test.end(), [](double v) {
    EXPECT_TRUE(std::isinf(stan::math::log_modified_bessel_second_kind_frac(v,0))) << "Infinity for z = 0";
  });
}

// struct TestValue {
//   double v;
//   double z;
//   double value;
//   double grad_v;
//   double grad_z;

//   TestValue(double _v, double _z, double _value, double _grad_v, double _grad_z)
//       : v(_v), z(_z), value(_value), grad_v(_grad_v), grad_z(_grad_z) {}
// };

// // Test values generated by Mathematica with the following code:
// //
// // Needs["NumericalCalculus`"]
// // vs = {1.11, 2.3,5.68, -5.42,14.1, 31.4};
// // zs=  {0.08, 1.35, 5.42, 15.61 };
// // Print["std::array<TestValue, ", Length[vs]*Length[zs], "> testValues = {"];
// // Block[{$MaxPrecision = 80, $MinPrecision = 40}, {
// //   For[i = 1, i <= Length[vs], i++, {
// //     For[j = 1, j <= Length[zs], j++, {
// //       v = vs[[i]];
// //       z = zs[[j]];
// //       val = N[Log[BesselK[v,z]]];
// //       derivV1 = ND[Log[BesselK[x,z]], x,v, WorkingPrecision -> 40,
// //          Terms -> 10];
// //       derivV2 =ND[Log[BesselK[x,z]], x, v, Mehod -> NIntegrate,
// //          WorkingPrecision -> 40,Terms -> 10];
// //       If[derivV1 != derivV2,Throw["Derivatives are not equal"], True];
// //
// //       derivZAnalytic = N[-(BesselK[v-1,z] / BesselK[v,z]) -v/z];
// //
// //       Print["  TestValue(",v,",",z,",",NumberForm[val,20],",",
// //          NumberForm[derivV1, 20],",",NumberForm[derivZAnalytic,20],"),"]
// //     }]
// //   }]
// // }];
// // Print["};"]
// std::array<TestValue, 25> testValues
//     = {TestValue(1.11, 0.08, 2.818523457573867, 2.829911601437736,
//                  -14.03658220085253),

//        TestValue(-6.1,1,8.44532282091819,-2.426667154237633,
//                  -6.196902250273741),

//        TestValue(1.11, 1.35, -0.998097434891474, 0.6152521730077149,
//                  -1.530923296686704),

//        TestValue(1.11, 5.42, -5.955948396465054, 0.1879584909843152,
//                  -1.106440428824266),

//        TestValue(1.11, 15.61, -16.72766451823992, 0.06892414600018968,
//                  -1.03392628816099),

//        TestValue(2.3, 0.08, 6.863227971882178, 3.819854901094816,
//                  -28.7806690845209),

//        TestValue(2.3, 1.35, 0.06987248240261751, 1.154660739860992,
//                  -2.093596088937962),

//        TestValue(2.3, 5.42, -5.61499495879317, 0.3836097936394976,
//                  -1.163530511255285),

//        TestValue(2.3, 15.61, -16.6018291360965, 0.1424766644498164,
//                  -1.041724687228245),

//        TestValue(5.68, 0.08, 21.84072789337136, 4.86529683432042,
//                  -71.00854621471304),

//        TestValue(5.68, 1.35, 5.694211980821007, 2.059362720336962,
//                  -4.348069083335933),

//        TestValue(5.68, 5.42, -3.461732931197324, 0.871361663968448,
//                  -1.494542640840436),

//        TestValue(5.68, 15.61, -15.7728086430974, 0.3465321506980504,
//                  -1.092248498203896),

//        TestValue(-5.42, 0.08, 20.58236563679835, -4.813975363003773,
//                  -67.75904881616597),

//        TestValue(-5.42, 1.35, 5.165103157433492, -2.010328485614283,
//                  -4.163267880991594),

//        TestValue(-5.42, 5.42, -3.68396561124711, -0.838004123015669,
//                  -1.462346691964963),

//        TestValue(-5.42, 15.61, -15.8609147084551, -0.331194933440333,
//                  -1.086946750262522),

//        TestValue(14.1, 0.08, 67.50570476690205, 5.829180002337242,
//                  -176.2530534042938),

//        TestValue(14.1, 1.35, 27.62684751052749, 3.005980765250504,
//                  -10.49582397274031),

//        TestValue(14.1, 5.42, 7.514275718863735, 1.653392699685558,
//                  -2.799615923544776),

//        TestValue(14.1, 15.61, -10.894914609875, 0.7957534956863927,
//                  -1.36543716043966),

//        TestValue(31.4, 0.08, 176.4074573511386, 6.649677370800394,
//                  -392.5013157871182),

//        TestValue(31.4, 1.35, 87.6613625015311, 3.824335037153727,
//                  -23.28145189933024),

//        TestValue(31.4, 5.42, 43.79002118380178, 2.441700751731519,
//                  -5.881782356655926),

//        TestValue(31.4, 15.61, 8.8739384452097, 1.436044516907569,
//                  -2.252872461451006)};

// double allowed_rel_error = 1e-6;

// ::testing::AssertionResult check_relative_error(double expected,
//                                                 double actual) {
//   ::testing::AssertionResult base_result = ::testing::AssertionSuccess();
//   double difference = fabs(expected - actual);
//   double threshold = fabs(allowed_rel_error * expected);
//   if (difference > threshold) {
//     base_result = ::testing::AssertionFailure();
//   }

//   return base_result << "Expected: " << expected << " actual: " << actual
//                      << " difference: " << difference
//                      << " threshold: " << threshold;
// }

// #define EXPECT_REL_ERROR(a, b) EXPECT_TRUE(check_relative_error(a, b))

// TEST(AgradRev, log_modified_bessel_second_kind_frac_double_double) {
//   std::for_each(
//       testValues.begin(), testValues.end(), [](TestValue test) {
//         double f1
//             = stan::math::log_modified_bessel_second_kind_frac(test.v, test.z);
//         EXPECT_REL_ERROR(test.value, f1);
//       });
// }

// TEST(AgradRev, log_modified_bessel_second_kind_frac_double_var) {
//   std::for_each(
//       testValues.begin(), testValues.end(), [](TestValue test) {
//         AVAR z(test.z);
//         try {
//           AVAR f = stan::math::log_modified_bessel_second_kind_frac(test.v, z);
//           EXPECT_REL_ERROR(test.value, f.val());

//           AVEC x = createAVEC(z);
//           VEC g;
//           f.grad(x, g);
//           EXPECT_REL_ERROR(test.grad_z, g[0]);
//         } catch (const std::domain_error& err) {
//           std::cout << "\nAt v = " << test.v << ", z = " << test.z << ":\n";
//           throw;
//         }
//       });
// }

// TEST(AgradRev, log_modified_bessel_second_kind_frac_var_double) {
//   std::for_each(
//       testValues.begin(), testValues.end(), [](TestValue test) {
//         try {
//           AVAR v(test.v);
//           AVAR f = stan::math::log_modified_bessel_second_kind_frac(v, test.z);
//           EXPECT_REL_ERROR(test.value, f.val());

//           AVEC x = createAVEC(v);
//           VEC g;
//           f.grad(x, g);
//           EXPECT_REL_ERROR(test.grad_v, g[0]);
//         } catch (const std::domain_error& err) {
//           std::cout << "\nAt v = " << test.v << ", z = " << test.z << ":\n";
//           throw;
//         }
//       });
// }

// TEST(AgradRev, log_modified_bessel_second_kind_frac_var_var) {
//   std::for_each(
//       testValues.begin(), testValues.end(), [](TestValue test) {
//         try {
//           AVAR v(test.v);
//           AVAR z(test.z);
//           AVAR f = stan::math::log_modified_bessel_second_kind_frac(v, z);
//           EXPECT_REL_ERROR(test.value, f.val());

//           AVEC x = createAVEC(v, z);
//           VEC g;
//           f.grad(x, g);
//           EXPECT_REL_ERROR(test.grad_v, g[0]);
//           EXPECT_REL_ERROR(test.grad_z, g[1]);
//         } catch (const std::domain_error& err) {
//           std::cout << "\nAt v = " << test.v << ", z = " << test.z << ":\n";
//           throw;
//         }
//       });
// }
