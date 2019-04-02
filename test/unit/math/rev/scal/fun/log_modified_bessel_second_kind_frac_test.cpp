#include <stan/math/rev/scal/fun/log_modified_bessel_second_kind_frac.hpp>
#include <stan/math/rev/scal/fun/modified_bessel_second_kind.hpp>
#include <stan/math/rev/core/set_zero_all_adjoints.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/scal/fun/nan_util.hpp>
#include <test/unit/math/rev/scal/util.hpp>
#include <iostream>
#include <fstream>
#include <limits>

// Set to true to write CSV file with recurrence test results for analysis
bool output_debug_csv = true;

using stan::math::LOG_2;
using stan::math::besselk_internal::ComputationType;
using stan::math::besselk_internal::choose_computation_type;
using stan::math::besselk_internal::rothwell_max_v;
using stan::math::besselk_internal::small_z_min_v;
using stan::math::log_diff_exp;
using stan::math::log_modified_bessel_second_kind_frac;
using stan::math::log_sum_exp;
using stan::math::recover_memory;
using stan::math::var;

// TODO(martinmodrak) add values close to decision boundaries
std::array<double, 33> v_to_test = {0,
                                    3.15e-7,
                                    2.62e-6,
                                    9.2e-5,
                                    0.0026,
                                    0.0843,
                                    0.17345,
                                    1,
                                    1.63,
                                    7.42,
                                    42.42424,
                                    86.5,
                                    113.8,
                                    148.7565,
                                    180.6,
                                    246.3,
                                    300.5,
                                    513.6,
                                    712.456,
                                    714.456,
                                    1235.6,
                                    8656,
                                    15330.75,
                                    37634.2,
                                    85323,
                                    rothwell_max_v - 1,
                                    rothwell_max_v - 1e-8,
                                    rothwell_max_v,
                                    rothwell_max_v + 1,
                                    small_z_min_v - 1,
                                    small_z_min_v - 1e-8,
                                    small_z_min_v,
                                    small_z_min_v + 1};

std::array<double, 19> z_to_test
    = {1.48e-7, 3.6e-6,   7.248e-5, 4.32e-4, 8.7e-3, 0.04523, 0.17532,
       1,       3,        11.32465, 105.6,   1038.4, 4236,    11457.6,
       62384,   105321.6, 158742.3, 196754,  1.98e6};

double allowed_recurrence_error = 1e-7;

const char* computation_type_to_string(ComputationType c) {
  switch (c) {
    case ComputationType::Rothwell:
      return "Rothwell";
    case ComputationType::Asymp_v:
      return "Asymp_v";
    case ComputationType::Asymp_z:
      return "Asymp_z";
    case ComputationType::Asymp_small_z_relative:
      return "Small_z";
    default:
      return "Unknown";
  }
}

// TEST(AgradRev, log_modified_bessel_second_kind_frac_can_compute_double) {
//   std::for_each(v_to_test.begin(), v_to_test.end(), [](double v) {
//     std::for_each(z_to_test.begin(), z_to_test.end(), [v](double z) {
//       try {
//           log_modified_bessel_second_kind_frac(v, z);
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
//         log_modified_bessel_second_kind_frac(v_var, z_var);
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
  std::ostream* debug_output = 0;
  if (output_debug_csv) {
    debug_output = new std::ofstream("log_besselk_test.csv");
    *debug_output << "v,z,method,ratio,grad_z,grad_v,value_m2,value_m1,"
                     "value,value_p1,value_p2"
                  << std::endl;
  }

  for (auto v_iter = v_to_test.begin(); v_iter != v_to_test.end(); ++v_iter) {
    for (auto z_iter = z_to_test.begin(); z_iter != z_to_test.end(); ++z_iter) {
      for (int sign = -1; sign <= 1; sign += 2) {
        double v = sign * (*v_iter);
        double z = *z_iter;
        AVAR v_var(v);
        AVAR z_var(z);

        try {
          AVAR left_hand = log_modified_bessel_second_kind_frac(v_var, z_var);
          AVAR right_hand;

          AVAR log_K_vm1
              = log_modified_bessel_second_kind_frac(v_var - 1, z_var);
          AVAR log_K_vm2
              = log_modified_bessel_second_kind_frac(v_var - 2, z_var);

          AVAR log_K_vp1
              = log_modified_bessel_second_kind_frac(v_var + 1, z_var);
          AVAR log_K_vp2
              = log_modified_bessel_second_kind_frac(v_var + 2, z_var);

          // Trying to find the most numerically stable way to compute the
          // recursive formula
          if (v > 0) {
            if (v < 1 + 1e-4) {
              if (z > 1) {
                right_hand = log_diff_exp(
                   log_K_vp2, LOG_2 + log(v_var + 1) - log(z_var) + log_K_vp1);
              } else {
                right_hand = log(z_var) - log(2) - log(v_var) +
                  log_diff_exp(log_K_vp1, log_K_vm1);
              }
            } else {
              right_hand = log_sum_exp(
                  log_K_vm2, LOG_2 + log(v_var - 1) - log(z_var) + log_K_vm1);
            }
          } else {
            if (v > -1 - 1e-4) {
              if (z > 1) {
                right_hand = log_diff_exp(log_K_vm2, 
                  LOG_2 + log(-v_var + 1) - log(z_var) + log_K_vm1);
              } else {
                if (v == 0) {
                  AVAR right_hand_base =
                    log(modified_bessel_second_kind(0, z_var));
                  right_hand = var(new stan::math::precomp_vv_vari(
                    value_of(right_hand_base), v_var.vi_, right_hand_base.vi_,
                      0, 1));
                } else {
                  right_hand = log(z_var) - log(2) - log(-v_var) +
                    log_diff_exp(log_K_vm1, log_K_vp1);
                }
              }
            } else {
              right_hand = log_sum_exp(
                  log_K_vp2, LOG_2 + log(-v_var - 1) - log(z_var) + log_K_vp1);
            }
          }

          AVAR ratio = left_hand / right_hand;

          EXPECT_NEAR(ratio.val(), 1.0, allowed_recurrence_error);

          AVEC x = createAVEC(v_var, z_var, left_hand, right_hand);
          VEC g;
          ratio.grad(x, g);
          EXPECT_NEAR(g[0], 0, allowed_recurrence_error);
          EXPECT_NEAR(g[1], 0, allowed_recurrence_error);
          if (debug_output != 0) {
            *debug_output << v << "," << z << ","
                          << computation_type_to_string(
                                 choose_computation_type(v, z))
                          << "," << ratio.val() << "," << g[0] << "," << g[1]
                          << "," << log_K_vm2 << "," << log_K_vm1 << ","
                          << left_hand << "," << log_K_vp1 << "," << log_K_vp2
                          << std::endl;
          }
          recover_memory();
        } catch (...) {
          std::cout << "\nAt v = " << v << ", z = " << z << " method = "
                    << computation_type_to_string(choose_computation_type(v, z))
                    << ":" << std::endl;
          throw;
        }
      }
    }
  }
}

struct fun {
  template <typename T_v, typename T_z>
  inline typename boost::math::tools::promote_args<T_v, T_z>::type operator()(
      const T_v& arg1, const T_z& arg2) const {
    return log_modified_bessel_second_kind_frac(arg1, arg2);
  }
};

TEST(AgradRev, log_modified_bessel_second_kind_frac_input) {
  fun f;
  test_nan(f, 1, 1, true, false);
  EXPECT_THROW(log_modified_bessel_second_kind_frac(1.0, -1.0),
               std::domain_error);
  EXPECT_THROW(log_modified_bessel_second_kind_frac(
                   std::numeric_limits<double>::infinity(), 1),
               std::domain_error);
  EXPECT_THROW(log_modified_bessel_second_kind_frac(
                   1.0, std::numeric_limits<double>::infinity()),
               std::domain_error);

  std::for_each(v_to_test.begin(), v_to_test.end(), [](double v) {
    EXPECT_TRUE(std::isinf(log_modified_bessel_second_kind_frac(v, 0)))
        << "Infinity for z = 0";
  });
}

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
// vs = {3*10^-8, 8.42*10^-6,-6.1, 2.3,5.68, -5.42,14.1, 31.4};
// zs=  {4*10^-8,2.56*10^-6,1, 1.35, 5.42, 15.61 };
// Print["std::array<TestValue, ", Length[vs]*Length[zs], "> testValues = {"];
//   Block[{$MaxPrecision = 80, $MinPrecision = 40}, {
//     For[i = 1, i <= Length[vs], i++, {
//       For[j = 1, j <= Length[zs], j++, {
//         v = vs[[i]];
//         z = zs[[j]];
//         val = N[Log[BesselK[v,z]]];
//         derivV1 = ND[Log[BesselK[x,z]], x,v,
//           WorkingPrecision -> 40, Terms -> 10];
//         derivV2 =ND[Log[BesselK[x,z]], x, v, Mehod -> NIntegrate,
//           WorkingPrecision -> 40,Terms -> 10];
//         If[derivV1 != derivV2,Throw["Derivatives are not equal"], True];

//         derivZAnalytic = N[-(BesselK[v-1,z] / BesselK[v,z]) -v/z];

//         Print["  TestValue(",CForm[v],",",CForm[z],",",CForm[val],",",
//           CForm[derivV1],",",CForm[derivZAnalytic],"),"]
//       }]
//     }]
//   }];
// Print["};"]
std::array<TestValue, 48> testValues = {
    TestValue(3.e-8, 4.e-8, 2.842016709795352,
              2.965102177697192015626573541625171847363e-6,
              -1.4576989270969522e6),

    TestValue(3.e-8, 2.56e-6, 2.564290279944925, 1.7334431640179288e-6,
              -30067.88746129086),

    TestValue(3.e-8, 1, -0.8650643989067878,
              2.19330054363531404298968678241041735235e-8, -1.429625398260402),

    TestValue(3.e-8, 1.35, -1.345899369014985, 1.7261838025937896e-8,
              -1.3274527125104543),

    TestValue(3.e-8, 5.42, -6.060513933302737, 5.1040390297908136e-9,
              -1.0886145077261176),

    TestValue(3.e-8, 15.61, -16.76593131015518, 1.862584828000698e-9,
              -1.031547715742656),

    TestValue(8.419999999999999e-6, 4.e-8, 2.8420167133274603,
              0.0008389671519545296, -1.4576989372317046e6),

    TestValue(8.419999999999999e-6, 2.56e-6, 2.5642902819953073,
              0.00048703094653465367, -30067.887581282954),

    TestValue(8.419999999999999e-6, 1, -0.8650643988808719,
              6.155863540697342e-6, -1.4296253982806608),

    TestValue(8.419999999999999e-6, 1.35, -1.3458993689945886,
              4.8448467223866964e-6, -1.327452712522753),

    TestValue(8.419999999999999e-6, 5.42, -6.060513933296705,
              1.4325944074403853e-6, -1.0886145077271503),

    TestValue(8.419999999999999e-6, 15.61, -16.76593131015298,
              5.232112423113283e-7, -1.0315477157427932),

    TestValue(-6.1, 4.e-8, 112.40381226433932, -19.451621523830077,
              -1.5249999999999988e8),

    TestValue(-6.1, 2.56e-6, 87.034625455845, -15.29273844045875,
              -2.3828125000002505e6),

    TestValue(-6.1, 1, 8.44532282091819, -2.4266671542376326,
              -6.196902250273741),

    TestValue(-6.1, 1.35, 6.57501751331719, -2.134067641376717,
              -4.648139310106594),

    TestValue(-6.1, 5.42, -3.0846894347387646, -0.9237850302846817,
              -1.5483090878535395),

    TestValue(-6.1, 15.61, -15.622087849008592, -0.3711535955458681,
              -1.1012821947806557),

    TestValue(2.3, 4.e-8, 40.234369470202246, 18.327573443762898,
              -5.750000000000001e7),

    TestValue(2.3, 2.56e-6, 30.668938378473744, 14.168690360397694,
              -898437.5000009845),

    TestValue(2.3, 1, 0.8839980660421263, 1.3959577532638583,
              -2.6154838303420496),

    TestValue(2.3, 1.35, 0.06987248240261751, 1.1546607398609916,
              -2.0935960889379617),

    TestValue(2.3, 5.42, -5.61499495879317, 0.3836097936394976,
              -1.1635305112552852),

    TestValue(2.3, 15.61, -16.6018291360965, 0.14247666444981644,
              -1.041724687228245),

    TestValue(5.68, 4.e-8, 104.25024571265277, 19.373881532593177, -1.42e8),

    TestValue(5.68, 2.56e-6, 80.62778979916949, 15.21499844922185,
              -2.2187500000002733e6),

    TestValue(5.68, 1, 7.441890184577832, 2.3506506024792277,
              -5.785344508161015),

    TestValue(5.68, 1.35, 5.694211980821007, 2.0593627203369618,
              -4.348069083335933),

    TestValue(5.68, 5.42, -3.4617329311973237, 0.8713616639684476,
              -1.4945426408404363),

    TestValue(5.68, 15.61, -15.772808643097397, 0.34653215069805043,
              -1.0922484982038956),

    TestValue(-5.42, 4.e-8, 99.21965255139615, -19.322551217399283,
              -1.3549999999999988e8),

    TestValue(-5.42, 2.56e-6, 76.67850623958635, -15.16366813402795,
              -2.1171875000002896e6),

    TestValue(-5.42, 1, 6.837171294885847, -2.300631355301549,
              -5.531331879830391),

    TestValue(-5.42, 1.35, 5.165103157433492, -2.0103284856142833,
              -4.163267880991594),

    TestValue(-5.42, 5.42, -3.68396561124711, -0.838004123015669,
              -1.4623466919649633),

    TestValue(-5.42, 15.61, -15.8609147084551, -0.331194933440333,
              -1.086946750262522),

    TestValue(14.1, 4.e-8, 272.0779010168817, 20.33782841756488, -3.525e8),

    TestValue(14.1, 2.56e-6, 213.4376495415102, 16.178945334179055,
              -5.507812500000098e6),

    TestValue(14.1, 1, 31.87398406459845, 3.304895298656512,
              -14.138107946970608),

    TestValue(14.1, 1.35, 27.626847510527487, 3.0059807652505035,
              -10.49582397274031),

    TestValue(14.1, 5.42, 7.514275718863735, 1.6533926996855581,
              -2.7996159235447755),

    TestValue(14.1, 15.61, -10.894914609874995, 0.7957534956863926,
              -1.36543716043966),

    TestValue(31.4, 4.e-8, 631.9793629723308, 21.158333378062316, -7.85e8),

    TestValue(31.4, 2.56e-6, 501.3904341548371, 16.999450294647104,
              -1.2265625000000043e7),

    TestValue(31.4, 1, 97.0914080168896, 4.1242173963614, -31.41644277047241),

    TestValue(31.4, 1.35, 87.66136250153109, 3.824335037153727,
              -23.281451899330236),

    TestValue(31.4, 5.42, 43.790021183801784, 2.4417007517315192,
              -5.881782356655926),

    TestValue(31.4, 15.61, 8.8739384452097, 1.436044516907569,
              -2.252872461451006),
};

double allowed_rel_error = 1e-8;

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
  std::for_each(testValues.begin(), testValues.end(), [](TestValue test) {
    double f1 = log_modified_bessel_second_kind_frac(test.v, test.z);
    EXPECT_REL_ERROR(test.value, f1)
        << "ratio at v = " << test.v << ", z = " << test.z;
  });
}

TEST(AgradRev, log_modified_bessel_second_kind_frac_double_var) {
  std::for_each(testValues.begin(), testValues.end(), [](TestValue test) {
    AVAR z(test.z);
    try {
      AVAR f = log_modified_bessel_second_kind_frac(test.v, z);
      EXPECT_REL_ERROR(test.value, f.val())
          << "ratio at v = " << test.v << ", z = " << test.z;

      AVEC x = createAVEC(z);
      VEC g;
      f.grad(x, g);
      EXPECT_REL_ERROR(test.grad_z, g[0])
          << "grad_z at v = " << test.v << ", z = " << test.z;
    } catch (const std::domain_error& err) {
      std::cout << "\nAt v = " << test.v << ", z = " << test.z << ":\n";
      throw;
    }
  });
}

TEST(AgradRev, log_modified_bessel_second_kind_frac_var_double) {
  std::for_each(testValues.begin(), testValues.end(), [](TestValue test) {
    try {
      AVAR v(test.v);
      AVAR f = log_modified_bessel_second_kind_frac(v, test.z);
      EXPECT_REL_ERROR(test.value, f.val())
          << "ratio at v = " << test.v << ", z = " << test.z;

      AVEC x = createAVEC(v);
      VEC g;
      f.grad(x, g);
      EXPECT_REL_ERROR(test.grad_v, g[0])
          << "grad_v at v = " << test.v << ", z = " << test.z;
    } catch (const std::domain_error& err) {
      std::cout << "\nAt v = " << test.v << ", z = " << test.z << ":\n";
      throw;
    }
  });
}

TEST(AgradRev, log_modified_bessel_second_kind_frac_var_var) {
  std::for_each(testValues.begin(), testValues.end(), [](TestValue test) {
    try {
      AVAR v(test.v);
      AVAR z(test.z);
      AVAR f = log_modified_bessel_second_kind_frac(v, z);
      EXPECT_REL_ERROR(test.value, f.val())
          << "ratio at v = " << test.v << ", z = " << test.z;

      AVEC x = createAVEC(v, z);
      VEC g;
      f.grad(x, g);
      EXPECT_REL_ERROR(test.grad_v, g[0])
          << "grad_v at v = " << test.v << ", z = " << test.z;
      EXPECT_REL_ERROR(test.grad_z, g[1])
          << "grad_z at v = " << test.v << ", z = " << test.z;
    } catch (const std::domain_error& err) {
      std::cout << "\nAt v = " << test.v << ", z = " << test.z << ":\n";
      throw;
    }
  });
}
