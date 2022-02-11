#ifndef TEST_PROB_TEST_FIXTURE_CDF_HPP
#define TEST_PROB_TEST_FIXTURE_CDF_HPP

#include <stan/math/rev.hpp>
#include <test/prob/utility.hpp>
#include <type_traits>

using Eigen::Dynamic;
using Eigen::Matrix;
using stan::is_constant_all;
using stan::is_vector;
using stan::scalar_type;
using stan::math::value_of;
using stan::math::var;
using std::vector;

class AgradCdfTest {
 public:
  virtual void valid_values(vector<vector<double>>& /*parameters*/,
                            vector<double>& /* cdf */) {
    throw std::runtime_error("valid_values() not implemented");
  }

  // don't need to list nan. checked by the test.
  virtual void invalid_values(vector<size_t>& /*index*/,
                              vector<double>& /*value*/) {
    throw std::runtime_error("invalid_values() not implemented");
  }

  virtual bool has_lower_bound() { return false; }

  virtual double lower_bound() {
    return -std::numeric_limits<double>::infinity();
  }

  virtual bool has_upper_bound() { return false; }

  virtual double upper_bound() {
    return std::numeric_limits<double>::infinity();
  }

  // also include 2 templated functions:
  /*
    template <typename T_y, typename T_loc, typename T_scale,
    typename T3, typename T4, typename T5,
    typename T6, typename T7, typename T8,
    typename T9>
    stan::return_type_t<T_y, T_loc, T_scale>
    cdf(const T_y& y, const T_loc& mu, const T_scale& sigma,
    const T3&, const T4&, const T5&, const T6&, const T7&, const T8&, const T9&)
    { return stan::math::normal_cdf(y, mu, sigma);
    }

    template <typename T_y, typename T_loc, typename T_scale,
    typename T3, typename T4, typename T5,
    typename T6, typename T7, typename T8,
    typename T9>
    stan::return_type_t<T_y, T_loc, T_scale>
    cdf_function(const T_y& y, const T_loc& mu, const T_scale& sigma,
    const T3&, const T4&, const T5&, const T6&, const T7&, const T8&, const T9&)
    { using stan::math::erf; return (0.5 + 0.5 * erf((y - mu) / (sigma *
    SQRT_TWO)));
    }
  */
};

template <class T>
class AgradCdfTestFixture : public ::testing::Test {
 public:
  std::tuple_element_t<0, T> TestClass;
  typedef std::tuple_element_t<1, T> ArgClass;
  typedef std::tuple_element_t<0, ArgClass> T0;
  typedef std::tuple_element_t<1, ArgClass> T1;
  typedef std::tuple_element_t<2, ArgClass> T2;
  typedef std::tuple_element_t<3, ArgClass> T3;
  typedef std::tuple_element_t<4, ArgClass> T4;
  typedef std::tuple_element_t<5, ArgClass> T5;

  typedef typename scalar_type<T0>::type Scalar0;
  typedef typename scalar_type<T1>::type Scalar1;
  typedef typename scalar_type<T2>::type Scalar2;
  typedef typename scalar_type<T3>::type Scalar3;
  typedef typename scalar_type<T4>::type Scalar4;
  typedef typename scalar_type<T5>::type Scalar5;

  typedef
      typename stan::math::fvar<stan::partials_return_t<T0, T1, T2, T3, T4, T5>>
          T_fvar_return;
  using T_return_type = stan::return_type_t<T0, T1, T2, T3, T4, T5>;

  void call_all_versions() {
    vector<double> cdf;
    vector<vector<double>> parameters;
    TestClass.valid_values(parameters, cdf);

    T0 p0 = get_params<T0>(parameters, 0);
    T1 p1 = get_params<T1>(parameters, 1);
    T2 p2 = get_params<T2>(parameters, 2);
    T3 p3 = get_params<T3>(parameters, 3);
    T4 p4 = get_params<T4>(parameters, 4);
    T5 p5 = get_params<T5>(parameters, 5);

    EXPECT_NO_THROW(({
      TestClass.template cdf<T0, T1, T2, T3, T4, T5>(p0, p1, p2, p3, p4, p5);
    })) << "Calling cdf throws exception with default parameters";
  }

  void test_valid_values() {
    vector<double> expected_cdf;
    vector<vector<double>> parameters;
    TestClass.valid_values(parameters, expected_cdf);

    for (size_t n = 0; n < parameters.size(); n++) {
      T0 p0 = get_params<T0>(parameters, n, 0);
      T1 p1 = get_params<T1>(parameters, n, 1);
      T2 p2 = get_params<T2>(parameters, n, 2);
      T3 p3 = get_params<T3>(parameters, n, 3);
      T4 p4 = get_params<T4>(parameters, n, 4);
      T5 p5 = get_params<T5>(parameters, n, 5);

      T_return_type cdf(0);
      EXPECT_NO_THROW(({
        cdf = TestClass.template cdf<T0, T1, T2, T3, T4, T5>(p0, p1, p2, p3, p4,
                                                             p5);
      })) << "Valid parameters failed at index: "
          << n << " -- " << parameters[n];
      EXPECT_TRUE(cdf >= 0)
          << "cdf value must be greater than or equal to 0. cdf value: " << cdf;
      EXPECT_TRUE(cdf <= 1)
          << "cdf value must be less than or equal to 1. cdf value: " << cdf;

      if (all_scalar<T0, T1, T2, T3, T4, T5>::value) {
        EXPECT_TRUE(stan::math::abs(expected_cdf[n] - cdf) < 1e-8)
            << "For all scalar inputs cdf should match the provided value. "
               "Failed at index: "
            << n << std::endl
            << "expected: " << expected_cdf[n] << std::endl
            << "actual:   " << cdf;
      }
    }
  }

  void test_nan_value(const vector<double>& parameters, const size_t n) {
    vector<double> invalid_params(parameters);
    invalid_params[n] = std::numeric_limits<double>::quiet_NaN();

    Scalar0 p0 = get_param<Scalar0>(invalid_params, 0);
    Scalar1 p1 = get_param<Scalar1>(invalid_params, 1);
    Scalar2 p2 = get_param<Scalar2>(invalid_params, 2);
    Scalar3 p3 = get_param<Scalar3>(invalid_params, 3);
    Scalar4 p4 = get_param<Scalar4>(invalid_params, 4);
    Scalar5 p5 = get_param<Scalar5>(invalid_params, 5);

    EXPECT_THROW(({
                   TestClass.template cdf<Scalar0, Scalar1, Scalar2, Scalar3,
                                          Scalar4, Scalar5>(p0, p1, p2, p3, p4,
                                                            p5);
                 }),
                 std::domain_error)
        << "NaN value at index " << n << " should have failed" << std::endl
        << invalid_params;
  }

  void test_invalid_values() {
    vector<double> parameters = this->first_valid_params();

    vector<size_t> index;
    vector<double> invalid_values;
    TestClass.invalid_values(index, invalid_values);

    for (size_t n = 0; n < index.size(); n++) {
      vector<double> invalid_params(parameters);
      invalid_params[index[n]] = invalid_values[n];

      Scalar0 p0 = get_param<Scalar0>(invalid_params, 0);
      Scalar1 p1 = get_param<Scalar1>(invalid_params, 1);
      Scalar2 p2 = get_param<Scalar2>(invalid_params, 2);
      Scalar3 p3 = get_param<Scalar3>(invalid_params, 3);
      Scalar4 p4 = get_param<Scalar4>(invalid_params, 4);
      Scalar5 p5 = get_param<Scalar5>(invalid_params, 5);

      EXPECT_THROW(({
                     TestClass.template cdf<Scalar0, Scalar1, Scalar2, Scalar3,
                                            Scalar4, Scalar5>(p0, p1, p2, p3,
                                                              p4, p5);
                   }),
                   std::domain_error)
          << "Invalid value " << n << " should have failed" << std::endl
          << invalid_params;
    }
    if (std::numeric_limits<Scalar0>::has_quiet_NaN && parameters.size() > 0)
      test_nan_value(parameters, 0);
    if (std::numeric_limits<Scalar1>::has_quiet_NaN && parameters.size() > 1)
      test_nan_value(parameters, 1);
    if (std::numeric_limits<Scalar2>::has_quiet_NaN && parameters.size() > 2)
      test_nan_value(parameters, 2);
    if (std::numeric_limits<Scalar3>::has_quiet_NaN && parameters.size() > 3)
      test_nan_value(parameters, 3);
    if (std::numeric_limits<Scalar4>::has_quiet_NaN && parameters.size() > 4)
      test_nan_value(parameters, 4);
    if (std::numeric_limits<Scalar5>::has_quiet_NaN && parameters.size() > 5)
      test_nan_value(parameters, 5);
  }

  void add_finite_diff_1storder(const vector<double>& params,
                                vector<double>& finite_diff, const size_t n) {
    const double e = 1e-8;
    const double e2 = 2 * e;

    vector<double> plus(6);
    vector<double> minus(6);
    for (size_t i = 0; i < 6; i++) {
      plus[i] = get_param<double>(params, i);
      minus[i] = get_param<double>(params, i);
    }
    plus[n] += e;
    minus[n] -= e;

    double cdf_plus
        = TestClass.cdf(plus[0], plus[1], plus[2], plus[3], plus[4], plus[5]);
    double cdf_minus = TestClass.cdf(minus[0], minus[1], minus[2], minus[3],
                                     minus[4], minus[5]);

    finite_diff.push_back((cdf_plus - cdf_minus) / e2);
  }

  void calculate_finite_diff(const vector<double>& params,
                             vector<double>& finite_diff) {
    if (!is_constant_all<Scalar0>::value && !is_empty<Scalar0>::value)
      add_finite_diff_1storder(params, finite_diff, 0);
    if (!is_constant_all<Scalar1>::value && !is_empty<Scalar1>::value)
      add_finite_diff_1storder(params, finite_diff, 1);
    if (!is_constant_all<Scalar2>::value && !is_empty<Scalar2>::value)
      add_finite_diff_1storder(params, finite_diff, 2);
    if (!is_constant_all<Scalar3>::value && !is_empty<Scalar3>::value)
      add_finite_diff_1storder(params, finite_diff, 3);
    if (!is_constant_all<Scalar4>::value && !is_empty<Scalar4>::value)
      add_finite_diff_1storder(params, finite_diff, 4);
    if (!is_constant_all<Scalar5>::value && !is_empty<Scalar5>::value)
      add_finite_diff_1storder(params, finite_diff, 5);
  }

  // works for <double>
  template <typename... Args>
  double calculate_gradients_1storder(vector<double>& grad, double& logprob,
                                      Args&... args) {
    return logprob;
  }
  template <typename... Args>
  double calculate_gradients_2ndorder(vector<double>& grad, double& logprob,
                                      Args&... x) {
    return logprob;
  }
  template <typename... Args>
  double calculate_gradients_3rdorder(vector<double>& grad, double& logprob,
                                      Args&... x) {
    return logprob;
  }

  // works for <var>
  template <typename... Args>
  double calculate_gradients_1storder(vector<double>& grad, var& logprob,
                                      Args&... args) {
    stan::math::set_zero_all_adjoints();
    logprob.grad();
    add_adjoints(grad, args...);
    return logprob.val();
  }
  template <typename... Args>
  double calculate_gradients_2ndorder(vector<double>& grad, var& logprob,
                                      Args&... args) {
    return logprob.val();
  }
  template <typename... Args>
  double calculate_gradients_3rdorder(vector<double>& grad, var& logprob,
                                      Args&... args) {
    return logprob.val();
  }

  // works for fvar<double>
  template <typename... Args>
  double calculate_gradients_1storder(vector<double>& grad,
                                      fvar<double>& logprob, Args&... args) {
    grad.push_back(logprob.d_);
    return logprob.val();
  }
  template <typename... Args>
  double calculate_gradients_2ndorder(vector<double>& grad,
                                      fvar<double>& logprob, Args&... args) {
    return logprob.val();
  }
  template <typename... Args>
  double calculate_gradients_3rdorder(vector<double>& grad,
                                      fvar<double>& logprob, Args&... args) {
    return logprob.val();
  }

  // works for fvar<fvar<double> >
  template <typename... Args>
  double calculate_gradients_1storder(vector<double>& grad,
                                      fvar<fvar<double>>& logprob,
                                      Args&... args) {
    grad.push_back(logprob.d_.val_);

    return logprob.val().val();
  }
  template <typename... Args>
  double calculate_gradients_2ndorder(vector<double>& grad,
                                      fvar<fvar<double>>& logprob,
                                      Args&... args) {
    grad.push_back(logprob.d_.d_);
    return logprob.val().val();
  }
  template <typename... Args>
  double calculate_gradients_3rdorder(vector<double>& grad,
                                      fvar<fvar<double>>& logprob,
                                      Args&... args) {
    return logprob.val().val();
  }

  // works for fvar<var>
  template <typename... Args>
  double calculate_gradients_1storder(vector<double>& grad, fvar<var>& logprob,
                                      Args&... args) {
    stan::math::set_zero_all_adjoints();
    logprob.val_.grad();
    add_adjoints(grad, args...);
    return logprob.val_.val();
  }
  template <typename... Args>
  double calculate_gradients_2ndorder(vector<double>& grad, fvar<var>& logprob,
                                      Args&... args) {
    stan::math::set_zero_all_adjoints();
    logprob.d_.grad();
    add_adjoints(grad, args...);
    return logprob.val_.val();
  }
  template <typename... Args>
  double calculate_gradients_3rdorder(vector<double>& grad, fvar<var>& logprob,
                                      Args&... args) {
    return logprob.val_.val();
  }

  // works for fvar<fvar<var> >
  template <typename... Args>
  double calculate_gradients_1storder(vector<double>& grad,
                                      fvar<fvar<var>>& logprob, Args&... args) {
    stan::math::set_zero_all_adjoints();
    logprob.val_.val_.grad();
    add_adjoints(grad, args...);
    return logprob.val_.val_.val();
  }
  template <typename... Args>
  double calculate_gradients_2ndorder(vector<double>& grad,
                                      fvar<fvar<var>>& logprob, Args&... args) {
    stan::math::set_zero_all_adjoints();
    logprob.d_.val_.grad();
    add_adjoints(grad, args...);

    return logprob.val_.val_.val();
  }
  template <typename... Args>
  double calculate_gradients_3rdorder(vector<double>& grad,
                                      fvar<fvar<var>>& logprob, Args&... args) {
    stan::math::set_zero_all_adjoints();
    logprob.d_.d_.grad();
    add_adjoints(grad, args...);
    return logprob.val_.val_.val();
  }

  void test_finite_diffs_equal(const vector<double>& parameters,
                               const vector<double>& finite_dif,
                               const vector<double>& gradients) {
    ASSERT_EQ(finite_dif.size(), gradients.size())
        << "Number of finite diff gradients and calculated gradients must "
           "match -- error in test fixture";
    for (size_t i = 0; i < finite_dif.size(); i++) {
      EXPECT_NEAR(finite_dif[i], gradients[i], 1e-4)
          << "Comparison of finite diff to calculated gradient failed for i="
          << i << ": " << parameters << std::endl
          << "  finite diffs: " << finite_dif << std::endl
          << "  grads:        " << gradients;
    }
  }

  void test_finite_diff() {
    if (all_constant<T0, T1, T2, T3, T4, T5>::value) {
      SUCCEED() << "No test for all double arguments";
      return;
    }
    if (any_vector<T0, T1, T2, T3, T4, T5>::value) {
      SUCCEED() << "No test for vector arguments";
      return;
    }
    vector<double> expected_cdf;
    vector<vector<double>> parameters;
    TestClass.valid_values(parameters, expected_cdf);

    for (size_t n = 0; n < parameters.size(); n++) {
      vector<double> finite_diffs;
      vector<double> gradients;

      if (!std::is_same<Scalar0, fvar<double>>::value
          && !std::is_same<Scalar0, fvar<fvar<double>>>::value
          && !std::is_same<Scalar1, fvar<double>>::value
          && !std::is_same<Scalar1, fvar<fvar<double>>>::value
          && !std::is_same<Scalar2, fvar<double>>::value
          && !std::is_same<Scalar2, fvar<fvar<double>>>::value
          && !std::is_same<Scalar3, fvar<double>>::value
          && !std::is_same<Scalar3, fvar<fvar<double>>>::value
          && !std::is_same<Scalar4, fvar<double>>::value
          && !std::is_same<Scalar4, fvar<fvar<double>>>::value
          && !std::is_same<Scalar5, fvar<double>>::value
          && !std::is_same<Scalar5, fvar<fvar<double>>>::value) {
        calculate_finite_diff(parameters[n], finite_diffs);
        Scalar0 p0 = get_param<Scalar0>(parameters[n], 0);
        Scalar1 p1 = get_param<Scalar1>(parameters[n], 1);
        Scalar2 p2 = get_param<Scalar2>(parameters[n], 2);
        Scalar3 p3 = get_param<Scalar3>(parameters[n], 3);
        Scalar4 p4 = get_param<Scalar4>(parameters[n], 4);
        Scalar5 p5 = get_param<Scalar5>(parameters[n], 5);

        T_return_type cdf
            = TestClass.template cdf<Scalar0, Scalar1, Scalar2, Scalar3,
                                     Scalar4, Scalar5>(p0, p1, p2, p3, p4, p5);

        calculate_gradients_1storder(gradients, cdf, p0, p1, p2, p3, p4, p5);

        test_finite_diffs_equal(parameters[n], finite_diffs, gradients);

        stan::math::recover_memory();
      }
    }
  }

  void test_gradients_equal(const vector<double>& expected_gradients,
                            const vector<double>& gradients) {
    ASSERT_EQ(expected_gradients.size(), gradients.size())
        << "Number of expected gradients and calculated gradients must match "
           "-- error in test fixture";
    for (size_t i = 0; i < expected_gradients.size(); i++) {
      EXPECT_NEAR(expected_gradients[i], gradients[i], 1e-4)
          << "Comparison of expected gradient to calculated gradient failed";
    }
  }

  void test_gradient_function() {
    if (all_constant<T0, T1, T2, T3, T4, T5>::value) {
      SUCCEED() << "No test for all double arguments";
      return;
    }
    if (any_vector<T0, T1, T2, T3, T4, T5>::value) {
      SUCCEED() << "No test for vector arguments";
      return;
    }
    vector<double> expected_cdf;
    vector<vector<double>> parameters;
    TestClass.valid_values(parameters, expected_cdf);

    for (size_t n = 0; n < parameters.size(); n++) {
      vector<double> expected_gradients1;
      vector<double> expected_gradients2;
      vector<double> expected_gradients3;
      vector<double> gradients1;
      vector<double> gradients2;
      vector<double> gradients3;

      Scalar0 p0 = get_param<Scalar0>(parameters[n], 0);
      Scalar1 p1 = get_param<Scalar1>(parameters[n], 1);
      Scalar2 p2 = get_param<Scalar2>(parameters[n], 2);
      Scalar3 p3 = get_param<Scalar3>(parameters[n], 3);
      Scalar4 p4 = get_param<Scalar4>(parameters[n], 4);
      Scalar5 p5 = get_param<Scalar5>(parameters[n], 5);

      T_return_type cdf
          = TestClass.template cdf<Scalar0, Scalar1, Scalar2, Scalar3, Scalar4,
                                   Scalar5>(p0, p1, p2, p3, p4, p5);

      T_return_type cdf_funct
          = TestClass.template cdf_function<Scalar0, Scalar1, Scalar2, Scalar3,
                                            Scalar4, Scalar5>(p0, p1, p2, p3,
                                                              p4, p5);

      calculate_gradients_1storder(expected_gradients1, cdf_funct, p0, p1, p2,
                                   p3, p4, p5);
      calculate_gradients_1storder(gradients1, cdf, p0, p1, p2, p3, p4, p5);
      calculate_gradients_2ndorder(expected_gradients2, cdf_funct, p0, p1, p2,
                                   p3, p4, p5);
      calculate_gradients_2ndorder(gradients2, cdf, p0, p1, p2, p3, p4, p5);
      calculate_gradients_3rdorder(expected_gradients3, cdf_funct, p0, p1, p2,
                                   p3, p4, p5);
      calculate_gradients_3rdorder(gradients3, cdf, p0, p1, p2, p3, p4, p5);

      test_gradients_equal(expected_gradients1, gradients1);
      test_gradients_equal(expected_gradients2, gradients2);
      test_gradients_equal(expected_gradients3, gradients3);

      stan::math::recover_memory();
    }
  }

  /**
   * Test that the vectorized functions work as expected when the elements
   * of the vector are the same
   *
   * For cdf this means cdf([a, a, a]) == cdf(a) * cdf(a) * cdf(a)
   */
  void test_repeat_as_vector() {
    using stan::math::pow;

    if (stan::is_any_var_matrix<T0, T1, T2, T3, T4, T5>::value) {
      // There is no way to do this test for a `var_value` matrix
      // because this is testing what happens when all elements of
      // the vector are the same thing. This works the PIMP var types
      // stored in other containers because every var can point at the
      // same vari. However, when a var_value of length N is allocated
      // there are N values and N adjoints and they are all separate.

      SUCCEED() << "No test for var_value<Eigen::Matrix<T, R, C>> arguments";
      return;
    }
    if (!any_vector<T0, T1, T2, T3, T4, T5>::value) {
      SUCCEED() << "No test for non-vector arguments";
      return;
    }

    const size_t N_REPEAT = 3;
    vector<double> expected_cdf;
    vector<vector<double>> parameters;
    TestClass.valid_values(parameters, expected_cdf);

    for (size_t n = 0; n < parameters.size(); n++) {
      vector<double> single_gradients1;
      vector<double> single_gradients2;
      vector<double> single_gradients3;

      Scalar0 p0_ = get_param<Scalar0>(parameters[n], 0);
      Scalar1 p1_ = get_param<Scalar1>(parameters[n], 1);
      Scalar2 p2_ = get_param<Scalar2>(parameters[n], 2);
      Scalar3 p3_ = get_param<Scalar3>(parameters[n], 3);
      Scalar4 p4_ = get_param<Scalar4>(parameters[n], 4);
      Scalar5 p5_ = get_param<Scalar5>(parameters[n], 5);
      std::vector<Scalar0> p0s_((is_vector<T0>::value) ? N_REPEAT : 1, p0_);
      std::vector<Scalar1> p1s_((is_vector<T1>::value) ? N_REPEAT : 1, p1_);
      std::vector<Scalar2> p2s_((is_vector<T2>::value) ? N_REPEAT : 1, p2_);
      std::vector<Scalar3> p3s_((is_vector<T3>::value) ? N_REPEAT : 1, p3_);
      std::vector<Scalar4> p4s_((is_vector<T4>::value) ? N_REPEAT : 1, p4_);
      std::vector<Scalar5> p5s_((is_vector<T5>::value) ? N_REPEAT : 1, p5_);

      T_return_type cdf
          = pow(TestClass.cdf(p0_, p1_, p2_, p3_, p4_, p5_), N_REPEAT);

      double single_cdf = calculate_gradients_1storder(
          single_gradients1, cdf, p0s_, p1s_, p2s_, p3s_, p4s_, p5s_);
      calculate_gradients_2ndorder(single_gradients2, cdf, p0s_, p1s_, p2s_,
                                   p3s_, p4s_, p5s_);
      calculate_gradients_3rdorder(single_gradients3, cdf, p0s_, p1s_, p2s_,
                                   p3s_, p4s_, p5s_);

      T0 p0 = get_repeated_params<T0>(parameters[n], 0, N_REPEAT);
      T1 p1 = get_repeated_params<T1>(parameters[n], 1, N_REPEAT);
      T2 p2 = get_repeated_params<T2>(parameters[n], 2, N_REPEAT);
      T3 p3 = get_repeated_params<T3>(parameters[n], 3, N_REPEAT);
      T4 p4 = get_repeated_params<T4>(parameters[n], 4, N_REPEAT);
      T5 p5 = get_repeated_params<T5>(parameters[n], 5, N_REPEAT);

      T_return_type multiple_cdf = TestClass.cdf(p0, p1, p2, p3, p4, p5);
      vector<double> multiple_gradients1;
      vector<double> multiple_gradients2;
      vector<double> multiple_gradients3;

      calculate_gradients_1storder(multiple_gradients1, multiple_cdf, p0, p1,
                                   p2, p3, p4, p5);
      calculate_gradients_2ndorder(multiple_gradients2, multiple_cdf, p0, p1,
                                   p2, p3, p4, p5);
      calculate_gradients_3rdorder(multiple_gradients3, multiple_cdf, p0, p1,
                                   p2, p3, p4, p5);

      stan::math::recover_memory();

      stan::test::expect_near_rel(
          "cdf with repeated vector input should match a multiple of cdf of "
          "single input",
          stan::math::value_of_rec(single_cdf),
          stan::math::value_of_rec(multiple_cdf));

      stan::test::expect_near_rel(
          "scalar and vectorized results should have the same first order "
          "gradients",
          single_gradients1, multiple_gradients1);
      stan::test::expect_near_rel(
          "scalar and vectorized results should have the same second order "
          "gradients",
          single_gradients2, multiple_gradients2);
      stan::test::expect_near_rel(
          "scalar and vectorized results should have the same third order "
          "gradients",
          single_gradients3, multiple_gradients3);
    }
  }

  /**
   * Test that the vectorized functions work as expected when the elements
   * of the vector are different
   *
   * For cdf this means cdf([a, b, c]) == cdf(a) * cdf(b) * cdf(c)
   */
  void test_as_scalars_vs_as_vector() {
    if (!any_vector<T0, T1, T2, T3, T4, T5>::value) {
      SUCCEED() << "No test for non-vector arguments";
      return;
    }
    vector<double> expected_cdf;
    vector<vector<double>> parameters;
    TestClass.valid_values(parameters, expected_cdf);

    vector<double> single_gradients1;
    vector<double> single_gradients2;
    vector<double> single_gradients3;

    vector<double> multiple_gradients1;
    vector<double> multiple_gradients2;
    vector<double> multiple_gradients3;

    T0 p0 = get_params<T0>(parameters, 0);
    T1 p1 = get_params<T1>(parameters, 1);
    T2 p2 = get_params<T2>(parameters, 2);
    T3 p3 = get_params<T3>(parameters, 3);
    T4 p4 = get_params<T4>(parameters, 4);
    T5 p5 = get_params<T5>(parameters, 5);

    vector<Scalar0> p0s = {get_param<Scalar0>(parameters[0], 0)};
    vector<Scalar1> p1s = {get_param<Scalar1>(parameters[0], 1)};
    vector<Scalar2> p2s = {get_param<Scalar2>(parameters[0], 2)};
    vector<Scalar3> p3s = {get_param<Scalar3>(parameters[0], 3)};
    vector<Scalar4> p4s = {get_param<Scalar4>(parameters[0], 4)};
    vector<Scalar5> p5s = {get_param<Scalar5>(parameters[0], 5)};

    T_return_type single_cdf = TestClass.cdf(
        p0s.back(), p1s.back(), p2s.back(), p3s.back(), p4s.back(), p5s.back());

    for (size_t n = 1; n < parameters.size(); n++) {
      if (is_vector<T0>::value)
        p0s.push_back(get_param<Scalar0>(parameters[n], 0));

      if (is_vector<T1>::value)
        p1s.push_back(get_param<Scalar1>(parameters[n], 1));

      if (is_vector<T2>::value)
        p2s.push_back(get_param<Scalar2>(parameters[n], 2));

      if (is_vector<T3>::value)
        p3s.push_back(get_param<Scalar3>(parameters[n], 3));

      if (is_vector<T4>::value)
        p4s.push_back(get_param<Scalar4>(parameters[n], 4));

      if (is_vector<T5>::value)
        p5s.push_back(get_param<Scalar5>(parameters[n], 5));

      single_cdf *= TestClass.cdf(p0s.back(), p1s.back(), p2s.back(),
                                  p3s.back(), p4s.back(), p5s.back());
    }

    calculate_gradients_1storder(single_gradients1, single_cdf, p0s, p1s, p2s,
                                 p3s, p4s, p5s);
    calculate_gradients_2ndorder(single_gradients2, single_cdf, p0s, p1s, p2s,
                                 p3s, p4s, p5s);
    calculate_gradients_3rdorder(single_gradients3, single_cdf, p0s, p1s, p2s,
                                 p3s, p4s, p5s);

    T_return_type multiple_cdf = TestClass.cdf(p0, p1, p2, p3, p4, p5);

    calculate_gradients_1storder(multiple_gradients1, multiple_cdf, p0, p1, p2,
                                 p3, p4, p5);
    calculate_gradients_2ndorder(multiple_gradients2, multiple_cdf, p0, p1, p2,
                                 p3, p4, p5);
    calculate_gradients_3rdorder(multiple_gradients3, multiple_cdf, p0, p1, p2,
                                 p3, p4, p5);

    stan::math::recover_memory();

    if (stan::math::is_inf(stan::math::value_of_rec(single_cdf))
        && stan::math::value_of_rec(single_cdf)
               == stan::math::value_of_rec(multiple_cdf)) {
      return;
    }

    stan::test::expect_near_rel(
        "product of scalar cdfs should match vectorized result",
        stan::math::value_of_rec(single_cdf),
        stan::math::value_of_rec(multiple_cdf));

    stan::test::expect_near_rel(
        "scalar and vectorized results should have the same first order "
        "gradients",
        single_gradients1, multiple_gradients1);
    stan::test::expect_near_rel(
        "scalar and vectorized results should have the same second order "
        "gradients",
        single_gradients2, multiple_gradients2);
    stan::test::expect_near_rel(
        "scalar and vectorized results should have the same third order "
        "gradients",
        single_gradients3, multiple_gradients3);
  }

  void test_lower_bound() {
    using stan::math::value_of;
    const size_t N_REPEAT = 3;
    vector<double> expected_cdf;
    vector<vector<double>> parameters;
    TestClass.valid_values(parameters, expected_cdf);

    if (!TestClass.has_lower_bound()) {
      if (!std::numeric_limits<Scalar0>::has_infinity) {
        for (size_t n = 0; n < parameters.size(); n++)
          parameters[n][0] = value_of(
              value_of(value_of(std::numeric_limits<Scalar0>::min())));
      } else {
        for (size_t n = 0; n < parameters.size(); n++)
          parameters[n][0] = -std::numeric_limits<double>::infinity();
      }
    } else {
      for (size_t n = 0; n < parameters.size(); n++)
        parameters[n][0] = TestClass.lower_bound();
    }

    for (size_t n = 0; n < parameters.size(); n++) {
      T0 p0 = get_repeated_params<T0>(parameters[n], 0, N_REPEAT);
      T1 p1 = get_repeated_params<T1>(parameters[n], 1, N_REPEAT);
      T2 p2 = get_repeated_params<T2>(parameters[n], 2, N_REPEAT);
      T3 p3 = get_repeated_params<T3>(parameters[n], 3, N_REPEAT);
      T4 p4 = get_repeated_params<T4>(parameters[n], 4, N_REPEAT);
      T5 p5 = get_repeated_params<T5>(parameters[n], 5, N_REPEAT);

      T_return_type cdf_at_lower_bound
          = TestClass.template cdf<T0, T1, T2, T3, T4, T5>(p0, p1, p2, p3, p4,
                                                           p5);
      EXPECT_TRUE(0.0 == cdf_at_lower_bound)
          << "CDF evaluated at lower bound should equal 0";
    }
  }

  void test_upper_bound() {
    using stan::math::value_of;
    const size_t N_REPEAT = 3;
    vector<double> expected_cdf;
    vector<vector<double>> parameters;
    TestClass.valid_values(parameters, expected_cdf);

    if (!TestClass.has_upper_bound()) {
      if (!std::numeric_limits<Scalar0>::has_infinity) {
        for (size_t n = 0; n < parameters.size(); n++)
          parameters[n][0] = value_of(
              value_of(value_of(std::numeric_limits<Scalar0>::max())));
      } else {
        for (size_t n = 0; n < parameters.size(); n++)
          parameters[n][0] = std::numeric_limits<double>::infinity();
      }
    } else {
      for (size_t n = 0; n < parameters.size(); n++)
        parameters[n][0] = TestClass.upper_bound();
    }

    for (size_t n = 0; n < parameters.size(); n++) {
      T0 p0 = get_repeated_params<T0>(parameters[n], 0, N_REPEAT);
      T1 p1 = get_repeated_params<T1>(parameters[n], 1, N_REPEAT);
      T2 p2 = get_repeated_params<T2>(parameters[n], 2, N_REPEAT);
      T3 p3 = get_repeated_params<T3>(parameters[n], 3, N_REPEAT);
      T4 p4 = get_repeated_params<T4>(parameters[n], 4, N_REPEAT);
      T5 p5 = get_repeated_params<T5>(parameters[n], 5, N_REPEAT);

      T_return_type cdf_at_upper_bound
          = TestClass.template cdf<T0, T1, T2, T3, T4, T5>(p0, p1, p2, p3, p4,
                                                           p5);
      EXPECT_TRUE(1.0 == cdf_at_upper_bound)
          << "CDF evaluated at upper bound is " << cdf_at_upper_bound
          << " but should equal 1";
    }
  }

  void test_length_0_vector() {
    if (!any_vector<T0, T1, T2, T3, T4, T5>::value) {
      SUCCEED() << "No test for non-vector arguments";
      return;
    }
    const size_t N_REPEAT = 0;
    vector<double> expected_cdf;
    vector<vector<double>> parameters;
    TestClass.valid_values(parameters, expected_cdf);

    T0 p0 = get_repeated_params<T0>(parameters[0], 0, N_REPEAT);
    T1 p1 = get_repeated_params<T1>(parameters[0], 1, N_REPEAT);
    T2 p2 = get_repeated_params<T2>(parameters[0], 2, N_REPEAT);
    T3 p3 = get_repeated_params<T3>(parameters[0], 3, N_REPEAT);
    T4 p4 = get_repeated_params<T4>(parameters[0], 4, N_REPEAT);
    T5 p5 = get_repeated_params<T5>(parameters[0], 5, N_REPEAT);

    T_return_type cdf = TestClass.template cdf<T0, T1, T2, T3, T4, T5>(
        p0, p1, p2, p3, p4, p5);

    EXPECT_TRUE(1.0 == cdf) << "cdf with an empty vector should return 1.0";
  }

  vector<double> first_valid_params() {
    vector<vector<double>> params;
    vector<double> cdf;

    TestClass.valid_values(params, cdf);
    return params[0];
  }
};

TYPED_TEST_SUITE_P(AgradCdfTestFixture);

TYPED_TEST_P(AgradCdfTestFixture, CallAllVersions) {
  this->call_all_versions();
}

TYPED_TEST_P(AgradCdfTestFixture, ValidValues) { this->test_valid_values(); }

TYPED_TEST_P(AgradCdfTestFixture, InvalidValues) {
  this->test_invalid_values();
}

TYPED_TEST_P(AgradCdfTestFixture, FiniteDiff) { this->test_finite_diff(); }

TYPED_TEST_P(AgradCdfTestFixture, Function) { this->test_gradient_function(); }

TYPED_TEST_P(AgradCdfTestFixture, RepeatAsVector) {
  this->test_repeat_as_vector();
}

TYPED_TEST_P(AgradCdfTestFixture, AsScalarsVsAsVector) {
  this->test_as_scalars_vs_as_vector();
}

TYPED_TEST_P(AgradCdfTestFixture, LowerBound) { this->test_lower_bound(); }

TYPED_TEST_P(AgradCdfTestFixture, UpperBound) { this->test_upper_bound(); }

TYPED_TEST_P(AgradCdfTestFixture, Length0Vector) {
  this->test_length_0_vector();
}

REGISTER_TYPED_TEST_SUITE_P(AgradCdfTestFixture, CallAllVersions, ValidValues,
                            InvalidValues, FiniteDiff, Function, RepeatAsVector,
                            LowerBound, UpperBound, AsScalarsVsAsVector,
                            Length0Vector);

GTEST_ALLOW_UNINSTANTIATED_PARAMETERIZED_TEST(AgradCdfTestFixture);

#endif
