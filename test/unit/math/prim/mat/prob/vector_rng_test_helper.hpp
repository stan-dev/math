#ifndef TEST_UNIT_MATH_PRIM_MAT_PROB_VECTOR_RNG_TEST_HELPER_HPP
#define TEST_UNIT_MATH_PRIM_MAT_PROB_VECTOR_RNG_TEST_HELPER_HPP

#include <gtest/gtest.h>
#include <boost/math/distributions.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <stan/math/prim/mat.hpp>
#include <test/unit/math/prim/scal/prob/util.hpp>
#include <test/unit/math/prim/scal/meta/apply_template_permutations.hpp>
#include <test/unit/math/prim/mat/prob/VectorRNGTestRig.hpp>
#include <algorithm>
#include <tuple>
#include <vector>

/*
 * Fill the vector-like variable params with values from the values argument.
 *
 * values can be shorter than params (in which cause multiple copies of values
 * are tiled over the params vector)
 *
 * @tparam T_param Type of params vector
 * @param params Parameter vector to write values to
 * @param params Values to copy into params
 */
template <typename T_param>
void assign_parameter_values(T_param& params,
                             const std::vector<double>& values) {
  if (values.size() == 0)
    return;

  for (int i = 0; i < params.size(); i++) {
    params[i] = values[i % values.size()];
  }
}

/*
 * Fill the vector params with values from the values argument.
 *
 * values can be shorter than params (in which cause multiple copies of values
 * are tiled over the params vector)
 *
 * @param params Parameter vector to write values to
 * @param params Values to copy into params
 */
void assign_parameter_values(std::vector<double>& params,
                             const std::vector<double>& values) {
  if (values.size() == 0)
    return;

  for (size_t i = 0; i < params.size(); i++) {
    params[i] = values[i % values.size()];
  }
}

/*
 * Fill the vector params with values from the values argument.
 *
 * values can be shorter than params (in which cause multiple copies of values
 * are tiled over the params vector)
 *
 * @param params Parameter vector to write values to
 * @param params Values to copy into params
 */
void assign_parameter_values(std::vector<int>& params,
                             const std::vector<int>& values) {
  if (values.size() == 0)
    return;

  for (size_t i = 0; i < params.size(); i++) {
    params[i] = values[i % values.size()];
  }
}

/*
 * Assign param the first value of values
 *
 * @param param Output parameter to write value to
 * @param params Vector with value to copy into param
 */
void assign_parameter_values(double& param, const std::vector<double>& values) {
  if (values.size() == 0)
    return;

  param = values[0];
}

/*
 * Assign param the first value of values
 *
 * @param param Output parameter to write value to
 * @param params Vector with value to copy into param
 */
void assign_parameter_values(int& param, const std::vector<int>& values) {
  if (values.size() == 0)
    return;

  param = values[0];
}

/*
 * Resize v to be length N
 *
 * @tparam T Type of v
 * @param v Variable to resize
 * @param N New size
 */
template<typename T>
void resize_if_vector(T& v, int N) {
  v.resize(N);
}

/*
 * For doubles, resize_if_vector does nothing
 */
template<>
void resize_if_vector(double& v, int N) {
}

/*
 * For ints, resize_if_vector does nothing
 */
template<>
void resize_if_vector(int& v, int N) {
}

/*
 * check_dist_throws feeds rig.generate_samples various
 * combinations of valid and invalid parameters (as defined by good_p1_, bad_p1,
 * good_p2_, bad_p2_, good_p3_, and bad_p3_). For all combinations of valid
 * (good) parameters, generate_samples should throw no errors. For all
 * combinations with an invalid (bad) parameter, generate_samples should throw
 * domain_errors.
 *
 * If rig.good_p2_ or rig.good_p3_ are empty, then it is assumed that those parameters
 * are unused and will not be tested.
 *
 * rig.generate_samples will also be passed various other guaranteed invalid
 * values like positive infinity, negative infinity, and NaNs (these should also
 * cause domain_errors).
 *
 * rig.generated_samples is also tested to reject incompatibly sized input
 * vectors. These should cause invalid_argument errors.
 *
 * @tparam T_param1 Type of first parameter
 * @tparam T_param2 Type of second parameter
 * @tparam T_param3 Type of third parameter
 * @tparam T_rig Type of test rig for random number generator
 * @param T_rig Test rig for random number generator
 */
struct check_dist_throws {
  template<typename T_param1, typename T_param2, typename T_param3,
           typename T_rig>
  void operator()(const T_rig& rig) const {
    boost::random::mt19937 rng;

    T_param1 p1;
    T_param2 p2;
    T_param3 p3;

    using T_scalar_param1 = typename stan::scalar_type<T_param1>::type;
    using T_scalar_param2 = typename stan::scalar_type<T_param2>::type;
    using T_scalar_param3 = typename stan::scalar_type<T_param3>::type;

    bool p2_is_used = rig.p2_is_used();
    bool p3_is_used = rig.p3_is_used();

    resize_if_vector(p1, 5);  // No-op if p1 is a scalar
    resize_if_vector(p2, 5);  // No-op if p2 is a scalar
    resize_if_vector(p3, 5);  // No-op if p3 is a scalar

    // Make copies of the input arguments so that we can randomly shuffle them
    // in the tests
    std::vector<T_scalar_param1> good_p1 =
      rig.template get_good_p1<T_scalar_param1>();
    std::vector<T_scalar_param1> bad_p1 =
      rig.template get_bad_p1<T_scalar_param1>();
    std::vector<T_scalar_param2> good_p2 =
      rig.template get_good_p2<T_scalar_param2>();
    std::vector<T_scalar_param2> bad_p2 =
      rig.template get_bad_p2<T_scalar_param2>();
    std::vector<T_scalar_param3> good_p3 =
      rig.template get_good_p3<T_scalar_param3>();
    std::vector<T_scalar_param3> bad_p3 =
      rig.template get_bad_p3<T_scalar_param3>();

    // Try a few combinations of parameters that should work
    for (int i = 0; i < 5; i++) {
      std::random_shuffle(good_p1.begin(), good_p1.end());
      std::random_shuffle(good_p2.begin(), good_p2.end());
      std::random_shuffle(good_p3.begin(), good_p3.end());
      assign_parameter_values(p1, good_p1);
      assign_parameter_values(p2, good_p2);
      assign_parameter_values(p3, good_p3);
      EXPECT_NO_THROW(rig.generate_samples(p1, p2, p3, rng));
    }

    // Now try putting incompatible values in first parameter
    for (auto bad_p1_value : bad_p1) {
      assign_parameter_values(p1, { bad_p1_value });
      assign_parameter_values(p2, good_p2);
      assign_parameter_values(p3, good_p3);
      EXPECT_THROW(rig.generate_samples(p1, p2, p3, rng), std::domain_error);
    }

    // Now try putting incompatible values in second parameter
    if (p2_is_used) {
      for (auto bad_p2_value : bad_p2) {
        assign_parameter_values(p1, good_p1);
        assign_parameter_values(p2, { bad_p2_value });
        assign_parameter_values(p3, good_p3);
        EXPECT_THROW(rig.generate_samples(p1, p2, p3, rng), std::domain_error);
      }
    }

    // Now try putting incompatible values in third parameter
    if (p3_is_used) {
      for (auto bad_p3_value : bad_p3) {
        assign_parameter_values(p1, good_p1);
        assign_parameter_values(p2, good_p2);
        assign_parameter_values(p3, { bad_p3_value });
        EXPECT_THROW(rig.generate_samples(p1, p2, p3, rng), std::domain_error);
      }
    }

    // Check to make sure _rng rejects vector arguments of different lengths for
    // all parameter pairs.

    // If p1 is a scalar or the only vector, this test is skipped
    resize_if_vector(p1, 3);  // No-op if p1 is a scalar
    resize_if_vector(p2, 4);  // No-op if p2 is a scalar
    resize_if_vector(p3, 4);  // No-op if p3 is a scalar
    if (stan::length(p1) != 1 &&
        ((p2_is_used && stan::length(p2) != 1) ||
         (p3_is_used && stan::length(p3) != 1))) {
      assign_parameter_values(p1, good_p1);
      assign_parameter_values(p2, good_p2);
      assign_parameter_values(p3, good_p3);
      EXPECT_THROW(rig.generate_samples(p1, p2, p3, rng),
                   std::invalid_argument);
    }

    // If p2 is a scalar or the only vector, this test is skipped
    resize_if_vector(p1, 4);  // No-op if p1 is a scalar
    resize_if_vector(p2, 3);  // No-op if p2 is a scalar
    resize_if_vector(p3, 4);  // No-op if p3 is a scalar
    if (p2_is_used && stan::length(p2) != 1 &&
        (stan::length(p1) != 1 || (p3_is_used && stan::length(p3) != 1))) {
      assign_parameter_values(p1, good_p1);
      assign_parameter_values(p2, good_p2);
      assign_parameter_values(p3, good_p3);
      EXPECT_THROW(rig.generate_samples(p1, p2, p3, rng),
                   std::invalid_argument);
    }

    // If p3 is a scalar or the only vector, this test is skipped
    resize_if_vector(p1, 4);  // No-op if p1 is a scalar
    resize_if_vector(p2, 4);  // No-op if p2 is a scalar
    resize_if_vector(p3, 3);  // No-op if p3 is a scalar
    if (p3_is_used && stan::length(p3) != 1 &&
        (stan::length(p1) != 1 || (p2_is_used && stan::length(p2) != 1))) {
      assign_parameter_values(p1, good_p1);
      assign_parameter_values(p2, good_p2);
      assign_parameter_values(p3, good_p3);
      EXPECT_THROW(rig.generate_samples(p1, p2, p3, rng),
                   std::invalid_argument);
    }
  }
};

/*
 * This function calls check_dist_throws with the given test_rig
 * and all combinations of double, std::vector<double>, Eigen::VectorXd,
 * and Eigen::RowVectorXd as template arguments
 *
 * @tparam T_rig Test rig type for random number generator
 * @param T_rig Test rig for random number generator
 */
template<typename T_rig>
void check_dist_throws_all_types(const T_rig& rig) {
  using Eigen::VectorXd;
  using Eigen::RowVectorXd;

  apply_template_permutations<std::tuple<int, std::vector<int>, double,
                                         std::vector<double>, VectorXd,
                                         RowVectorXd> >
    (check_dist_throws{}, rig);
}

/*
 * Convert a scalar to a length one vector
 *
 * @tparam T Scalar type
 * @param v Input scalar
 * @return vector of length 1 with value v
 */
template<typename T>
std::vector<T> promote_to_vector(T v) {
  return std::vector<T>(1, v);
}

/*
 * For arguments that are already vectors, just copy. Probably more efficient
 * just using std::move but cpplint complained about use of unapproved Rvalues.
 */
template<typename T>
std::vector<T> promote_to_vector(std::vector<T> v) {
  return v;
}

/*
 * check_quantiles uses assert_matches_quantiles to check random numbers
 * generated by rig.generate_samples against the reference quantiles computed
 * with rig.generate_quantiles.
 *
 * If rig.good_p2_ or rig.good_p3_ are empty, then it is assumed that those
 * parameters are unused and will not be tested.
 *
 * @tparam T_param1 Type of first parameter
 * @tparam T_param2 Type of second parameter
 * @tparam T_param3 Type of third parameter
 * @tparam T_rig Type of test rig for random number generator
 * @param T_rig Test rig for random number generator
 */
struct check_quantiles {
  template<typename T_param1, typename T_param2, typename T_param3,
           typename T_rig>
  void operator()(const T_rig& rig) const {
    boost::random::mt19937 rng;
    T_param1 p1;
    T_param2 p2;
    T_param3 p3;
    resize_if_vector(p1, rig.M_);  // No-op if p1 is scalar
    resize_if_vector(p2, rig.M_);  // No-op if p2 is scalar
    resize_if_vector(p3, rig.M_);  // No-op if p3 is scalar

    assign_parameter_values(p1, rig.template get_good_p1
                            <typename stan::scalar_type<T_param1>::type>());
    assign_parameter_values(p2, rig.template get_good_p2
                            <typename stan::scalar_type<T_param2>::type>());
    assign_parameter_values(p3, rig.template get_good_p3
                            <typename stan::scalar_type<T_param3>::type>());

    bool p2_is_used = rig.p2_is_used();
    bool p3_is_used = rig.p3_is_used();

    int M = std::max({ stan::length(p1),
          (p2_is_used) ? stan::length(p2) : 1,
          (p3_is_used) ? stan::length(p3) : 1 });

    stan::scalar_seq_view<T_param1> p1_vec(p1);
    stan::scalar_seq_view<T_param2> p2_vec(p2);
    stan::scalar_seq_view<T_param3> p3_vec(p3);

    std::vector<std::vector<double> > samples_to_test_transpose;
    for (int n = 0; n < rig.N_; ++n) {
      // If p1, p2, and p3 are scalars, the output is a scalar. Need to promote
      // it to a std::vector
      samples_to_test_transpose.
        push_back(promote_to_vector(rig.generate_samples(p1, p2, p3, rng)));
    }

    for (int m = 0; m < M; ++m) {
      std::vector<double> samples_to_test;
      for (int n = 0; n < rig.N_; ++n) {
        samples_to_test.push_back(samples_to_test_transpose[n][m]);
      }
      std::vector<double> quantiles =
        rig.generate_quantiles(p1_vec[m], p2_vec[m], p3_vec[m]);

      assert_matches_quantiles(samples_to_test, quantiles, 1e-6);
    }
  }
};

/*
 * Call check_quantiles on the given test rig with all combinations of
 * double, std::vector<double>, Eigen::VectorXd, and Eigen::RowVectorXd
 *
 * @tparam T_rig Type of test rig for random number generator
 * @param T_rig Test rig for random number generator
 */
template<typename T_rig>
void check_quantiles_all_types(const T_rig& rig) {
  using Eigen::VectorXd;
  using Eigen::RowVectorXd;

  apply_template_permutations<std::tuple<int, double, std::vector<int>,
                                         std::vector<double>, VectorXd,
                                         RowVectorXd > >
    (check_quantiles{}, rig);
}

#endif
