#ifndef TEST_UNIT_MATH_PRIM_PROB_VECTOR_RNG_TEST_HELPER_HPP
#define TEST_UNIT_MATH_PRIM_PROB_VECTOR_RNG_TEST_HELPER_HPP

#include <gtest/gtest.h>
#include <boost/math/distributions.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <stan/math/prim.hpp>
#include <test/unit/math/prim/meta/apply_template_permutations.hpp>
#include <test/unit/math/prim/prob/util.hpp>
#include <test/unit/math/prim/prob/VectorRealRNGTestRig.hpp>
#include <test/unit/math/prim/prob/VectorIntRNGTestRig.hpp>
#include <algorithm>
#include <map>
#include <random>
#include <tuple>
#include <vector>

namespace internal {
template <class C>
void shuffle_container(C& x) {
  std::random_device rng;
  std::mt19937 twister(rng());
  std::shuffle(x.begin(), x.end(), twister);
}
}  // namespace internal

using ArgumentTypes
    = std::tuple<int, double, std::vector<int>, std::vector<double>,
                 Eigen::VectorXd, Eigen::RowVectorXd>;

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
template <typename T>
void resize_if_vector(T& v, int N) {
  v.resize(N);
}

/*
 * For doubles, resize_if_vector does nothing
 */
template <>
void resize_if_vector(double& v, int N) {}

/*
 * For ints, resize_if_vector does nothing
 */
template <>
void resize_if_vector(int& v, int N) {}

/*
 * check_dist_throws feeds rig.generate_samples various
 * combinations of valid and invalid parameters (as defined by good_p1_, bad_p1,
 * good_p2_, bad_p2_, good_p3_, and bad_p3_). For all combinations of valid
 * (good) parameters, generate_samples should throw no errors. For all
 * combinations with an invalid (bad) parameter, generate_samples should throw
 * domain_errors.
 *
 * If rig.good_p1_, rig.good_p2_ or rig.good_p3_ are empty, then it is assumed
 * that those parameters are unused and will not be tested.
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
  template <typename T_param1, typename T_param2, typename T_param3,
            typename T_rig>
  void operator()(const T_rig& rig) const {
    boost::random::mt19937 rng;

    T_param1 p1;
    T_param2 p2;
    T_param3 p3;

    using T_scalar_param1 = typename stan::scalar_type<T_param1>::type;
    using T_scalar_param2 = typename stan::scalar_type<T_param2>::type;
    using T_scalar_param3 = typename stan::scalar_type<T_param3>::type;

    bool p1_is_used = rig.p1_is_used();
    bool p2_is_used = rig.p2_is_used();
    bool p3_is_used = rig.p3_is_used();

    resize_if_vector(p1, 5);  // No-op if p1 is a scalar
    resize_if_vector(p2, 5);  // No-op if p2 is a scalar
    resize_if_vector(p3, 5);  // No-op if p3 is a scalar

    // Make copies of the input arguments so that we can randomly shuffle them
    // in the tests
    std::vector<T_scalar_param1> good_p1
        = rig.template get_good_p1<T_scalar_param1>();
    std::vector<T_scalar_param1> bad_p1
        = rig.template get_bad_p1<T_scalar_param1>();
    std::vector<T_scalar_param2> good_p2
        = rig.template get_good_p2<T_scalar_param2>();
    std::vector<T_scalar_param2> bad_p2
        = rig.template get_bad_p2<T_scalar_param2>();
    std::vector<T_scalar_param3> good_p3
        = rig.template get_good_p3<T_scalar_param3>();
    std::vector<T_scalar_param3> bad_p3
        = rig.template get_bad_p3<T_scalar_param3>();

    // Try a few combinations of parameters that should work
    for (int i = 0; i < 5; i++) {
      internal::shuffle_container(good_p1);
      internal::shuffle_container(good_p2);
      internal::shuffle_container(good_p3);
      assign_parameter_values(p1, good_p1);
      assign_parameter_values(p2, good_p2);
      assign_parameter_values(p3, good_p3);
      EXPECT_NO_THROW(rig.generate_samples(p1, p2, p3, rng));
    }

    // Now try putting incompatible values in first parameter
    if (p1_is_used) {
      for (auto bad_p1_value : bad_p1) {
        assign_parameter_values(p1, {bad_p1_value});
        assign_parameter_values(p2, good_p2);
        assign_parameter_values(p3, good_p3);
        EXPECT_THROW(rig.generate_samples(p1, p2, p3, rng), std::domain_error);
      }
    }

    // Now try putting incompatible values in second parameter
    if (p2_is_used) {
      for (auto bad_p2_value : bad_p2) {
        assign_parameter_values(p1, good_p1);
        assign_parameter_values(p2, {bad_p2_value});
        assign_parameter_values(p3, good_p3);
        EXPECT_THROW(rig.generate_samples(p1, p2, p3, rng), std::domain_error);
      }
    }

    // Now try putting incompatible values in third parameter
    if (p3_is_used) {
      for (auto bad_p3_value : bad_p3) {
        assign_parameter_values(p1, good_p1);
        assign_parameter_values(p2, good_p2);
        assign_parameter_values(p3, {bad_p3_value});
        EXPECT_THROW(rig.generate_samples(p1, p2, p3, rng), std::domain_error);
      }
    }

    // Check to make sure _rng rejects vector arguments of different lengths for
    // all parameter pairs.

    // If p1 is a scalar or the only vector, this test is skipped
    resize_if_vector(p1, 3);  // No-op if p1 is a scalar
    resize_if_vector(p2, 4);  // No-op if p2 is a scalar
    resize_if_vector(p3, 4);  // No-op if p3 is a scalar
    if (stan::math::size(p1) != 1
        && ((p2_is_used && stan::math::size(p2) != 1)
            || (p3_is_used && stan::math::size(p3) != 1))) {
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
    if (p2_is_used && stan::math::size(p2) != 1
        && (stan::math::size(p1) != 1
            || (p3_is_used && stan::math::size(p3) != 1))) {
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
    if (p3_is_used && stan::math::size(p3) != 1
        && (stan::math::size(p1) != 1
            || (p2_is_used && stan::math::size(p2) != 1))) {
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
 * and all combinations of int, std::vector<int>, double, std::vector<double>,
 * Eigen::VectorXd, and Eigen::RowVectorXd as template arguments
 *
 * @tparam T_rig Test rig type for random number generator
 * @param T_rig Test rig for random number generator
 */
template <typename T_rig>
void check_dist_throws_all_types(const T_rig& rig) {
  apply_template_permutations<ArgumentTypes, ArgumentTypes, ArgumentTypes>(
      check_dist_throws{}, rig);
}

/*
 * This function calls check_dist_throws with the given test_rig
 * where the first argument can only be an int or a std::vector<int> and
 * the other arguments can be any combination of int, std::vector<int>,
 * double, std::vector<double>, Eigen::VectorXd, or Eigen::RowVectorXd
 *
 * @tparam T_rig Test rig type for random number generator
 * @param T_rig Test rig for random number generator
 */
template <typename T_rig>
void check_dist_throws_int_first_argument(const T_rig& rig) {
  apply_template_permutations<std::tuple<int, std::vector<int>>, ArgumentTypes,
                              ArgumentTypes>(check_dist_throws{}, rig);
}

/*
 * This function calls check_dist_throws with the given test_rig where
 * the first argument can only be a double or a std::vector<double>
 * and the other arguments can be any combination of int,
 * std::vector<int>, double, std::vector<double>, Eigen::VectorXd, or
 * Eigen::RowVectorXd
 *
 * @tparam T_rig Test rig type for random number generator
 * @param T_rig Test rig for random number generator
 */
template <typename T_rig>
void check_dist_throws_real_first_argument(const T_rig& rig) {
  apply_template_permutations<std::tuple<double, std::vector<double>>,
                              ArgumentTypes, ArgumentTypes>(check_dist_throws{},
                                                            rig);
}

/*
 * Convert a scalar to a length one vector
 *
 * @tparam T Scalar type
 * @param v Input scalar
 * @return vector of length 1 with value v
 */
template <typename T>
std::vector<T> promote_to_vector(T v) {
  return std::vector<T>(1, v);
}

/*
 * For arguments that are already vectors, just copy. Probably more efficient
 * just using std::move but cpplint complained about use of unapproved Rvalues.
 */
template <typename T>
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
  template <typename T_param1, typename T_param2, typename T_param3,
            typename T_rig>
  void operator()(const T_rig& rig) const {
    boost::random::mt19937 rng;
    T_param1 p1;
    T_param2 p2;
    T_param3 p3;
    resize_if_vector(p1, rig.M_);  // No-op if p1 is scalar
    resize_if_vector(p2, rig.M_);  // No-op if p2 is scalar
    resize_if_vector(p3, rig.M_);  // No-op if p3 is scalar

    assign_parameter_values(
        p1,
        rig.template get_good_p1<typename stan::scalar_type<T_param1>::type>());
    assign_parameter_values(
        p2,
        rig.template get_good_p2<typename stan::scalar_type<T_param2>::type>());
    assign_parameter_values(
        p3,
        rig.template get_good_p3<typename stan::scalar_type<T_param3>::type>());

    bool p1_is_used = rig.p1_is_used();
    bool p2_is_used = rig.p2_is_used();
    bool p3_is_used = rig.p3_is_used();

    int M = std::max({(p1_is_used) ? stan::math::size(p1) : 1,
                      (p2_is_used) ? stan::math::size(p2) : 1,
                      (p3_is_used) ? stan::math::size(p3) : 1});

    stan::scalar_seq_view<T_param1> p1_vec(p1);
    stan::scalar_seq_view<T_param2> p2_vec(p2);
    stan::scalar_seq_view<T_param3> p3_vec(p3);

    std::vector<std::vector<double>> samples_to_test_transpose;
    for (int n = 0; n < rig.N_; ++n) {
      // If p1, p2, and p3 are scalars, the output is a scalar. Need to promote
      // it to a std::vector
      samples_to_test_transpose.push_back(
          promote_to_vector(rig.generate_samples(p1, p2, p3, rng)));
    }

    for (int m = 0; m < M; ++m) {
      std::vector<double> samples_to_test;
      for (int n = 0; n < rig.N_; ++n) {
        samples_to_test.push_back(samples_to_test_transpose[n][m]);
      }
      std::vector<double> quantiles
          = rig.generate_quantiles(p1_vec[m], p2_vec[m], p3_vec[m]);

      assert_matches_quantiles(samples_to_test, quantiles, 1e-6);
    }
  }
};

/*
 * Call check_quantiles on the given test rig for a continuous distribution that
 * has no parameters.
 *
 * @tparam T_rig Type of test rig for random number generator
 * @param T_rig Test rig for random number generator
 */
template <typename T_rig>
void check_quantiles_no_params(const T_rig& rig) {
  apply_template_permutations<std::tuple<double>, std::tuple<double>,
                              std::tuple<double>>(check_quantiles{}, rig);
}

/*
 * Call check_quantiles on the given test rig for a continuous distribution that
 * has one parameter that can be an int, std::vector<int>,
 * double, std::vector<double>, Eigen::VectorXd, or an Eigen::RowVectorXd
 *
 * @tparam T_rig Type of test rig for random number generator
 * @param T_rig Test rig for random number generator
 */
template <typename T_rig>
void check_quantiles_real(const T_rig& rig) {
  apply_template_permutations<ArgumentTypes, std::tuple<double>,
                              std::tuple<double>>(check_quantiles{}, rig);
}

/*
 * Call check_quantiles on the given test rig for a continuous distribution that
 * has two parameters that can each be an int, std::vector<int>,
 * double, std::vector<double>, Eigen::VectorXd, or an Eigen::RowVectorXd
 *
 * @tparam T_rig Type of test rig for random number generator
 * @param T_rig Test rig for random number generator
 */
template <typename T_rig>
void check_quantiles_real_real(const T_rig& rig) {
  apply_template_permutations<ArgumentTypes, ArgumentTypes, std::tuple<double>>(
      check_quantiles{}, rig);
}

/*
 * Call check_quantiles on the given test rig for a continuous
 * distribution that has two parameters where the first can only be a
 * double or a std::vector<double> and the other argument can be int,
 * std::vector<int>, double, std::vector<double>, Eigen::VectorXd, or
 * an Eigen::RowVectorXd
 *
 * @tparam T_rig Type of test rig for random number generator
 * @param T_rig Test rig for random number generator
 */
template <typename T_rig>
void check_quantiles_real_first_argument(const T_rig& rig) {
  apply_template_permutations<std::tuple<double, std::vector<double>>,
                              ArgumentTypes, std::tuple<double>>(
      check_quantiles{}, rig);
}

/*
 * Call check_quantiles on the given test rig for a continuous distribution that
 * has three parameters that can each be an int, std::vector<int>,
 * double, std::vector<double>, Eigen::VectorXd, or an Eigen::RowVectorXd
 *
 * @tparam T_rig Type of test rig for random number generator
 * @param T_rig Test rig for random number generator
 */
template <typename T_rig>
void check_quantiles_real_real_real(const T_rig& rig) {
  apply_template_permutations<ArgumentTypes, ArgumentTypes, ArgumentTypes>(
      check_quantiles{}, rig);
}

/*
 * check_counts uses assert_chi_squared to check distributions of discrete
 * random numbers generated by rig.generate_samples against the expected counts
 * computed with with rig.pmf.
 *
 * It isn't possible to do this when the random number has infinite support,
 * so the check is actually performed on a transformed random variable f(x):
 *   f(x) = { x ,              if x is in rig.test_points_
 *            not_test_point,  if x not in rig.test_points_ }
 *
 * The transformed random variable has the support
 *   { rig.test_points_[0], rig.test_points_[1], ..., rig.test_points_[M - 1],
 *     not_test_point }
 *
 *   pmf of f = { pmf(f),        if f(x) is in the support of the tested
 *                               distribution
 *                0.0,           if f(x) is not in the support of the tested
 *                               distribution (this should be enforced by the
 *                               implementation of rig.pmf)
 *                the remainder, if f(x) == not_test_point }
 *
 * After generating for each set of parameters rig.N_ random numbers, the exact
 * and empirical pdfs with assert_chi_squared.
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
struct check_counts {
  template <typename T_param1, typename T_param2, typename T_param3,
            typename T_rig>
  void operator()(const T_rig& rig) const {
    boost::random::mt19937 rng;
    T_param1 p1;
    T_param2 p2;
    T_param3 p3;
    resize_if_vector(p1, rig.M_);  // No-op if p1 is scalar
    resize_if_vector(p2, rig.M_);  // No-op if p2 is scalar
    resize_if_vector(p3, rig.M_);  // No-op if p3 is scalar

    assign_parameter_values(
        p1,
        rig.template get_good_p1<typename stan::scalar_type<T_param1>::type>());
    assign_parameter_values(
        p2,
        rig.template get_good_p2<typename stan::scalar_type<T_param2>::type>());
    assign_parameter_values(
        p3,
        rig.template get_good_p3<typename stan::scalar_type<T_param3>::type>());

    bool p2_is_used = rig.p2_is_used();
    bool p3_is_used = rig.p3_is_used();

    int M = std::max({stan::math::size(p1),
                      (p2_is_used) ? stan::math::size(p2) : 1,
                      (p3_is_used) ? stan::math::size(p3) : 1});

    stan::scalar_seq_view<T_param1> p1_vec(p1);
    stan::scalar_seq_view<T_param2> p2_vec(p2);
    stan::scalar_seq_view<T_param3> p3_vec(p3);

    std::vector<std::vector<int>> samples_to_test_transpose;
    for (int n = 0; n < rig.N_; ++n) {
      // If p1, p2, and p3 are scalars, the output is a scalar. Need to promote
      // it to a std::vector
      samples_to_test_transpose.push_back(
          promote_to_vector(rig.generate_samples(p1, p2, p3, rng)));
    }

    for (int m = 0; m < M; ++m) {
      std::vector<int> samples_to_test;
      for (int n = 0; n < rig.N_; ++n) {
        samples_to_test.push_back(samples_to_test_transpose[n][m]);
      }

      // Generated expected number of counts of transformed random variable
      std::vector<double> epmf;
      double total = 0.0;
      for (size_t n = 0; n < rig.test_points_.size(); ++n) {
        double e
            = rig.N_
              * rig.pmf(rig.test_points_[n], p1_vec[m], p2_vec[m], p3_vec[m]);
        epmf.push_back(e);
        total += e;
      }
      epmf.push_back(rig.N_ - total);

      // Generate samples of transformed random variable
      std::map<int, int> count_map;
      int remainder = 0;
      for (size_t n = 0; n < rig.test_points_.size(); ++n)
        count_map[rig.test_points_[n]] = 0;
      for (size_t n = 0; n < samples_to_test.size(); ++n) {
        int sample = samples_to_test[n];
        if (count_map.find(sample) != count_map.end()) {
          count_map[sample] += 1;
        } else {
          remainder += 1;
        }
      }

      // Transform to vector
      std::vector<int> counts;
      for (size_t n = 0; n < rig.test_points_.size(); ++n) {
        counts.push_back(count_map[rig.test_points_[n]]);
      }
      counts.push_back(remainder);

      // Trim the zero probability outputs
      std::vector<int> counts_trimmed;
      std::vector<double> epmf_trimmed;
      for (size_t n = 0; n < counts.size(); ++n) {
        if (epmf[n] > 0.0) {
          counts_trimmed.push_back(counts[n]);
          epmf_trimmed.push_back(epmf[n]);
        }
      }

      // If there's only one non-zero probability output, just sanity check it
      if (counts_trimmed.size() == 1) {
        EXPECT_EQ(static_cast<double>(counts_trimmed[0]), epmf_trimmed[0]);
      } else {
        assert_chi_squared(counts_trimmed, epmf_trimmed, 1e-6);
      }
    }
  }
};

/*
 * Call check_counts on the given test rig for a pmf that takes one
 * argument that can be an int, std::vector<int>, double,
 * std::vector<double>, Eigen::VectorXd, or an Eigen::RowVectorXd
 *
 * @tparam T_rig Type of test rig for random number generator
 * @param T_rig Test rig for random number generator
 */
template <typename T_rig>
void check_counts_real(const T_rig& rig) {
  apply_template_permutations<ArgumentTypes, std::tuple<double>,
                              std::tuple<double>>(check_counts{}, rig);
}

/*
 * Call check_counts on the given test rig for a pmf that takes two
 * arguments that can each be an int, std::vector<int>, double,
 * std::vector<double>, Eigen::VectorXd, or an Eigen::RowVectorXd
 *
 * @tparam T_rig Type of test rig for random number generator
 * @param T_rig Test rig for random number generator
 */
template <typename T_rig>
void check_counts_real_real(const T_rig& rig) {
  apply_template_permutations<ArgumentTypes, ArgumentTypes, std::tuple<double>>(
      check_counts{}, rig);
}

/*
 * Call check_counts on the given test rig for a pmf that takes three
 * arguments that can each be an int, std::vector<int>, double,
 * std::vector<double>, Eigen::VectorXd, or an Eigen::RowVectorXd
 *
 * @tparam T_rig Type of test rig for random number generator
 * @param T_rig Test rig for random number generator
 */
template <typename T_rig>
void check_counts_real_real_real(const T_rig& rig) {
  apply_template_permutations<ArgumentTypes, ArgumentTypes, ArgumentTypes>(
      check_counts{}, rig);
}

/*
 * Call check_counts on the given test rig for a pmf that takes for the first
 * argument either an int or a std::vector<int> and then for its second
 * argument an int, std::vector<int>, double, std::vector<double>,
 * Eigen::VectorXd, or an Eigen::RowVectorXd
 *
 * @tparam T_rig Type of test rig for random number generator
 * @param T_rig Test rig for random number generator
 */
template <typename T_rig>
void check_counts_int_real(const T_rig& rig) {
  apply_template_permutations<std::tuple<int, std::vector<int>>, ArgumentTypes,
                              std::tuple<double>>(check_counts{}, rig);
}

/*
 * Call check_counts on the given test rig for a pmf that takes for the first
 * argument either an int or a std::vector<int> and then for each of its second
 * or third arguments an int, std::vector<int>, double, std::vector<double>,
 * Eigen::VectorXd, or an Eigen::RowVectorXd
 *
 * @tparam T_rig Type of test rig for random number generator
 * @param T_rig Test rig for random number generator
 */
template <typename T_rig>
void check_counts_int_real_real(const T_rig& rig) {
  apply_template_permutations<std::tuple<int, std::vector<int>>, ArgumentTypes,
                              ArgumentTypes>(check_counts{}, rig);
}

#endif
