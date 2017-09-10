#ifndef TEST_UNIT_MATH_PRIM_MAT_PROB_VECTOR_RNG_TEST_HELPER_HPP
#define TEST_UNIT_MATH_PRIM_MAT_PROB_VECTOR_RNG_TEST_HELPER_HPP

#include <gtest/gtest.h>
#include <boost/math/distributions.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <stan/math/prim/mat.hpp>
#include <test/unit/math/prim/scal/prob/util.hpp>
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
  if(values.size() == 0)
    return;

  // Static cast the size here because sometimes it is unsigned (if T_param
  // is a std::vector) and sometimes it is signed (when it is an Eigen type)
  for (size_t i = 0; i < static_cast<size_t>(params.size()); i++) {
    params[i] = values[i % values.size()];
  }
}

/*
 * Assign param the first value of values
 *
 * @param param Output parameter to write value to
 * @param params Vector with value to copy into param
 */
template <>
void assign_parameter_values(double& param, const std::vector<double>& values) {
  if(values.size() == 0)
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
 * check_dist_throws feeds the given generate_samples callable various
 * combinations of valid and invalid parameters (as defined by good_p1_, bad_p1,
 * good_p2_, bad_p2_, good_p3_, and bad_p3_). For all combinations of valid
 * (good) parameters, generate_samples should throw no errors. For all
 * combinations with an invalid (bad) parameter, generate_samples should throw
 * domain_errors.
 *
 * If good_p2_ or good_p3_ are empty, then it is assumed that those parameters
 * are unused and will not be tested.
 *
 * The generate_samples callable is also passed various other guaranteed invalid
 * values like positive infinity, negative infinity, and NaNs (these should also
 * cause domain_errors).
 *
 * The generated_samples callable is also tested to reject incompatibly sized
 * input vectors. These should cause invalid_argument errors.
 *
 * The generate_samples callable must have the signature:
 *   T_out generate_samples(T_param1 p1, T_param2 p2, T_param3, T_rng)
 * where T_param1, T_param2, and T_param3 can be any combination of:
 *   double, std::vector<double>, VectorXd, or RowVectorXd
 * and T_rng is a boost random number generator. p2 and p3 are required for the
 * function signature, but they can be ignored if they are not needed.
 *
 * The output of generate_samples must be of
 *   stan::length(out) == stan::length(p1) if only p1 is used
 *   stan::length(out) == stan::max_size(p1, p2) if p1 and p2 are used
 *   stan::length(out) == stan::max_size(p1, p2, p3) if all parameters are used
 *
 * Each element of the output, out[i], should correspond to a random variate
 * generated according to parameters p1[i], p2[i], and p3[i]
 *
 * @tparam T_param1 Type of first parameter
 * @tparam T_param2 Type of second parameter
 * @tparam T_param3 Type of third parameter
 * @tparam T_generate_samples Type of generate_samples functor
 * @param generate_samples generate_samples functor to test
 * @param good_p1_ Valid values of first parameter
 * @param bad_p1_ Invalid values of first parameter
 * @param good_p2_ Valid values of second parameter
 * (leave empty if p2 is unused)
 * @param bad_p2_ Invalid values of second parameter
 * @param good_p3_ Valid values of third parameter
 * (leave empty if p3 is unused)
 * @param bad_p3_ Invalid values of third parameter
 */
template<typename T_param1, typename T_param2, typename T_param3,
         typename T_generate_samples>
void check_dist_throws(T_generate_samples generate_samples,
                       const std::vector<double>& good_p1_,
                       const std::vector<double>& bad_p1_,
                       const std::vector<double>& good_p2_,
                       const std::vector<double>& bad_p2_,
                       const std::vector<double>& good_p3_,
                       const std::vector<double>& bad_p3_) {
  boost::random::mt19937 rng;

  T_param1 p1;
  T_param2 p2;
  T_param2 p3;

  bool p2_is_used = good_p2_.size() > 0;
  bool p3_is_used = good_p3_.size() > 0;

  resize_if_vector(p1, 5);  // No-op if p1 is a scalar
  resize_if_vector(p2, 5);  // No-op if p2 is a scalar
  resize_if_vector(p3, 5);  // No-op if p3 is a scalar

  // Make copies of the input arguments so that we can randomly shuffle them
  // in the tests
  std::vector<double> good_p1 = good_p1_;
  std::vector<double> bad_p1 = bad_p1_;
  std::vector<double> good_p2 = good_p2_;
  std::vector<double> bad_p2 = bad_p2_;
  std::vector<double> good_p3 = good_p3_;
  std::vector<double> bad_p3 = bad_p3_;

  // Try a few combinations of parameters that should work
  for (int i = 0; i < 5; i++) {
    std::random_shuffle(good_p1.begin(), good_p1.end());
    std::random_shuffle(good_p2.begin(), good_p2.end());
    std::random_shuffle(good_p3.begin(), good_p3.end());
    assign_parameter_values(p1, good_p1);
    assign_parameter_values(p2, good_p2);
    assign_parameter_values(p3, good_p3);
    EXPECT_NO_THROW(generate_samples(p1, p2, p3, rng));
  }

  // Now try putting incompatible values in first parameter
  for (auto bad_p1_value : bad_p1) {
    assign_parameter_values(p1, { bad_p1_value });
    assign_parameter_values(p2, good_p2);
    assign_parameter_values(p3, good_p3);
    EXPECT_THROW(generate_samples(p1, p2, p3, rng), std::domain_error);
  }

  // Now try putting incompatible values in second parameter
  for (auto bad_p2_value : bad_p2) {
    assign_parameter_values(p1, good_p1);
    assign_parameter_values(p2, { bad_p2_value });
    assign_parameter_values(p3, good_p3);
    EXPECT_THROW(generate_samples(p1, p2, p3, rng), std::domain_error);
  }

  // Now try putting incompatible values in third parameter
  for (auto bad_p3_value : bad_p3) {
    assign_parameter_values(p1, good_p1);
    assign_parameter_values(p2, good_p2);
    assign_parameter_values(p3, { bad_p3_value });
    EXPECT_THROW(generate_samples(p1, p2, p3, rng), std::domain_error);
  }

  // Make sure sampler throws errors with these values
  std::vector<double> bad_values = { stan::math::positive_infinity(),
                                     stan::math::negative_infinity(),
                                     stan::math::not_a_number() };

  for (auto bad_value : bad_values) {
    assign_parameter_values(p1, { bad_value });
    assign_parameter_values(p2, good_p2);
    assign_parameter_values(p3, good_p3);
    EXPECT_THROW(generate_samples(p1, p2, p3, rng), std::domain_error);

    if(p2_is_used) {
      assign_parameter_values(p1, good_p1);
      assign_parameter_values(p2, { bad_value });
      EXPECT_THROW(generate_samples(p1, p2, p3, rng), std::domain_error);
    }

    if(p3_is_used) {
      assign_parameter_values(p2, good_p2);
      assign_parameter_values(p3, { bad_value });
      EXPECT_THROW(generate_samples(p1, p2, p3, rng), std::domain_error);
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
    EXPECT_THROW(generate_samples(p1, p2, p3, rng), std::invalid_argument);
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
    EXPECT_THROW(generate_samples(p1, p2, p3, rng), std::invalid_argument);
  }

  // If p3 is a scalar or the only vector, this test is skipped
  resize_if_vector(p1, 4);  // No-op if p1 is a scalar
  resize_if_vector(p2, 4);  // No-op if p2 is a scalar
  resize_if_vector(p3, 3);  // No-op if p3 is a scalar
  if (p3_is_used && stan::length(p3) != 1 &&
      (stan::length(p1) != 1 || (p2_is_used && stan::length(p2)))) {
    assign_parameter_values(p1, good_p1);
    assign_parameter_values(p2, good_p2);
    assign_parameter_values(p3, good_p3);
    EXPECT_THROW(generate_samples(p1, p2, p3, rng), std::invalid_argument);
  }
}

/*
 * This function calls check_dist_throws with the given generate_samples
 * callable and all combinations of double, std::vector<double>,
 * Eigen::VectorXd, and Eigen::RowVectorXd as template arguments
 *
 * @tparam T_generate_samples Type of generate_samples functor
 * @param generate_samples generate_samples functor to test
 * @param good_p1 Valid values of first parameter
 * @param bad_p1 Invalid values of first parameter
 * @param good_p2 Valid values of second parameter
 * @param bad_p2 Invalid values of second parameter
 * @param good_p3 Valid values of third parameter
 * @param bad_p3 Invalid values of third parameter
 */
template<typename T_generate_samples>
void check_dist_throws_all_types(T_generate_samples generate_samples,
                                 const std::vector<double>& good_p1,
                                 const std::vector<double>& bad_p1,
                                 const std::vector<double>& good_p2,
                                 const std::vector<double>& bad_p2,
                                 const std::vector<double>& good_p3,
                                 const std::vector<double>& bad_p3) {
  using Eigen::VectorXd;
  using Eigen::RowVectorXd;

  check_dist_throws<double, double, double>
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<double, std::vector<double>, double >
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<double, VectorXd, double>
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<double, RowVectorXd, double>
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<std::vector<double>, double, double>
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<std::vector<double>, std::vector<double>, double >
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<std::vector<double>, VectorXd, double>
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<std::vector<double>, RowVectorXd, double>
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<VectorXd, double, double>
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<VectorXd, std::vector<double>, double >
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<VectorXd, VectorXd, double>
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<VectorXd, RowVectorXd, double>
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<RowVectorXd, double, double>
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<RowVectorXd, std::vector<double>, double >
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<RowVectorXd, VectorXd, double>
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<RowVectorXd, RowVectorXd, double>
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);

  check_dist_throws<double, double, std::vector<double> >
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<double, std::vector<double>, std::vector<double> >
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<double, VectorXd, std::vector<double> >
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<double, RowVectorXd, std::vector<double> >
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<std::vector<double>, double, std::vector<double> >
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<std::vector<double>, std::vector<double>, std::vector<double> >
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<std::vector<double>, VectorXd, std::vector<double> >
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<std::vector<double>, RowVectorXd, std::vector<double> >
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<VectorXd, double, std::vector<double> >
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<VectorXd, std::vector<double>, std::vector<double> >
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<VectorXd, VectorXd, std::vector<double> >
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<VectorXd, RowVectorXd, std::vector<double> >
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<RowVectorXd, double, std::vector<double> >
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<RowVectorXd, std::vector<double>, std::vector<double> >
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<RowVectorXd, VectorXd, std::vector<double> >
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<RowVectorXd, RowVectorXd, std::vector<double> >
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);

  check_dist_throws<double, double, VectorXd>
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<double, std::vector<double>, VectorXd>
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<double, VectorXd, VectorXd>
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<double, RowVectorXd, VectorXd>
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<std::vector<double>, double, VectorXd>
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<std::vector<double>, std::vector<double>, VectorXd>
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<std::vector<double>, VectorXd, VectorXd>
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<std::vector<double>, RowVectorXd, VectorXd>
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<VectorXd, double, VectorXd>
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<VectorXd, std::vector<double>, VectorXd>
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<VectorXd, VectorXd, VectorXd>
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<VectorXd, RowVectorXd, VectorXd>
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<RowVectorXd, double, VectorXd>
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<RowVectorXd, std::vector<double>, VectorXd>
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<RowVectorXd, VectorXd, VectorXd>
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<RowVectorXd, RowVectorXd, VectorXd>
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);

  check_dist_throws<double, double, RowVectorXd>
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<double, std::vector<double>, RowVectorXd>
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<double, VectorXd, RowVectorXd>
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<double, RowVectorXd, RowVectorXd>
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<std::vector<double>, double, RowVectorXd>
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<std::vector<double>, std::vector<double>, RowVectorXd>
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<std::vector<double>, VectorXd, RowVectorXd>
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<std::vector<double>, RowVectorXd, RowVectorXd>
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<VectorXd, double, RowVectorXd>
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<VectorXd, std::vector<double>, RowVectorXd>
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<VectorXd, VectorXd, RowVectorXd>
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<VectorXd, RowVectorXd, RowVectorXd>
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<RowVectorXd, double, RowVectorXd>
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<RowVectorXd, std::vector<double>, RowVectorXd>
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<RowVectorXd, VectorXd, RowVectorXd>
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
  check_dist_throws<RowVectorXd, RowVectorXd, RowVectorXd>
    (generate_samples, good_p1, bad_p1, good_p2, bad_p2, good_p3, bad_p3);
}

/*
 * Convert a scalar to a length one vector
 *
 * @param v Input scalar
 * @return vector of length 1 with value v
 */
std::vector<double> promote_to_vector(double v) {
  return std::vector<double>(1, v);
}

/*
 * For arguments that are already vectors, just copy. Probably more efficient
 * just using std::move but cpplint complained about use of unapproved Rvalues.
 */
std::vector<double> promote_to_vector(std::vector<double> v) {
  return v;
}

/*
 * check_quantiles uses assert_matches_quantiles to check random numbers
 * generated by the callable generate_samples against the reference
 * quantiles computed with the given callable generate_quantiles.
 *
 * If good_p2 or good_p3 are empty, then it is assumed that those parameters
 * are unused and will not be tested.
 *
 * The generate_samples callable must have the signature:
 *   T_out generate_samples(T_param1 p1, T_param2 p2, T_param3, p3, T_rng)
 * where T_param1, T_param2, and T_param3 can be any combination of:
 *   double, std::vector<double>, VectorXd, or RowVectorXd
 * and T_rng is a boost random number generator.  p2 and p3 are required for the
 * function signature, but they can be ignored if they are not needed.
 *
 * The output of generate_samples must be of
 *   stan::length(out) == stan::length(p1) if only p1 is used
 *   stan::length(out) == stan::max_size(p1, p2) if p1 and p2 are used
 *   stan::length(out) == stan::max_size(p1, p2, p3) if all parameters are used
 *
 * Each element of the output, out[i], should correspond to a random variate
 * generated according to parameters p1[i], p2[i], and p3[i]
 *
 * The generate_quantiles callable must have the signature:
 *  std::vector<double> generate_samples(int N, double p1, double p2, double p3)
 * (similarly to above, ignore parameters p2 and p3 if they aren't necessary)
 *
 * N is the number of samples that the quantiles will be compared
 * against and p1, p2, and p3 are the values of the parameters that these
 * quantiles correspond to.
 *
 * The output of generate_quantiles should be the quantiles used for comparing
 * empirical distributions used by assert_matches_quantiles
 *
 * Parameters for each of these callables are copied from the vectors good_p1,
 * good_p2, and good_p3.
 *
 * @tparam T_param1 Type of first parameter
 * @tparam T_param2 Type of second parameter
 * @tparam T_param3 Type of third parameter
 * @tparam T_generate_samples Type of callable for generating samples
 * @tparam T_generate_quantiles Type of callable for generating quantiles
 * @param sampler sampler functor to test
 * @param good_p1 Valid values of first parameter
 * @param good_p2 Valid values of second parameter
 * @param good_p3 Valid values of third parameter
 */
template<typename T_param1, typename T_param2, typename T_param3,
         typename T_generate_samples, typename T_generate_quantiles>
void check_quantiles(int N, int M_vec,
                     T_generate_samples generate_samples,
                     T_generate_quantiles generate_quantiles,
                     const std::vector<double>& good_p1,
                     const std::vector<double>& good_p2,
                     const std::vector<double>& good_p3) {
  boost::random::mt19937 rng;
  T_param1 p1;
  T_param2 p2;
  T_param3 p3;
  resize_if_vector(p1, M_vec);  // No-op if p1 is scalar
  resize_if_vector(p2, M_vec);  // No-op if p2 is scalar
  resize_if_vector(p3, M_vec);  // No-op if p3 is scalar

  assign_parameter_values(p1, good_p1);
  assign_parameter_values(p2, good_p2);
  assign_parameter_values(p3, good_p3);

  bool p2_is_used = good_p2.size() > 0;
  bool p3_is_used = good_p3.size() > 0;

  int M = std::max({ stan::length(p1),
                     (p2_is_used) ? stan::length(p2) : 1,
                     (p3_is_used) ? stan::length(p3) : 1 });

  stan::scalar_seq_view<T_param1> p1_vec(p1);
  stan::scalar_seq_view<T_param2> p2_vec(p2);
  stan::scalar_seq_view<T_param3> p3_vec(p3);

  std::vector<std::vector<double> > samples_to_test_transpose;
  for (int n = 0; n < N; ++n) {
    // If p1, p2, and p3 are scalars, the output is a scalar. Need to promote it
    // to a std::vector
    samples_to_test_transpose.
      push_back(promote_to_vector(generate_samples(p1, p2, p3, rng)));
  }

  for (int m = 0; m < M; ++m) {
    std::vector<double> samples_to_test;
    for (int n = 0; n < N; ++n) {
      samples_to_test.push_back(samples_to_test_transpose[n][m]);
    }
    std::vector<double> quantiles =
      generate_quantiles(N, p1_vec[m], p2_vec[m], p3_vec[m]);

    assert_matches_quantiles(samples_to_test, quantiles, 1e-6);
  }
}

/*
 * Call check_quantiles with the given arguments and all combinations of
 * double, std::vector<double>, Eigen::VectorXd, and Eigen::RowVectorXd
 *
 * @tparam T_generate_samples Type of callable for generating samples
 * @tparam T_generate_quantiles Type of callable for generating quantiles
 * @param sampler sampler functor to test
 * @param good_p1 Valid values of first parameter
 * @param good_p2 Valid values of second parameter
 */
template<typename T_generate_samples, typename T_generate_quantiles>
void check_quantiles_all_types(int N, int M,
                               T_generate_samples generate_samples,
                               T_generate_quantiles generate_quantiles,
                               const std::vector<double>& good_p1,
                               const std::vector<double>& good_p2,
                               const std::vector<double>& good_p3) {
  using Eigen::VectorXd;
  using Eigen::RowVectorXd;

  check_quantiles<double, double, double>
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<double, std::vector<double>, double>
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<double, VectorXd, double>
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<double, RowVectorXd, double>
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<std::vector<double>, double, double>
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<std::vector<double>, std::vector<double>, double>
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<std::vector<double>, VectorXd, double>
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<std::vector<double>, RowVectorXd, double>
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<VectorXd, double, double>
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<VectorXd, std::vector<double>, double>
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<VectorXd, VectorXd, double>
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<VectorXd, RowVectorXd, double>
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<RowVectorXd, double, double>
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<RowVectorXd, std::vector<double>, double>
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<RowVectorXd, VectorXd, double>
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<RowVectorXd, RowVectorXd, double>
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);

  check_quantiles<double, double, std::vector<double> >
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<double, std::vector<double>, std::vector<double> >
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<double, VectorXd, std::vector<double> >
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<double, RowVectorXd, std::vector<double> >
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<std::vector<double>, double, std::vector<double> >
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<std::vector<double>, std::vector<double>, std::vector<double> >
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<std::vector<double>, VectorXd, std::vector<double> >
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<std::vector<double>, RowVectorXd, std::vector<double> >
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<VectorXd, double, std::vector<double> >
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<VectorXd, std::vector<double>, std::vector<double> >
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<VectorXd, VectorXd, std::vector<double> >
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<VectorXd, RowVectorXd, std::vector<double> >
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<RowVectorXd, double, std::vector<double> >
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<RowVectorXd, std::vector<double>, std::vector<double> >
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<RowVectorXd, VectorXd, std::vector<double> >
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<RowVectorXd, RowVectorXd, std::vector<double> >
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);

  check_quantiles<double, double, VectorXd>
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<double, std::vector<double>, VectorXd>
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<double, VectorXd, VectorXd>
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<double, RowVectorXd, VectorXd>
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<std::vector<double>, double, VectorXd>
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<std::vector<double>, std::vector<double>, VectorXd>
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<std::vector<double>, VectorXd, VectorXd>
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<std::vector<double>, RowVectorXd, VectorXd>
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<VectorXd, double, VectorXd>
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<VectorXd, std::vector<double>, VectorXd>
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<VectorXd, VectorXd, VectorXd>
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<VectorXd, RowVectorXd, VectorXd>
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<RowVectorXd, double, VectorXd>
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<RowVectorXd, std::vector<double>, VectorXd>
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<RowVectorXd, VectorXd, VectorXd>
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<RowVectorXd, RowVectorXd, VectorXd>
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);

  check_quantiles<double, double, RowVectorXd>
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<double, std::vector<double>, RowVectorXd>
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<double, VectorXd, RowVectorXd>
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<double, RowVectorXd, RowVectorXd>
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<std::vector<double>, double, RowVectorXd>
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<std::vector<double>, std::vector<double>, RowVectorXd>
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<std::vector<double>, VectorXd, RowVectorXd>
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<std::vector<double>, RowVectorXd, RowVectorXd>
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<VectorXd, double, RowVectorXd>
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<VectorXd, std::vector<double>, RowVectorXd>
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<VectorXd, VectorXd, RowVectorXd>
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<VectorXd, RowVectorXd, RowVectorXd>
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<RowVectorXd, double, RowVectorXd>
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<RowVectorXd, std::vector<double>, RowVectorXd>
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<RowVectorXd, VectorXd, RowVectorXd>
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
  check_quantiles<RowVectorXd, RowVectorXd, RowVectorXd>
    (N, M, generate_samples, generate_quantiles, good_p1, good_p2, good_p3);
}

#endif
