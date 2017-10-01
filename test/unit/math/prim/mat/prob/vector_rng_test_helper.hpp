#ifndef TEST_UNIT_MATH_PRIM_MAT_PROB_VECTOR_RNG_TEST_HELPER_HPP
#define TEST_UNIT_MATH_PRIM_MAT_PROB_VECTOR_RNG_TEST_HELPER_HPP

#include <gtest/gtest.h>
#include <boost/math/distributions.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <stan/math/prim/mat.hpp>
#include <test/unit/math/prim/scal/prob/util.hpp>
#include <vector>

/*
 * VectorRNGTestRig is an Abstract class for wrapping up random number
 * generators in a way so that they can be tested to
 *  1) provide interfaces to all the necessary types
 *  2) generate the correct distributions
 *
 * A test rig that inherits from VectorRNGTestRig must implement
 *  1) generate_samples
 *  2) generate_quantiles
 * and must initialize good_p1_, bad_p1_, good_p2_, bad_p2, good_p3_, bad_p3_
 * (depending on how many parameters the wrapped random number generator uses)
 *
 * generate_samples is a wrapper around the actual random number generator. It
 * is expected to have the signature:
 *   T_out generate_samples(T_param1 p1, T_param2 p2, T_param3, T_rng)
 * where T_param1, T_param2, and T_param3 can be any combination of:
 *   double, std::vector<double>, VectorXd, or RowVectorXd
 * and T_rng is a boost random number generator. p2 and p3 are required for the
 * function signature, but they can be ignored if they are not needed.
 *
 * Each element of the output, out[i], should correspond to a random variate
 * generated according to parameters p1[i], p2[i], and p3[i]
 *
 * generate_samples is expected to throw domain_errors when passed invalid
 * values like positive infinity, negative infinity, and NaNs.
 *
 * generated_samples is also expected to reject incompatibly sized
 * input vectors. These should cause invalid_argument errors.
 *
 * The output of generate_samples must be of
 *   stan::length(out) == stan::length(p1) if only p1 is used
 *   stan::length(out) == stan::max_size(p1, p2) if p1 and p2 are used
 *   stan::length(out) == stan::max_size(p1, p2, p3) if all parameters are used
 * It must be defined to take three arguments, but not all must be used.
 *
 * The generate_quantiles callable must have the signature:
 *  std::vector<double> generate_samples(double p1, double p2, double p3)
 * (similarly to above, ignore parameters p2 and p3 if they aren't necessary)
 *
 * generate_quantiles should compute the quantiles that will be compared against
 * (using assert_matches_quantiles) empirical distributions generated with
 * generate_samples (which will be of length N_ for each parameter combination).
 *
 * p1, p2, and p3 are the values of the parameters that these quantiles
 * correspond to.
 *
 * good_p1_ and bad_p1_ should be initialized to lists of valid and invalid
 * parameters for the first parameter of the tested RNG
 *
 * Likewise, good_p2_, bad_p2_, and good_p3_, bad_p3_ should be initialized.
 *
 * If testing a distribution that uses only two parameters, leave good_p3_ as
 * a length zero vector and it will not be tested. Similarly, if good_p2_ is
 * length zero, the second parameter will not be tested either
 */
class VectorRNGTestRig {
public:
  int N_; // Number of samples used in the quantiles tests
  int M_; // Length of vectors for the vectorization tests

  std::vector<double> good_p1_;
  std::vector<double> bad_p1_;
  std::vector<double> good_p2_;
  std::vector<double> bad_p2_;
  std::vector<double> good_p3_;
  std::vector<double> bad_p3_;

  /*
   * This function wraps up the random number generator for testing.
   *
   * The tested rng can have up to three parameters. Any unused parameters can
   * be ignored.
   *
   * It *must* be implemented in the child class (it isn't virtual here because
   * C++ doesn't allow templated virtual member functions)
   */
  template<typename T1, typename T2, typename T3, typename T_rng>
  auto generate_samples(const T1& p1, const T2& p2, const T3& p3,
                        T_rng& rng) const;

  /*
   * This function builds the quantiles that we will supply to
   * assert_matches_quantiles to test the random number generator implemented in
   * generate_samples
   */
  virtual std::vector<double> generate_quantiles(double p1, double p2,
                                                 double p3) const = 0;

  VectorRNGTestRig(int N, int M,
                   std::vector<double> good_p1,
                   std::vector<double> bad_p1,
                   std::vector<double> good_p2,
                   std::vector<double> bad_p2,
                   std::vector<double> good_p3,
                   std::vector<double> bad_p3) : N_(N), M_(M),
                                                 good_p1_(good_p1),
                                                 bad_p1_(bad_p1),
                                                 good_p2_(good_p2),
                                                 bad_p2_(bad_p2),
                                                 good_p3_(good_p3),
                                                 bad_p3_(bad_p3) {
  }

  VectorRNGTestRig(int N, int M,
                   std::vector<double> good_p1,
                   std::vector<double> bad_p1,
                   std::vector<double> good_p2,
                   std::vector<double> bad_p2) : N_(N), M_(M),
                                                 good_p1_(good_p1),
                                                 bad_p1_(bad_p1),
                                                 good_p2_(good_p2),
                                                 bad_p2_(bad_p2){
  }
  
  VectorRNGTestRig(int N, int M,
                   std::vector<double> good_p1,
                   std::vector<double> bad_p1) : N_(N), M_(M),
                                                 good_p1_(good_p1),
                                                 bad_p1_(bad_p1) {
  }
};

/*
 * call_permutations_helper is the primary function for the recursive template
 * function call_permutations. It calls func<type[I], type[J], type[K]>(rig),
 * and continues the recursion by decrementing K.
 *
 * tparam T_typelist Tuple templated with list of types to iterate through
 * tparam T_functor Templated functor
 * tparam T_rig Random number generator test rig
 * tparam I Loop index for first type
 * tparam J Loop index for second type
 * tparam K Loop index for third type
 */
template<typename T_typelist, typename T_functor, typename T_rig, int I, int J, int K>
struct call_permutations_helper {
  void operator()(const T_functor& func, const T_rig& rig) const {
    func.template operator()<typename std::tuple_element<I, T_typelist>::type,
                             typename std::tuple_element<J, T_typelist>::type,
                             typename std::tuple_element<K, T_typelist>::type>
      (rig);

    call_permutations_helper<T_typelist, T_functor, T_rig, I, J, K - 1>{}(func,
                                                                          rig);
  }
};

/*
 * This edge case catches when the right-most argument has completed one
 * iteration through all types. It computes func<type[I], type[J], type[0]>,
 * carries from the second argument, and continues the recursion.
 *
 * tparam T_typelist Tuple templated with list of types to iterate through
 * tparam T_functor Templated functor
 * tparam T_rig Random number generator test rig
 * tparam I Loop index for first type
 * tparam J Loop index for second type
 * tparam K Loop index for third type
 */
template<typename T_typelist, typename T_functor, typename T_rig, int I, int J>
struct call_permutations_helper<T_typelist, T_functor, T_rig, I, J, 0> {
  void operator()(const T_functor& func, const T_rig& rig) const {
    func.template operator()<typename std::tuple_element<I, T_typelist>::type,
                             typename std::tuple_element<J, T_typelist>::type,
                             typename std::tuple_element<0, T_typelist>::type>
      (rig);
    
    call_permutations_helper<T_typelist, T_functor, T_rig, I, J - 1,
                             std::tuple_size<T_typelist>::value - 1>{}(func, rig);
  }
};

/*
 * This edge case catches when the second argument has completed one iteration
 * through all types. It computes func<type[I], type[0], type[0]>, carries from
 * the left-most argument and continues the recursion.
 *
 * tparam T_typelist Tuple templated with list of types to iterate through
 * tparam T_functor Templated functor
 * tparam T_rig Random number generator test rig
 * tparam I Loop index for first type
 */
template<typename T_typelist, typename T_functor, typename T_rig, int I>
struct call_permutations_helper<T_typelist, T_functor, T_rig, I, 0, 0> {
  void operator()(const T_functor& func, const T_rig& rig) const {
    func.template operator()<typename std::tuple_element<I, T_typelist>::type,
                             typename std::tuple_element<0, T_typelist>::type,
                             typename std::tuple_element<0, T_typelist>::type>
      (rig);
    
    call_permutations_helper<T_typelist, T_functor, T_rig, I - 1,
                             std::tuple_size<T_typelist>::value - 1,
                             std::tuple_size<T_typelist>::value - 1>{}(func,
                                                                       rig);
  }
};

/*
 * This edge case catches when all types have been iterated through for each
 * template argument. It computes func<type[0], type[0], type[0]> and ends the
 * recursion.
 *
 * tparam T_typelist Tuple templated with list of types to iterate through
 * tparam T_functor Templated functor
 * tparam T_rig Random number generator test rig
 */
template<typename T_typelist, typename T_functor, typename T_rig>
struct call_permutations_helper<T_typelist, T_functor, T_rig, 0, 0, 0> {
  void operator()(const T_functor& func, const T_rig& rig) const {
    func.template operator()<typename std::tuple_element<0, T_typelist>::type,
                             typename std::tuple_element<0, T_typelist>::type,
                             typename std::tuple_element<0, T_typelist>::type>
      (rig);
  }
};


/*
 * call_permutations uses template recursion to call the functor func
 * with all 3-element permutations of the types in T_typelist
 *
 * Roughly, this corresponds to:
 *  for (i in (I - 1):0) // indexes go backwards, like: I - 1, I - 2, ... 0
 *   for (j in (J - 1):0)
 *    for (k in (K - 1):0)
 *     func<type[i], type[j], type[k]>(rig)
 *
 * tparam T_typelist Tuple templated with list of types to iterate through
 * tparam T_functor Templated functor
 * tparam T_rig Random number generator test rig
 */
template<typename T_typelist, typename T_functor, typename T_rig>
void call_permutations(const T_functor& func, const T_rig& rig) {
  call_permutations_helper
    <T_typelist, T_functor, T_rig,
     std::tuple_size<T_typelist>::value - 1,
     std::tuple_size<T_typelist>::value - 1,
     std::tuple_size<T_typelist>::value - 1>{}(func, rig);
}

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
 * Fill the vector-like variable params with values from the values argument.
 *
 * values can be shorter than params (in which cause multiple copies of values
 * are tiled over the params vector)
 *
 * @param params Parameter vector to write values to
 * @param params Values to copy into params
 */
template <>
void assign_parameter_values(std::vector<double>& params,
                             const std::vector<double>& values) {
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
template <>
void assign_parameter_values(double& param, const std::vector<double>& values) {
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
  template<typename T_param1, typename T_param2, typename T_param3, typename T_rig>
  void operator()(const T_rig& rig) const {
    boost::random::mt19937 rng;

    T_param1 p1;
    T_param2 p2;
    T_param3 p3;

    bool p2_is_used = rig.good_p2_.size() > 0;
    bool p3_is_used = rig.good_p3_.size() > 0;

    resize_if_vector(p1, 5);  // No-op if p1 is a scalar
    resize_if_vector(p2, 5);  // No-op if p2 is a scalar
    resize_if_vector(p3, 5);  // No-op if p3 is a scalar

    // Make copies of the input arguments so that we can randomly shuffle them
    // in the tests
    std::vector<double> good_p1 = rig.good_p1_;
    std::vector<double> bad_p1 = rig.bad_p1_;
    std::vector<double> good_p2 = rig.good_p2_;
    std::vector<double> bad_p2 = rig.bad_p2_;
    std::vector<double> good_p3 = rig.good_p3_;
    std::vector<double> bad_p3 = rig.bad_p3_;

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
    for (auto bad_p2_value : bad_p2) {
      assign_parameter_values(p1, good_p1);
      assign_parameter_values(p2, { bad_p2_value });
      assign_parameter_values(p3, good_p3);
      EXPECT_THROW(rig.generate_samples(p1, p2, p3, rng), std::domain_error);
    }

    // Now try putting incompatible values in third parameter
    for (auto bad_p3_value : bad_p3) {
      assign_parameter_values(p1, good_p1);
      assign_parameter_values(p2, good_p2);
      assign_parameter_values(p3, { bad_p3_value });
      EXPECT_THROW(rig.generate_samples(p1, p2, p3, rng), std::domain_error);
    }

    // Make sure sampler throws errors with these values
    std::vector<double> bad_values = { stan::math::positive_infinity(),
                                       stan::math::negative_infinity(),
                                       stan::math::not_a_number() };

    for (auto bad_value : bad_values) {
      assign_parameter_values(p1, { bad_value });
      assign_parameter_values(p2, good_p2);
      assign_parameter_values(p3, good_p3);
      EXPECT_THROW(rig.generate_samples(p1, p2, p3, rng), std::domain_error);

      if (p2_is_used) {
        assign_parameter_values(p1, good_p1);
        assign_parameter_values(p2, { bad_value });
        EXPECT_THROW(rig.generate_samples(p1, p2, p3, rng), std::domain_error);
      }

      if (p3_is_used) {
        assign_parameter_values(p2, good_p2);
        assign_parameter_values(p3, { bad_value });
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
      EXPECT_THROW(rig.generate_samples(p1, p2, p3, rng), std::invalid_argument);
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
      EXPECT_THROW(rig.generate_samples(p1, p2, p3, rng), std::invalid_argument);
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
      EXPECT_THROW(rig.generate_samples(p1, p2, p3, rng), std::invalid_argument);
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

  call_permutations<std::tuple<double, std::vector<double>, VectorXd,
                               RowVectorXd> >(check_dist_throws{}, rig);
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

    assign_parameter_values(p1, rig.good_p1_);
    assign_parameter_values(p2, rig.good_p2_);
    assign_parameter_values(p3, rig.good_p3_);

    bool p2_is_used = rig.good_p2_.size() > 0;
    bool p3_is_used = rig.good_p3_.size() > 0;

    int M = std::max({ stan::length(p1),
          (p2_is_used) ? stan::length(p2) : 1,
          (p3_is_used) ? stan::length(p3) : 1 });

    stan::scalar_seq_view<T_param1> p1_vec(p1);
    stan::scalar_seq_view<T_param2> p2_vec(p2);
    stan::scalar_seq_view<T_param3> p3_vec(p3);

    std::vector<std::vector<double> > samples_to_test_transpose;
    for (int n = 0; n < rig.N_; ++n) {
      // If p1, p2, and p3 are scalars, the output is a scalar. Need to promote it
      // to a std::vector
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

  call_permutations<std::tuple<double, std::vector<double>, VectorXd,
                               RowVectorXd> >(check_quantiles{}, rig);
}

#endif
