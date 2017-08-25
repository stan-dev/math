#include <gtest/gtest.h>
#include <boost/math/distributions.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <stan/math/prim/mat.hpp>
#include <test/unit/math/prim/scal/prob/util.hpp>

/*
 * Objects of this class can be used to generate sequences of non-random
 * numbers.
 */
class NonRandomGenerator {
private:
  int last_used_idx_;

  const std::vector<double> real_values_ =
    { -5.0, -4.0, -3.5, -1.0, -0.5, 0.0, 1.5, 2.8, 3.0, 5.0, 7.0 };

public:
  NonRandomGenerator() : last_used_idx_(0) {}

  size_t size() const {
    return real_values_.size();
  }

  /*
   * Get next value from non-random generator. Increment internal index to the
   * next valid value (wrapping around to the beginning).
   *
   * @return A non-random number
   */
  double operator()() {
    last_used_idx_ = (last_used_idx_ + 1) % real_values_.size();
    return real_values_[last_used_idx_];
  }

  /*
   * Test to make sure that it is possible to generate at least one value with
   * the NonRandomGenerator that satisfies the given constraint.
   *
   * @tparam T_check Callable that takes in one double and returns a bool
   * @param check Callable that checks constraint
   */
  template<typename T_check>
  static void CanSatisfyConstraint(T_check check) {
    NonRandomGenerator non_random_generator;

    bool at_least_one_value_works = false;
    for(size_t i = 0; i < non_random_generator.size(); i++)
      at_least_one_value_works |= check(non_random_generator());
    if(!at_least_one_value_works)
      throw std::domain_error("None of the test values statisfy the constraint");
  }
};

/*
 * Build a length N vector-like variable of type T1. Fill it with numbers
 * generated from the callable non_random_generator that satisfy the constraint
 * check.
 *
 * Before calling this functor, check that the non_random_generator can
 * actually generate values that satisfy the check constraint, otherwise this
 * function will loop endlessly.
 *
 * This is a functor instead of a function because it takes advantage of partial
 * template specialization.
 *
 * @tparam T1 Type of v
 * @tparam T_check Callable type that takes in one double and returns a bool
 * @tparam T_generator Callable type that takes no arguments and returns a double
 * @param N Number of elements that v will be resized to
 * @param check Callable for validating p
 * @param non_random_generator Callable for generating possible parameters
 */
template<typename T1>
struct GenerateValues {
  template <typename T_check, typename T_generator>
  T1 operator()(int N, T_check check, T_generator &non_random_generator) {
    T1 vec(N);
    
    for (int i = 0; i < vec.size(); i++) {
      double v;
    
      do {
        v = non_random_generator();
      } while(!check(v));
    
      vec[i] = v;
    }

    return vec;
  }
};

/*
 * Generate a value from non_random_generator that satisfies the constraint
 * check.
 *
 * Before calling this functor, check that the non_random_generator can
 * actually generate values that satisfy the check constraint, otherwise this
 * function will loop endlessly.
 *
 * This is a functor instead of a function because it takes advantage of partial
 * template specialization.
 *
 * @tparam T1 Type of v
 * @tparam T_check Callable type that takes in one double and returns a bool
 * @tparam T_generator Callable type that takes no arguments and returns a
 * double
 * @param N Unused
 * @param check Callable for validating p
 * @param non_random_generator Callable for generating possible parameters
 */
template<>
struct GenerateValues<double> {
  template <typename T_check, typename T_generator>
  double operator()(int N, T_check check, T_generator &non_random_generator) {
    double v;
    
    do {
      v = non_random_generator();
    } while(!check(v));

    return v;
  }
};

/*
 * Sets first element of v to a given value
 *
 * @tparam T1 type of v
 * @param v Vector to assign values to
 * @param value Value to assign
 */
template <typename T1>
void SetValue(T1& v, double value) {
  v[0] = value;
}

/*
 * Sets v equal to value
 *
 * @param v Variable to assign values to
 * @param value Value to assign
 */
template<>
void SetValue(double& v, double value) {
  v = value;
}

/*
 * Check that every element of p (if p is a vector) passes check. If p is a
 * scalar, then just make sure it passes the check.
 *
 * @tparam T_check Callable type that takes in one double and returns a bool
 * @tparam T_param Any type compatible with stan::scalar_seq_view
 * @param check Callable for validating (values of) p
 * @param p Parameter to check
 */
template<typename T_check, typename T_param>
bool VectorizeCheck(T_check check, T_param p) {
  stan::scalar_seq_view<T_param> p_vec(p);

  bool out = true;
  for (int i = 0; i < p_vec.size(); i++) {
    out &= check(p_vec[i]);
  }

  return out;
}

/*
 * CheckDistThrows feeds the given generate_samples callable various
 * combinations of valid and invalid parameters (as defined by the two
 * constraints callables, check_p1 and check_p2). For generated parameters
 * that cause check_p1 or check_p2 to return false, the generate_samples
 * callable is expected to throw domain_errors. If both checks pass, the sampler
 * is expected to throw no errors.
 *
 * The generate_samples callable is also passed various other guaranteed invalid
 * values like positive infinity, negative infinity, and NaNs (these should also
 * cause domain_errors).
 *
 * The generated_samples callable is also tested to reject incompatibly sized
 * input vectors. These should cause invalid_argument errors.
 *
 * The generate_samples callable must have the signature:
 *   T_out generate_samples(T_param1 p1, T_param2 p2, T_rng)
 * where T_param1 and T_param2 can be any combination of:
 *   double, std::vector<double>, VectorXd, or RowVectorXd
 * and T_rng is a boost random number generator
 *
 * The output of generate_samples must be of
 *   stan::length(out) == stan::max_size(p1, p2)
 *
 * Each element of the output, out[i], should correspond to a random variate
 * generated according to parameters p1[i] and p2[i]
 *
 * @tparam T_param1 Type of first parameter
 * @tparam T_param2 Type of second parameter
 * @tparam T_generate_samples Type of generate_samples functor
 * @tparam T_check_p1 Callable type that takes in one double and returns a bool
 * @tparam T_check_p2 Callable type that takes in one double and returns a bool
 * @param generate_samples generate_samples functor to test
 * @param check_p1 Callable for validating first parameter
 * @param check_p2 Callable for validating second parameter
 */
template<typename T_param1, typename T_param2, typename T_generate_samples,
         typename T_check_p1, typename T_check_p2>
void CheckDistThrows(T_generate_samples generate_samples,
                     T_check_p1 check_p1, T_check_p2 check_p2) {
  boost::random::mt19937 rng;
  NonRandomGenerator non_random_generator;
  
  T_param1 p1;
  T_param2 p2;

  // Try a few combinations of parameters that should work
  for (int i = 0; i < 5; i++) {
    p1 = GenerateValues<T_param1>{}(5, check_p1, non_random_generator);
    p2 = GenerateValues<T_param2>{}(5, check_p2, non_random_generator);
    
    EXPECT_NO_THROW(generate_samples(p1, p2, rng));
  }

  // Now just try putting anything for the first parameter and checking the
  // throws
  for (int i = 0; i < 5; i++) {
    p1 = GenerateValues<T_param1>{}(5, [](double mean) { return true; },
                                    non_random_generator);
    p2 = GenerateValues<T_param2>{}(5, check_p2, non_random_generator);

    if (VectorizeCheck(check_p1, p1) && VectorizeCheck(check_p2, p2)) {
      EXPECT_NO_THROW(generate_samples(p1, p2, rng));
    } else {
      EXPECT_THROW(generate_samples(p1, p2, rng), std::domain_error);
    }
  }

  // Now just try putting anything for the second parameter and checking the
  // throws
  for (int i = 0; i < 5; i++) {
    p1 = GenerateValues<T_param1>{}(5, check_p1, non_random_generator);
    p2 = GenerateValues<T_param2>{}(5, [](double sd) { return sd > 0.0; },
                                    non_random_generator);

    if (VectorizeCheck(check_p1, p1) && VectorizeCheck(check_p2, p2)) {
      EXPECT_NO_THROW(generate_samples(p1, p2, rng));
    } else {
      EXPECT_THROW(generate_samples(p1, p2, rng), std::domain_error);
    }
  }

  // Make sure sampler throws errors with these values
  std::vector<double> bad_values = { stan::math::positive_infinity(),
                                     stan::math::negative_infinity(),
                                     stan::math::not_a_number() };

  for (auto bad_value : bad_values) {
    p1 = GenerateValues<T_param1>{}(5, check_p1, non_random_generator);
    p2 = GenerateValues<T_param2>{}(5, check_p2, non_random_generator);
    SetValue(p1, bad_value); // p1 is written here
    EXPECT_THROW(generate_samples(p1, p2, rng), std::domain_error);

    p1 = GenerateValues<T_param1>{}(5, check_p1, non_random_generator);
    SetValue(p2, bad_value); // p2 is written here
    EXPECT_THROW(generate_samples(p1, p2, rng), std::domain_error);
  }

  // Check to make sure _rng rejects vector arguments of different lengths
  // If either p1 or p2 are scalars, this test is skipped
  p1 = GenerateValues<T_param1>{}(5, check_p1, non_random_generator);
  p2 = GenerateValues<T_param2>{}(3, check_p2, non_random_generator);
  if (stan::length(p1) != 1 && stan::length(p2) != 1)
    EXPECT_THROW(generate_samples(p1, p2, rng), std::invalid_argument);
}

/*
 * This function calls CheckDistThrows with the given generate_samples
 * callable and all combinations of double, std::vector<double>,
 * Eigen::VectorXd, and Eigen::RowVectorXd as template arguments
 *
 * @tparam T_generate_samples Type of generate_samples functor
 * @tparam T_check_p1 Callable type that takes in one double and returns a bool
 * @tparam T_check_p2 Callable type that takes in one double and returns a bool
 * @param generate_samples generate_samples functor to test
 * @param c1 Callable for validating first parameter
 * @param c2 Callable for validating second parameter
 */
template<typename T_generate_samples, typename T_check_p1, typename T_check_p2>
void CheckDistThrowsAllTypes(T_generate_samples generate_samples,
                             T_check_p1 check_p1, T_check_p2 check_p2) {
  using namespace Eigen;
  // Check that given constraint can be satisfied, otherwise tests will deadlock
  NonRandomGenerator::CanSatisfyConstraint(check_p1);
  NonRandomGenerator::CanSatisfyConstraint(check_p2);
  
  CheckDistThrows<double, double>
    (generate_samples, check_p1, check_p2);
  CheckDistThrows<double, std::vector<double> >
    (generate_samples, check_p1, check_p2);
  CheckDistThrows<double, VectorXd>
    (generate_samples, check_p1, check_p2);
  CheckDistThrows<double, RowVectorXd>
    (generate_samples, check_p1, check_p2);
  CheckDistThrows<std::vector<double>, double>
    (generate_samples, check_p1, check_p2);
  CheckDistThrows<std::vector<double>, std::vector<double> >
    (generate_samples, check_p1, check_p2);
  CheckDistThrows<std::vector<double>, VectorXd>
    (generate_samples, check_p1, check_p2);
  CheckDistThrows<std::vector<double>, RowVectorXd>
    (generate_samples, check_p1, check_p2);
  CheckDistThrows<VectorXd, double>
    (generate_samples, check_p1, check_p2);
  CheckDistThrows<VectorXd, std::vector<double> >
    (generate_samples, check_p1, check_p2);
  CheckDistThrows<VectorXd, VectorXd>
    (generate_samples, check_p1, check_p2);
  CheckDistThrows<VectorXd, RowVectorXd>
    (generate_samples, check_p1, check_p2);
  CheckDistThrows<RowVectorXd, double>
    (generate_samples, check_p1, check_p2);
  CheckDistThrows<RowVectorXd, std::vector<double> >
    (generate_samples, check_p1, check_p2);
  CheckDistThrows<RowVectorXd, VectorXd>
    (generate_samples, check_p1, check_p2);
  CheckDistThrows<RowVectorXd, RowVectorXd>
    (generate_samples, check_p1, check_p2);
}

/*
 * Convert a scalar to a length one vector
 *
 * @param v Input scalar
 * @return vector of length 1 with value v
 */
std::vector<double> PromoteToVector(double v) {
  return std::vector<double>(1, v);
}

/*
 * For arguments that are already vectors, just return the input (via std::move)
 */
std::vector<double>&& PromoteToVector(std::vector<double>&& v) {
  return std::move(v);
}

/*
 * CheckQuantiles uses assert_matches_quantiles to check random numbers
 * generated by the callable generate_samples against the reference
 * quantiles computed with the given callable generate_quantiles.
 *
 * The generate_samples callable must have the signature:
 *   T_out generate_samples(T_param1 p1, T_param2 p2, T_rng)
 * where T_param1 and T_param2 can be any combination of:
 *   double, std::vector<double>, VectorXd, or RowVectorXd
 * and T_rng is a boost random number generator
 *
 * The output of generate_samples must be of
 *   stan::length(out) == stan::max_size(p1, p2)
 *
 * Each element of the output, out[i], should correspond to a random variate
 * generated according to parameters p1[i] and p2[i]
 *
 * The generate_quantiles callable must have the signature:
 *   std::vector<double> generate_samples(int N, double p1, double p2)
 *
 * N is the number of samples that the quantiles will be compared
 * against and p1 and p2 the values of the parameters that these quantiles
 * correspond to.
 *
 * The output of generate_quantiles should be the quantiles used for comparing
 * empirical distributions used by assert_matches_quantiles
 *
 * Parameters for each of these callables are generated according to constraints
 * callables check_p1 and check_p2.
 *
 * @tparam T_param1 Type of first parameter
 * @tparam T_param2 Type of second parameter
 * @tparam T_generate_samples Type of callable for generating samples
 * @tparam T_generate_quantiles Type of callable for generating quantiles
 * @tparam T_check_p1 Callable type that takes in one double and returns a bool
 * @tparam T_check_p2 Callable type that takes in one double and returns a bool
 * @param sampler sampler functor to test
 * @param check_p1 Callable for validating first parameter
 * @param check_p2 Callable for validating second parameter
 */
template<typename T_param1, typename T_param2,
         typename T_generate_samples, typename T_generate_quantiles,
         typename T_check_p1, typename T_check_p2>
void CheckQuantiles(int N, int M_vec,
                    T_generate_samples generate_samples,
                    T_generate_quantiles generate_quantiles,
                    T_check_p1 check_p1, T_check_p2 check_p2) {
  boost::random::mt19937 rng;
  NonRandomGenerator non_random_generator;
  T_param1 p1;
  T_param2 p2;
  int M = stan::max_size(M_vec, M_vec);
  p1 = GenerateValues<T_param1>{}(M, check_p1, non_random_generator);
  p2 = GenerateValues<T_param2>{}(M, check_p2, non_random_generator);
  stan::scalar_seq_view<T_param1> p1_vec(p1);
  stan::scalar_seq_view<T_param2> p2_vec(p2);

  std::vector<std::vector<double> > samples_to_test_transpose;
  for (int n = 0; n < N; ++n) {
    // If p1 and p2 are scalars, the output is a scalar. Need to promote it
    // to a std::vector
    samples_to_test_transpose.push_back(
      PromoteToVector(generate_samples(p1, p2, rng))
    );
  }
  for (int m = 0; m < M; ++m) {
    std::vector<double> samples_to_test;
    for (int n = 0; n < N; ++n) {
      samples_to_test.push_back(samples_to_test_transpose[n][m]);
    }
    std::vector<double> quantiles = generate_quantiles(N, p1_vec[m], p2_vec[m]);

    assert_matches_quantiles(samples_to_test, quantiles, 1e-6);
  }
}

/*
 * Call CheckQuantiles with the given arguments and all combinations of
 * double, std::vector<double>, Eigen::VectorXd, and Eigen::RowVectorXd
 *
 * @tparam T_generate_samples Type of callable for generating samples
 * @tparam T_generate_quantiles Type of callable for generating quantiles
 * @tparam T_check_p1 Callable type that takes in one double and returns a bool
 * @tparam T_check_p2 Callable type that takes in one double and returns a bool
 * @param sampler sampler functor to test
 * @param c1 Callable for validating first parameter
 * @param c2 Callable for validating second parameter
 */
template<typename T_generate_samples, typename T_generate_quantiles,
         typename T_check_p1, typename T_check_p2>
void CheckQuantilesAllTypes(int N, int M,
                            T_generate_samples generate_samples,
                            T_generate_quantiles generate_quantiles,
                            T_check_p1 check_p1, T_check_p2 check_p2) {
  using namespace Eigen;
  // Check that given constraint can be satisfied, otherwise tests will deadlock
  NonRandomGenerator::CanSatisfyConstraint(check_p1);
  NonRandomGenerator::CanSatisfyConstraint(check_p2);
  
  CheckQuantiles<double, double>
    (N, M, generate_samples, generate_quantiles, check_p1, check_p2);
  CheckQuantiles<double, std::vector<double> >
    (N, M, generate_samples, generate_quantiles, check_p1, check_p2);
  CheckQuantiles<double, VectorXd >
    (N, M, generate_samples, generate_quantiles, check_p1, check_p2);
  CheckQuantiles<double, RowVectorXd >
    (N, M, generate_samples, generate_quantiles, check_p1, check_p2);
  CheckQuantiles<std::vector<double>, double>
    (N, M, generate_samples, generate_quantiles, check_p1, check_p2);
  CheckQuantiles<std::vector<double>, std::vector<double> >
    (N, M, generate_samples, generate_quantiles, check_p1, check_p2);
  CheckQuantiles<std::vector<double>, VectorXd >
    (N, M, generate_samples, generate_quantiles, check_p1, check_p2);
  CheckQuantiles<std::vector<double>, RowVectorXd >
    (N, M, generate_samples, generate_quantiles, check_p1, check_p2);  
  CheckQuantiles<VectorXd, double>
    (N, M, generate_samples, generate_quantiles, check_p1, check_p2);
  CheckQuantiles<VectorXd, std::vector<double> >
    (N, M, generate_samples, generate_quantiles, check_p1, check_p2);
  CheckQuantiles<VectorXd, VectorXd>
    (N, M, generate_samples, generate_quantiles, check_p1, check_p2);
  CheckQuantiles<VectorXd, RowVectorXd>
    (N, M, generate_samples, generate_quantiles, check_p1, check_p2);  
  CheckQuantiles<RowVectorXd, double>
    (N, M, generate_samples, generate_quantiles, check_p1, check_p2);
  CheckQuantiles<RowVectorXd, std::vector<double> >
    (N, M, generate_samples, generate_quantiles, check_p1, check_p2);
  CheckQuantiles<RowVectorXd, VectorXd>
    (N, M, generate_samples, generate_quantiles, check_p1, check_p2);
  CheckQuantiles<RowVectorXd, RowVectorXd>
    (N, M, generate_samples, generate_quantiles, check_p1, check_p2);
}
