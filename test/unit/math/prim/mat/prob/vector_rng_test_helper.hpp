#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <boost/math/distributions.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <test/unit/math/prim/scal/prob/util.hpp>

enum class Constraint { None, Positive, NonNegative };

/*
 * Resize v to be length N and fill v with values that satisfy the given constraint c.
 *
 * @tparam T1 Type of v
 * @param v Vector to fill with values
 * @param c Constraint that values must fulfill
 * @param N Number of elements that v will be resized to
 */
template <typename T1>
void generateValues(T1& v, Constraint c, int N) {
  v.resize(N);
  
  if (c == Constraint::None) {
    for (int i = 0; i < N; i++)
      v[i] = i - 1;
  } else if (c == Constraint::Positive) {
    for (int i = 0; i < N; i++)
      v[i] = i + 1;
  } else if (c == Constraint::NonNegative) {
    for (int i = 0; i < N; i++)
      v[i] = i;
  }
}

/*
 * Set v to a value that satisfies the constraint (and as much as possible
 * does not satisfy other constraints -- but this isn't totally possible)
 *
 * @tparam T1 Type of v
 * @param v Variable to assign value to
 * @param c Constraint that value must fulfill
 * @param N Unused
 */
template <>
void generateValues(double& v, Constraint c, int N) {
  if (c == Constraint::None) {
    v = -1;
  } else if (c == Constraint::Positive) {
    v = 1;
  } else if (c == Constraint::NonNegative) {
    v = 0;
  }
}

/*
 * Sets first element of v to a given value
 *
 * @tparam T1 type of v
 * @param v Vector to assign values to
 * @param value Value to assign
 */
template <typename T1>
void setValue(T1& v, double value) {
  v[0] = value;
}

/*
 * Sets v equal to value
 *
 * @param v Variable to assign values to
 * @param value Value to assign
 */
template<>
void setValue(double& v, double value) {
  v = value;
}

/*
 * Return true if variables defined with constraint c1 do not violate constraint
 * c2
 *
 * @param c1 First constraint
 * @param c2 Second constraint
 * @return True if values with constraint c1 and compatible with c2
 */
bool compatible(Constraint c1, Constraint c2) {
  if (c1 == Constraint::NonNegative && c2 == Constraint::Positive)
    return false;

  if (c1 == Constraint::None && c2 == Constraint::NonNegative)
    return false;

  if (c1 == Constraint::None && c2 == Constraint::Positive)
    return false;

  return true;
}

/*
 * twoParamCheckTrows feeds the given generate_samples callable various
 * combinations of valid and invalid parameters (as defined by the two
 * constraints, c1 and c2). For generated constraints that are not compatible
 * with the supplied constraints, the generate_samples callable is expected
 * to throw domain_errors.
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
 * @param generate_samples generate_samples functor to test
 * @param c1 Constraint of first parameter
 * @param c2 Constraint of second parameter
 */
template<typename T_param1, typename T_param2, typename T_generate_samples>
void twoParamDistCheckThrows(T_generate_samples generate_samples, Constraint c1,
                             Constraint c2) {
  boost::random::mt19937 rng;

  T_param1 p1;
  T_param2 p2;

  std::vector<Constraint> constraints = { Constraint::Positive,
                                          Constraint::NonNegative,
                                          Constraint::None };

  for (auto it1 = constraints.begin(); it1 != constraints.end(); it1++) {
    for (auto it2 = constraints.begin(); it2 != constraints.end(); it2++) {
      generateValues(p1, *it1, 3);
      generateValues(p2, *it2, 3);

      if (compatible(*it1, c1) && compatible(*it2, c2)) {
        EXPECT_NO_THROW(generate_samples(p1, p2, rng));
      } else {
        EXPECT_THROW(generate_samples(p1, p2, rng), std::domain_error);
      }
    }
  }

  std::vector<double> badValues = { stan::math::positive_infinity(),
                                    stan::math::negative_infinity(),
                                    stan::math::not_a_number() };

  for (auto it = badValues.begin(); it != badValues.end(); it++) {
    generateValues(p1, c1, 3);
    generateValues(p2, c2, 3);
    setValue(p1, *it);
    EXPECT_THROW(generate_samples(p1, p2, rng), std::domain_error);

    generateValues(p1, c1, 3);
    setValue(p2, *it);
    EXPECT_THROW(generate_samples(p1, p2, rng), std::domain_error);
  }

  generateValues(p1, c1, 3);
  generateValues(p2, c2, 4);
  if (stan::length(p1) != 1 && stan::length(p2) != 1)
    EXPECT_THROW(generate_samples(p1, p2, rng), std::invalid_argument);
}

/*
 * Call twoParamDistCheckThrows with the given arguments and all combinations of
 * double, std::vector<double>, Eigen::VectorXd, and Eigen::RowVectorXd
 *
 * @tparam T_generate_samples Type of generate_samples functor
 * @param generate_samples generate_samples functor to test
 * @param c1 Constraint of first parameter
 * @param c2 Constraint of second parameter
 */
template<typename T_generate_samples>
void twoParamDistCheckThrowsAllTypes(T_generate_samples generate_samples, Constraint c1,
                                     Constraint c2) {
  using namespace Eigen;
  
  twoParamDistCheckThrows<double, double>
    (generate_samples, c1, c2);
  twoParamDistCheckThrows<double, std::vector<double> >
    (generate_samples, c1, c2);
  twoParamDistCheckThrows<double, VectorXd>
    (generate_samples, c1, c2);
  twoParamDistCheckThrows<double, RowVectorXd>
    (generate_samples, c1, c2);
  twoParamDistCheckThrows<std::vector<double>, double>
    (generate_samples, c1, c2);
  twoParamDistCheckThrows<std::vector<double>, std::vector<double> >
    (generate_samples, c1, c2);
  twoParamDistCheckThrows<std::vector<double>, VectorXd>
    (generate_samples, c1, c2);
  twoParamDistCheckThrows<std::vector<double>, RowVectorXd>
    (generate_samples, c1, c2);
  twoParamDistCheckThrows<VectorXd, double>
    (generate_samples, c1, c2);
  twoParamDistCheckThrows<VectorXd, std::vector<double> >
    (generate_samples, c1, c2);
  twoParamDistCheckThrows<VectorXd, VectorXd>
    (generate_samples, c1, c2);
  twoParamDistCheckThrows<VectorXd, RowVectorXd>
    (generate_samples, c1, c2);
  twoParamDistCheckThrows<RowVectorXd, double>
    (generate_samples, c1, c2);
  twoParamDistCheckThrows<RowVectorXd, std::vector<double> >
    (generate_samples, c1, c2);
  twoParamDistCheckThrows<RowVectorXd, VectorXd>
    (generate_samples, c1, c2);
  twoParamDistCheckThrows<RowVectorXd, RowVectorXd>
    (generate_samples, c1, c2);
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
 * For arguments that are already vectors, just return the input (via std::move)
 */
std::vector<double>&& promote_to_vector(std::vector<double>&& v) {
  return std::move(v);
}

/*
 * chiSquareTest uses assert_matches_quantiles to check random numbers generated
 * by a given callable generate_samples against the reference quantiles computed
 * with the given callable generate_quantiles.
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
 * c1 and c2.
 *
 * @tparam T_param1 Type of first parameter
 * @tparam T_param2 Type of second parameter
 * @tparam T_generate_samples Type of callable for generating samples
 * @tparam T_generate_quantiles Type of callable for generating quantiles
 * @param sampler sampler functor to test
 * @param c1 Constraint of first parameter
 * @param c2 Constraint of second parameter
 */
template<typename T_param1, typename T_param2,
         typename T_generate_samples, typename T_generate_quantiles>
void chiSquareTest(int N, int M_vec,
                   T_generate_samples generate_samples,
                   T_generate_quantiles generate_quantiles,
                   Constraint c1,
                   Constraint c2) {
  boost::random::mt19937 rng;
  T_param1 p1;
  T_param2 p2;
  int M = stan::max_size(M_vec, M_vec);
  generateValues(p1, c1, M);
  generateValues(p2, c2, M);
  stan::scalar_seq_view<T_param1> p1_vec(p1);
  stan::scalar_seq_view<T_param2> p2_vec(p2);

  std::vector<std::vector<double> > samples_to_test_transpose;
  for (int n = 0; n < N; ++n) {
    // If p1 and p2 are scalars, the output is a scalar. Need to promote it
    // to a std::vector
    samples_to_test_transpose.push_back(
      promote_to_vector(generate_samples(p1, p2, rng))
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
 * Call chiSquareTest with the given arguments and all combinations of
 * double, std::vector<double>, Eigen::VectorXd, and Eigen::RowVectorXd
 *
 * @tparam T_generate_samples Type of callable for generating samples
 * @tparam T_generate_quantiles Type of callable for generating quantiles
 * @param sampler sampler functor to test
 * @param c1 Constraint of first parameter
 * @param c2 Constraint of second parameter
 */
template<typename T_generate_samples, typename T_generate_quantiles>
void chiSquareTestAllTypes(int N, int M,
                           T_generate_samples generate_samples,
                           T_generate_quantiles generate_quantiles,
                           Constraint c1,
                           Constraint c2) {
  using namespace Eigen;
  
  chiSquareTest<double, double>
    (N, M, generate_samples, generate_quantiles, c1, c2);
  chiSquareTest<double, std::vector<double> >
    (N, M, generate_samples, generate_quantiles, c1, c2);
  chiSquareTest<double, VectorXd >
    (N, M, generate_samples, generate_quantiles, c1, c2);
  chiSquareTest<double, RowVectorXd >
    (N, M, generate_samples, generate_quantiles, c1, c2);
  chiSquareTest<std::vector<double>, double>
    (N, M, generate_samples, generate_quantiles, c1, c2);
  chiSquareTest<std::vector<double>, std::vector<double> >
    (N, M, generate_samples, generate_quantiles, c1, c2);
  chiSquareTest<std::vector<double>, VectorXd >
    (N, M, generate_samples, generate_quantiles, c1, c2);
  chiSquareTest<std::vector<double>, RowVectorXd >
    (N, M, generate_samples, generate_quantiles, c1, c2);  
  chiSquareTest<VectorXd, double>
    (N, M, generate_samples, generate_quantiles, c1, c2);
  chiSquareTest<VectorXd, std::vector<double> >
    (N, M, generate_samples, generate_quantiles, c1, c2);
  chiSquareTest<VectorXd, VectorXd>
    (N, M, generate_samples, generate_quantiles, c1, c2);
  chiSquareTest<VectorXd, RowVectorXd>
    (N, M, generate_samples, generate_quantiles, c1, c2);  
  chiSquareTest<RowVectorXd, double>
    (N, M, generate_samples, generate_quantiles, c1, c2);
  chiSquareTest<RowVectorXd, std::vector<double> >
    (N, M, generate_samples, generate_quantiles, c1, c2);
  chiSquareTest<RowVectorXd, VectorXd>
    (N, M, generate_samples, generate_quantiles, c1, c2);
  chiSquareTest<RowVectorXd, RowVectorXd>
    (N, M, generate_samples, generate_quantiles, c1, c2);
}
