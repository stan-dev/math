#ifndef TEST_UNIT_MATH_PRIM_PROB_VECTOR_RNG_TEST_RIG_HPP
#define TEST_UNIT_MATH_PRIM_PROB_VECTOR_RNG_TEST_RIG_HPP

#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <vector>

/*
 * VectorRNGTestRig is a class for wrapping up random number
 * generators in a way so that they can be tested to
 *  1) provide interfaces to all the necessary types
 *  2) generate the correct distributions
 *
 * VectorRNGTestRig only allows for interfaces tests. To actually test that a
 * random number generator is generating the right distributions of numbers,
 * VectorRealRNGTestRig and VectorIntRNGTestRig should be used depending on
 * whether or not the random number generator to be tested generates reals or
 * integers
 *
 * A test rig that inherits from VectorRNGTestRig must implement the function
 *  "generate_samples"
 * and must initialize good_p1_, good_p1_int_, bad_p1_, bad_p1_int_, good_p2_,
 * good_p2_int_, bad_p2_, bad_p2_int_, good_p3_, good_p3_int_, and bad_p3_
 * bad_p3_int_ (depending on how many parameters the wrapped random number
 * generator uses)
 *
 * generate_samples is a wrapper around the actual random number generator. It
 * is expected to have the signature:
 *   T_out generate_samples(T_param1 p1, T_param2 p2, T_param3, T_rng)
 * where T_param1, T_param2, and T_param3 can be any combination of:
 *   int, double, std::vector<int>, std::vector<double>, VectorXd, & RowVectorXd
 * and T_rng is a boost random number generator. p1, p2 and p3 are required for
 * the function signature, but they can be ignored if they are not needed.
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
 *   stan::math::size(out) == stan::math::size(p1) if only p1 is used
 *   stan::math::size(out) == stan::math::max_size(p1, p2) if p1 and p2 are used
 *   stan::math::size(out) == stan::math::max_size(p1, p2, p3) if all parameters
 * are used It must be defined to take three arguments, but not all must be
 * used.
 *
 * good_p1_ and bad_p1_ should be initialized to lists of valid and invalid
 * floating point parameters for the first parameter of the tested RNG
 *
 * good_p1_int_ and bad_p1_int_ should be initialized to lists of valid and
 * invalid integer parameters.
 *
 * good_p2_, good_p2_int_, bad_p2_, bad_p2_int_, good_p3_, good_p3_int_, and
 * bad_p3_int_ should be initialized similarly.
 *
 * If testing a distribution that uses only two parameters, leave good_p3_ and
 * good_p3_int as length zero vectors and it will not be tested. Similarly, if
 * good_p2_ and good_p2_int_ have length zero, the second parameter will not be
 * tested either.
 */
class VectorRNGTestRig {
 public:
  int N_;  // Number of samples used in the quantiles tests
  int M_;  // Length of vectors for the vectorization tests

  std::vector<int> test_points_;

  std::vector<double> good_p1_;
  std::vector<int> good_p1_int_;
  std::vector<double> bad_p1_;
  std::vector<int> bad_p1_int_;
  std::vector<double> good_p2_;
  std::vector<int> good_p2_int_;
  std::vector<double> bad_p2_;
  std::vector<int> bad_p2_int_;
  std::vector<double> good_p3_;
  std::vector<int> good_p3_int_;
  std::vector<double> bad_p3_;
  std::vector<int> bad_p3_int_;

  std::vector<double> always_bad_values_
      = {stan::math::positive_infinity(), stan::math::negative_infinity(),
         stan::math::not_a_number()};
  /*
   * This function wraps up the random number generator for testing.
   *
   * The tested rng can have up to three parameters. Any unused parameters can
   * be ignored.
   *
   * It *must* be implemented in the child class (it isn't virtual here because
   * C++ doesn't allow templated virtual member functions)
   */
  template <typename T1, typename T2, typename T3, typename T_rng>
  auto generate_samples(const T1& p1, const T2& p2, const T3& p3,
                        T_rng& rng) const;

  VectorRNGTestRig(int N, int M, std::vector<double> good_p1,
                   std::vector<int> good_p1_int, std::vector<double> bad_p1,
                   std::vector<int> bad_p1_int, std::vector<double> good_p2,
                   std::vector<int> good_p2_int, std::vector<double> bad_p2,
                   std::vector<int> bad_p2_int, std::vector<double> good_p3,
                   std::vector<int> good_p3_int, std::vector<double> bad_p3,
                   std::vector<int> bad_p3_int)
      : N_(N),
        M_(M),
        good_p1_(good_p1),
        good_p1_int_(good_p1_int),
        bad_p1_(bad_p1),
        bad_p1_int_(bad_p1_int),
        good_p2_(good_p2),
        good_p2_int_(good_p2_int),
        bad_p2_(bad_p2),
        bad_p2_int_(bad_p2_int),
        good_p3_(good_p3),
        good_p3_int_(good_p3_int),
        bad_p3_(bad_p3),
        bad_p3_int_(bad_p3_int) {
    if (good_p1.size() > 0 && good_p1_int.size() == 0)
      throw std::domain_error(
          "good_p1 has non-zero length, but good_p1_int "
          "still has length zero (good_p1_int must also "
          "have non-zero length)");

    if (good_p2.size() > 0 && good_p2_int.size() == 0)
      throw std::domain_error(
          "good_p2 has non-zero length, but good_p2_int "
          "still has length zero (good_p2_int must also "
          "have non-zero length)");

    if (good_p3.size() > 0 && good_p3_int.size() == 0)
      throw std::domain_error(
          "good_p3 has non-zero length, but good_p3_int "
          "still has length zero (good_p3_int must also "
          "have non-zero length)");
  }

  VectorRNGTestRig(int N, int M, std::vector<double> good_p1,
                   std::vector<int> good_p1_int, std::vector<double> bad_p1,
                   std::vector<int> bad_p1_int, std::vector<double> good_p2,
                   std::vector<int> good_p2_int, std::vector<double> bad_p2,
                   std::vector<int> bad_p2_int)
      : VectorRNGTestRig(N, M, good_p1, good_p1_int, bad_p1, bad_p1_int,
                         good_p2, good_p2_int, bad_p2, bad_p2_int, {}, {}, {},
                         {}) {}

  VectorRNGTestRig(int N, int M, std::vector<double> good_p1,
                   std::vector<int> good_p1_int, std::vector<double> bad_p1,
                   std::vector<int> bad_p1_int)
      : VectorRNGTestRig(N, M, good_p1, good_p1_int, bad_p1, bad_p1_int, {}, {},
                         {}, {}, {}, {}, {}, {}) {}

  VectorRNGTestRig(int N, int M)
      : VectorRNGTestRig(N, M, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}) {
  }

  /*
   * If no good values for p1, p2 or p3 are provided, it is assumed that those
   * parameters are unused
   */
  bool p1_is_used() const {
    return good_p1_.size() > 0 || good_p1_int_.size() > 0;
  }

  bool p2_is_used() const {
    return good_p2_.size() > 0 || good_p2_int_.size() > 0;
  }

  bool p3_is_used() const {
    return good_p3_.size() > 0 || good_p3_int_.size() > 0;
  }

  /*
   * All the accessors work the same. Based on the input template type, they
   * return the list of valid/invalid values for whichever parameter
   * is accessed.
   *
   * The bad value accessors always append on the list of always_bad_values
   * to whatever was already specified (always_bad_values contains infs and
   * nans and stuff that should always cause errors).
   *
   * @tparam T Template parameter determining which list of parameters to return
   * (can be double or int)
   * @return List of parameter values
   */
  template <typename T>
  std::vector<T> get_good_p1() const {
    return good_p1_;
  }

  template <typename T>
  std::vector<T> get_bad_p1() const {
    return stan::math::append_array(bad_p1_, always_bad_values_);
  }

  template <typename T>
  std::vector<T> get_good_p2() const {
    return good_p2_;
  }

  template <typename T>
  std::vector<T> get_bad_p2() const {
    return stan::math::append_array(bad_p2_, always_bad_values_);
  }

  template <typename T>
  std::vector<T> get_good_p3() const {
    return good_p3_;
  }

  template <typename T>
  std::vector<T> get_bad_p3() const {
    return stan::math::append_array(bad_p3_, always_bad_values_);
  }
};

template <>
std::vector<int> VectorRNGTestRig::get_good_p1<int>() const {
  return good_p1_int_;
}

template <>
std::vector<int> VectorRNGTestRig::get_bad_p1<int>() const {
  return bad_p1_int_;
}

template <>
std::vector<int> VectorRNGTestRig::get_good_p2<int>() const {
  return good_p2_int_;
}

template <>
std::vector<int> VectorRNGTestRig::get_bad_p2<int>() const {
  return bad_p2_int_;
}

template <>
std::vector<int> VectorRNGTestRig::get_good_p3<int>() const {
  return good_p3_int_;
}

template <>
std::vector<int> VectorRNGTestRig::get_bad_p3<int>() const {
  return bad_p3_int_;
}

#endif
