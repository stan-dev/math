#ifndef TEST_UNIT_MATH_PRIM_PROB_VECTOR_REAL_RNG_TEST_RIG_HPP
#define TEST_UNIT_MATH_PRIM_PROB_VECTOR_REAL_RNG_TEST_RIG_HPP

#include <test/unit/math/prim/prob/VectorRNGTestRig.hpp>
#include <vector>

/*
 * VectorRealRNGTestRig is a subclass of VectorRNGTestRig specialized for
 * testing random numbers generated from continuous distributions
 *
 * A test rig that inherits from VectorRealRNGTestRig must implement the
 * function "generate_quantiles" on top of the requirements listed in the docs
 * for VectorRNGTestRig.
 *
 * The generate_quantiles callable must have the signature:
 *  std::vector<double> generate_samples(double p1, double p2, double p3)
 * (similarly to generate_samples in VectorRNGTestRig, ignore parameters p2 and
 * p3 if they aren't necessary)
 *
 * generate_quantiles should compute the quantiles that will be compared against
 * (using assert_matches_quantiles) empirical distributions generated with
 * generate_samples (which will be of length N_ for each parameter combination).
 *
 * p1, p2, and p3 are the values of the parameters that these quantiles
 * correspond to.
 */
class VectorRealRNGTestRig : public VectorRNGTestRig {
 public:
  /*
   * This function builds the quantiles that we will supply to
   * assert_matches_quantiles to test the random number generator implemented in
   * generate_samples
   */
  virtual std::vector<double> generate_quantiles(double p1, double p2,
                                                 double p3) const = 0;

  VectorRealRNGTestRig(int N, int M, std::vector<double> good_p1,
                       std::vector<int> good_p1_int, std::vector<double> bad_p1,
                       std::vector<int> bad_p1_int, std::vector<double> good_p2,
                       std::vector<int> good_p2_int, std::vector<double> bad_p2,
                       std::vector<int> bad_p2_int, std::vector<double> good_p3,
                       std::vector<int> good_p3_int, std::vector<double> bad_p3,
                       std::vector<int> bad_p3_int)
      : VectorRNGTestRig(N, M, good_p1, good_p1_int, bad_p1, bad_p1_int,
                         good_p2, good_p2_int, bad_p2, bad_p2_int, good_p3,
                         good_p3_int, bad_p3, bad_p3_int) {}

  VectorRealRNGTestRig(int N, int M, std::vector<double> good_p1,
                       std::vector<int> good_p1_int, std::vector<double> bad_p1,
                       std::vector<int> bad_p1_int, std::vector<double> good_p2,
                       std::vector<int> good_p2_int, std::vector<double> bad_p2,
                       std::vector<int> bad_p2_int)
      : VectorRNGTestRig(N, M, good_p1, good_p1_int, bad_p1, bad_p1_int,
                         good_p2, good_p2_int, bad_p2, bad_p2_int, {}, {}, {},
                         {}) {}

  VectorRealRNGTestRig(int N, int M, std::vector<double> good_p1,
                       std::vector<int> good_p1_int, std::vector<double> bad_p1,
                       std::vector<int> bad_p1_int)
      : VectorRNGTestRig(N, M, good_p1, good_p1_int, bad_p1, bad_p1_int, {}, {},
                         {}, {}, {}, {}, {}, {}) {}

  VectorRealRNGTestRig(int N, int M)
      : VectorRNGTestRig(N, M, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}) {
  }
};

#endif
