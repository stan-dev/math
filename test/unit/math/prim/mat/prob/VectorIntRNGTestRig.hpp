#ifndef TEST_UNIT_MATH_PRIM_MAT_PROB_VECTOR_INT_RNG_TEST_RIG_HPP
#define TEST_UNIT_MATH_PRIM_MAT_PROB_VECTOR_INT_RNG_TEST_RIG_HPP

#include <test/unit/math/prim/mat/prob/VectorRNGTestRig.hpp>
#include <vector>

/*
 * VectorIntRNGTestRig is a subclass of VectorRNGTestRig specialized for
 * testing random numbers generated from continuous distributions.
 *
 * A test rig that inherits from VectorIntRNGTestRig must implement the
 * function "pmf" and supply values for test_points_ on top of the requirements
 * listed in the docs for VectorRNGTestRig.
 *
 * The pmf callable should wrap the probability mass function of the
 * distribution to be test. It is possible that pmf will be called with values
 * outside the support of the test distribution. In this case, pmf should return
 * a zero.
 *
 * test_points_ should be initialized to a list of values where the distribution
 * will be tested.
 */
class VectorIntRNGTestRig : public VectorRNGTestRig {
 public:
  std::vector<int> test_points_;

  /*
   * The first argument is templated to allow for distributions that take an
   * integer first argument (the various binomial distributions) as well as ones
   * that don't (the Poisson distribution)
   *
   * It *must* be implemented in the child class (it isn't virtual here because
   * C++ doesn't allow templated virtual member functions)
   */
  template <typename T1>
  double pmf(int y, T1 p1, double p2, double p3) const;

  VectorIntRNGTestRig(int N, int M, std::vector<int> test_points,
                      std::vector<double> good_p1, std::vector<int> good_p1_int,
                      std::vector<double> bad_p1, std::vector<int> bad_p1_int,
                      std::vector<double> good_p2, std::vector<int> good_p2_int,
                      std::vector<double> bad_p2, std::vector<int> bad_p2_int,
                      std::vector<double> good_p3, std::vector<int> good_p3_int,
                      std::vector<double> bad_p3, std::vector<int> bad_p3_int)
      : VectorRNGTestRig(N, M, good_p1, good_p1_int, bad_p1, bad_p1_int,
                         good_p2, good_p2_int, bad_p2, bad_p2_int, good_p3,
                         good_p3_int, bad_p3, bad_p3_int),
        test_points_(test_points) {}

  VectorIntRNGTestRig(int N, int M, std::vector<int> test_points,
                      std::vector<double> good_p1, std::vector<int> good_p1_int,
                      std::vector<double> bad_p1, std::vector<int> bad_p1_int,
                      std::vector<double> good_p2, std::vector<int> good_p2_int,
                      std::vector<double> bad_p2, std::vector<int> bad_p2_int)
      : VectorIntRNGTestRig(N, M, test_points, good_p1, good_p1_int, bad_p1,
                            bad_p1_int, good_p2, good_p2_int, bad_p2,
                            bad_p2_int, {}, {}, {}, {}) {}

  VectorIntRNGTestRig(int N, int M, std::vector<int> test_points,
                      std::vector<double> good_p1, std::vector<int> good_p1_int,
                      std::vector<double> bad_p1, std::vector<int> bad_p1_int)
      : VectorIntRNGTestRig(N, M, test_points, good_p1, good_p1_int, bad_p1,
                            bad_p1_int, {}, {}, {}, {}, {}, {}, {}, {}) {}
};

#endif
