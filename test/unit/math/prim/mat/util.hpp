#ifndef TEST_MATH_UNIT_PRIM_MAT_UTIL_HPP
#define TEST_MATH_UNIT_PRIM_MAT_UTIL_HPP

#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>

namespace stan {
namespace test {
namespace unit {

/**
 * Run a test that fails if the specified square matrix is not
 * symmetric.
 *
 * @param[in] a Matrix to test.
 */
void expect_symmetric(const Eigen::MatrixXd& a) {
  for (int j = 1; j < a.cols(); ++j)
    for (int i = 0; i < j; ++i)
      EXPECT_EQ(a(i, j), a(j, i)) << "failed symmetry at " << i << ", " << j;
}

/**
 * Return a randomly generated symmetric, positive-definite
 * matrix of the specified dimensionality using the specified
 * rng.
 *
 * @tparam RNG Class of random number generator.
 * @param[in] k Number of rows and columns in generated matrix.
 * @param[in, out] rng Random number generator.
 * @return Random k x k symmetric, positive-definite matrix.
 */
template <typename RNG>
Eigen::MatrixXd spd_rng(int k, RNG& rng) {
  using Eigen::MatrixXd;
  using stan::math::normal_rng;
  MatrixXd sigma = MatrixXd::Zero(k, k);
  for (int j = 0; j < k; ++j)
    for (int i = 0; i <= j; ++i)
      sigma(i, j) = normal_rng(0, 1, rng);
  for (int i = 0; i < k; ++i)
    sigma(i, i) *= sigma(i, i);               // pos. diagonal
  sigma = sigma.transpose() * sigma;          // reconstruct full matrix
  sigma = 0.5 * (sigma + sigma.transpose());  // symmetrize
  for (int i = 0; i < k; ++i)
    sigma(i, i) += 5;  // condition
  return sigma;
}
}  // namespace unit
}  // namespace test
}  // namespace stan
#endif
