#include <stan/math/prim/mat.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(probTransform, factorCovMatrix) {
  using Eigen::Array;
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::factor_cov_matrix;

  Matrix<double, Dynamic, Dynamic> L(3, 3);
  L << 1.7, 0, 0, -2.9, 14.2, 0, .2, -.5, 1.3;

  Matrix<double, Dynamic, Dynamic> Sigma = L.transpose() * L;

  Array<double, Dynamic, 1> CPCs(3);  // must be sized coming in
  Array<double, Dynamic, 1> sds(3);   // must be sized coming in

  // just check it doesn't bomb
  factor_cov_matrix(Sigma, CPCs, sds);

  // example of sizing for K=2
  L.resize(2, 2);
  L << 1.7, 0, -2.3, 0.5;
  Sigma = L.transpose() * L;
  CPCs.resize(1);
  sds.resize(2);
  factor_cov_matrix(Sigma, CPCs, sds);
}
