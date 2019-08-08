#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>


TEST(MathMatrix, inverse_sparse) {
  using stan::math::inverse;
  using triplet_d = Eigen::Triplet<double>;
  using sparse_mat_d = Eigen::SparseMatrix<double>;
  std::vector<triplet_d> tripletList;
  tripletList.reserve(4);

  tripletList.push_back(triplet_d(0, 0, 0.1));
  tripletList.push_back(triplet_d(0, 1, 0.2));
  tripletList.push_back(triplet_d(1, 0, 0.2));
  tripletList.push_back(triplet_d(1, 1, 1.0));

  sparse_mat_d mat(2, 2);
  mat.setFromTriplets(tripletList.begin(), tripletList.end());
  sparse_mat_d output = inverse(mat);
  EXPECT_FLOAT_EQ(16.6666666666, output.coeff(0, 0));
  EXPECT_FLOAT_EQ(-3.3333333333, output.coeff(1, 0));
  EXPECT_FLOAT_EQ(-3.3333333333, output.coeff(0, 1));
  EXPECT_FLOAT_EQ(1.66666666666, output.coeff(1, 1));
  std::cout << output;
}
