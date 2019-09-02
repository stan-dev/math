#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/scal/meta/is_eigen.hpp>
#include <Eigen/Sparse>
#include <gtest/gtest.h>
#include <vector>

TEST(MathMeta, primitive) {
  using stan::is_eigen;
  EXPECT_FALSE((is_eigen<bool>::value));
  EXPECT_FALSE((is_eigen<double>::value));
  EXPECT_FALSE((is_eigen<int>::value));

  EXPECT_FALSE((is_eigen<std::vector<double>>::value));

  EXPECT_FALSE((is_eigen<Eigen::Matrix<double, -1, -1>>::value));
  EXPECT_FALSE((is_eigen<Eigen::SparseMatrix<double>>::value));
  EXPECT_FALSE(
      (is_eigen<Eigen::MatrixBase<Eigen::Matrix<double, -1, -1>>>::value));
}
