#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/scal/meta/is_eigen.hpp>
#include <stan/math/prim/scal.hpp>
#include <Eigen/Sparse>
#include <gtest/gtest.h>
#include <vector>

TEST(MathMetaPrimScal, primitive) {
  using stan::is_eigen;
  EXPECT_FALSE((is_eigen<bool>::value));
  EXPECT_FALSE((is_eigen<double>::value));
  EXPECT_FALSE((is_eigen<int>::value));
}

TEST(MathMetaPrimArr, primitive) {
  using stan::is_eigen;
  EXPECT_FALSE((is_eigen<bool>::value));
  EXPECT_FALSE((is_eigen<double>::value));
  EXPECT_FALSE((is_eigen<int>::value));

  EXPECT_FALSE((is_eigen<std::vector<double>>::value));
}
