#include <stan/math/prim/mat/meta/operands_and_partials.hpp>
#include <gtest/gtest.h>

TEST(AgradPartialsVari, OperandsAndPartialsPrimMatZero) {
  using stan::math::detail::zero_vec_or_mat;
  using Eigen::Matrix;
  using Eigen::VectorXd;
  using Eigen::VectorXi;
  using Eigen::RowVectorXd;
  using Eigen::RowVectorXi;

  typedef Matrix<double, -1, -1> Mouble;
  typedef Matrix<int, -1, -1> Mint;

  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 5; ++j) {
      Mouble m(i, j);
      Mint p = zero_vec_or_mat<int, Mouble, -1, -1>::zero(m);
      for (int x = 0; x < i; ++x) {
        for (int y = 0; y < j; ++y) {
          EXPECT_EQ(0, p(x, y));
        }
      }
    }
  }

  VectorXd vd(5);
  VectorXi vi = zero_vec_or_mat<int, VectorXd, -1, 1>::zero(vd);
  for (int i = 0; i < 5; ++i) {
    EXPECT_EQ(0, vi(i));
  }

  RowVectorXd rvd(5);
  RowVectorXi rvi = zero_vec_or_mat<int, RowVectorXd, 1, -1>::zero(vd);
  for (int i = 0; i < 5; ++i) {
    EXPECT_EQ(0, rvi(i));
  }
}
