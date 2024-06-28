#ifdef STAN_OPENCL
#include <stan/math/rev/core.hpp>
#include <stan/math/opencl/rev/vari.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(AgradRev, matrix_cl_vari_block) {
  using stan::math::vari_value;
  Eigen::MatrixXd a = Eigen::MatrixXd::Random(3, 3);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(3, 3);
  stan::math::matrix_cl<double> a_cl(a);
  stan::math::matrix_cl<double> b_cl(b);

  vari_value<stan::math::matrix_cl<double>> A(a);
  EXPECT_MATRIX_EQ(a.block(0, 1, 2, 2),
                   stan::math::from_matrix_cl(A.block(0, 1, 2, 2).val_));
  vari_value<stan::math::matrix_cl<double>> B(a_cl);
  B.adj_ = b_cl;
  EXPECT_MATRIX_EQ(a.block(0, 1, 2, 2),
                   stan::math::from_matrix_cl(B.block(0, 1, 2, 2).val_));
  EXPECT_MATRIX_EQ(b.block(0, 1, 2, 2),
                   stan::math::from_matrix_cl(B.block(0, 1, 2, 2).adj_));
  vari_value<stan::math::matrix_cl<double>> C(a_cl, a_cl);
  EXPECT_MATRIX_EQ(a.block(0, 1, 2, 2),
                   stan::math::from_matrix_cl(C.block(0, 1, 2, 2).val_));
  EXPECT_MATRIX_EQ(a.block(0, 1, 2, 2),
                   stan::math::from_matrix_cl(C.block(0, 1, 2, 2).adj_));
}

#endif
