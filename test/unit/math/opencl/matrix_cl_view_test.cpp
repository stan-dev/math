#ifdef STAN_OPENCL
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <gtest/gtest.h>

TEST(matrix_cl_view, either) {
  using stan::math::matrix_cl_view;
  EXPECT_EQ(matrix_cl_view::Lower,
            either(matrix_cl_view::Lower, matrix_cl_view::Lower));
  EXPECT_EQ(matrix_cl_view::Entire,
            either(matrix_cl_view::Lower, matrix_cl_view::Upper));
  EXPECT_EQ(matrix_cl_view::Lower,
            either(matrix_cl_view::Lower, matrix_cl_view::Diagonal));

  EXPECT_EQ(matrix_cl_view::Upper,
            either(matrix_cl_view::Upper, matrix_cl_view::Upper));
  EXPECT_EQ(matrix_cl_view::Upper,
            either(matrix_cl_view::Upper, matrix_cl_view::Diagonal));

  EXPECT_EQ(matrix_cl_view::Entire,
            either(matrix_cl_view::Entire, matrix_cl_view::Upper));
  EXPECT_EQ(matrix_cl_view::Entire,
            either(matrix_cl_view::Entire, matrix_cl_view::Lower));
  EXPECT_EQ(matrix_cl_view::Entire,
            either(matrix_cl_view::Entire, matrix_cl_view::Diagonal));
  EXPECT_EQ(matrix_cl_view::Entire,
            either(matrix_cl_view::Entire, matrix_cl_view::Entire));
}

TEST(matrix_cl_view, operator_multiply) {
  using stan::math::matrix_cl_view;
  EXPECT_EQ(matrix_cl_view::Diagonal,
            both(matrix_cl_view::Upper, matrix_cl_view::Diagonal));
  EXPECT_EQ(matrix_cl_view::Diagonal,
            both(matrix_cl_view::Entire, matrix_cl_view::Diagonal));
  EXPECT_EQ(matrix_cl_view::Diagonal,
            both(matrix_cl_view::Lower, matrix_cl_view::Diagonal));
  EXPECT_EQ(matrix_cl_view::Diagonal,
            both(matrix_cl_view::Lower, matrix_cl_view::Upper));

  EXPECT_EQ(matrix_cl_view::Lower,
            both(matrix_cl_view::Lower, matrix_cl_view::Lower));
  EXPECT_EQ(matrix_cl_view::Lower,
            both(matrix_cl_view::Entire, matrix_cl_view::Lower));

  EXPECT_EQ(matrix_cl_view::Upper,
            both(matrix_cl_view::Upper, matrix_cl_view::Upper));
  EXPECT_EQ(matrix_cl_view::Upper,
            both(matrix_cl_view::Entire, matrix_cl_view::Upper));
  EXPECT_EQ(matrix_cl_view::Entire,
            both(matrix_cl_view::Entire, matrix_cl_view::Entire));
}

TEST(matrix_cl_view, contains_nonzero) {
  using stan::math::contains_nonzero;
  using stan::math::matrix_cl_view;
  EXPECT_EQ(false,
            contains_nonzero(matrix_cl_view::Upper, matrix_cl_view::Diagonal));
  EXPECT_EQ(false,
            contains_nonzero(matrix_cl_view::Entire, matrix_cl_view::Diagonal));
  EXPECT_EQ(false,
            contains_nonzero(matrix_cl_view::Lower, matrix_cl_view::Diagonal));
  EXPECT_EQ(false,
            contains_nonzero(matrix_cl_view::Lower, matrix_cl_view::Upper));

  EXPECT_EQ(true,
            contains_nonzero(matrix_cl_view::Lower, matrix_cl_view::Lower));
  EXPECT_EQ(true,
            contains_nonzero(matrix_cl_view::Entire, matrix_cl_view::Lower));

  EXPECT_EQ(true,
            contains_nonzero(matrix_cl_view::Upper, matrix_cl_view::Upper));
  EXPECT_EQ(true,
            contains_nonzero(matrix_cl_view::Entire, matrix_cl_view::Upper));
  EXPECT_EQ(true,
            contains_nonzero(matrix_cl_view::Entire, matrix_cl_view::Entire));
}

TEST(matrix_cl_view, transpose) {
  using stan::math::matrix_cl_view;
  using stan::math::transpose;
  EXPECT_EQ(matrix_cl_view::Lower, transpose(matrix_cl_view::Upper));
  EXPECT_EQ(matrix_cl_view::Upper, transpose(matrix_cl_view::Lower));
  EXPECT_EQ(matrix_cl_view::Diagonal, transpose(matrix_cl_view::Diagonal));
  EXPECT_EQ(matrix_cl_view::Entire, transpose(matrix_cl_view::Entire));
}

TEST(matrix_cl_view, invert) {
  using stan::math::invert;
  using stan::math::matrix_cl_view;
  EXPECT_EQ(matrix_cl_view::Lower, invert(matrix_cl_view::Upper));
  EXPECT_EQ(matrix_cl_view::Upper, invert(matrix_cl_view::Lower));
  EXPECT_EQ(matrix_cl_view::Entire, invert(matrix_cl_view::Diagonal));
  EXPECT_EQ(matrix_cl_view::Diagonal, invert(matrix_cl_view::Entire));
}

TEST(matrix_cl_view, from_eigen_uplo_type) {
  using stan::math::from_eigen_uplo_type;
  using stan::math::matrix_cl_view;
  EXPECT_EQ(matrix_cl_view::Lower, from_eigen_uplo_type(Eigen::Lower));
  EXPECT_EQ(matrix_cl_view::Upper, from_eigen_uplo_type(Eigen::Upper));
  EXPECT_EQ(matrix_cl_view::Entire, from_eigen_uplo_type(Eigen::SelfAdjoint));
  EXPECT_EQ(matrix_cl_view::Entire, from_eigen_uplo_type(Eigen::UnitDiag));
}

#endif
