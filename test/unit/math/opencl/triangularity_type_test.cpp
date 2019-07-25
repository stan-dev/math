#ifdef STAN_OPENCL
#include <stan/math/opencl/partial_types.hpp>
#include <gtest/gtest.h>

TEST(PartialViewCL, operator_add) {
  using stan::math::PartialViewCL;
  EXPECT_EQ(PartialViewCL::Lower, PartialViewCL::Lower + PartialViewCL::Lower);
  EXPECT_EQ(PartialViewCL::Entire, PartialViewCL::Lower + PartialViewCL::Upper);
  EXPECT_EQ(PartialViewCL::Lower,
            PartialViewCL::Lower + PartialViewCL::Diagonal);

  EXPECT_EQ(PartialViewCL::Upper, PartialViewCL::Upper + PartialViewCL::Upper);
  EXPECT_EQ(PartialViewCL::Upper,
            PartialViewCL::Upper + PartialViewCL::Diagonal);

  EXPECT_EQ(PartialViewCL::Entire,
            PartialViewCL::Entire + PartialViewCL::Upper);
  EXPECT_EQ(PartialViewCL::Entire,
            PartialViewCL::Entire + PartialViewCL::Lower);
  EXPECT_EQ(PartialViewCL::Entire,
            PartialViewCL::Entire + PartialViewCL::Diagonal);
  EXPECT_EQ(PartialViewCL::Entire,
            PartialViewCL::Entire + PartialViewCL::Entire);
}

TEST(PartialViewCL, operator_multiply) {
  using stan::math::PartialViewCL;
  EXPECT_EQ(PartialViewCL::Diagonal,
            PartialViewCL::Upper * PartialViewCL::Diagonal);
  EXPECT_EQ(PartialViewCL::Diagonal,
            PartialViewCL::Entire * PartialViewCL::Diagonal);
  EXPECT_EQ(PartialViewCL::Diagonal,
            PartialViewCL::Lower * PartialViewCL::Diagonal);
  EXPECT_EQ(PartialViewCL::Diagonal,
            PartialViewCL::Lower * PartialViewCL::Upper);

  EXPECT_EQ(PartialViewCL::Lower, PartialViewCL::Lower * PartialViewCL::Lower);
  EXPECT_EQ(PartialViewCL::Lower, PartialViewCL::Entire * PartialViewCL::Lower);

  EXPECT_EQ(PartialViewCL::Upper, PartialViewCL::Upper * PartialViewCL::Upper);
  EXPECT_EQ(PartialViewCL::Upper, PartialViewCL::Entire * PartialViewCL::Upper);
  EXPECT_EQ(PartialViewCL::Entire,
            PartialViewCL::Entire * PartialViewCL::Entire);
}

TEST(PartialViewCL, is_not_diagonal) {
  using stan::math::PartialViewCL;
  using stan::math::is_not_diagonal;
  EXPECT_EQ(false,
            is_not_diagonal(PartialViewCL::Upper, PartialViewCL::Diagonal));
  EXPECT_EQ(false,
            is_not_diagonal(PartialViewCL::Entire, PartialViewCL::Diagonal));
  EXPECT_EQ(false,
            is_not_diagonal(PartialViewCL::Lower, PartialViewCL::Diagonal));
  EXPECT_EQ(false, is_not_diagonal(PartialViewCL::Lower, PartialViewCL::Upper));

  EXPECT_EQ(true, is_not_diagonal(PartialViewCL::Lower, PartialViewCL::Lower));
  EXPECT_EQ(true, is_not_diagonal(PartialViewCL::Entire, PartialViewCL::Lower));

  EXPECT_EQ(true, is_not_diagonal(PartialViewCL::Upper, PartialViewCL::Upper));
  EXPECT_EQ(true, is_not_diagonal(PartialViewCL::Entire, PartialViewCL::Upper));
  EXPECT_EQ(true,
            is_not_diagonal(PartialViewCL::Entire, PartialViewCL::Entire));
}

TEST(PartialViewCL, transpose) {
  using stan::math::PartialViewCL;
  using stan::math::transpose;
  EXPECT_EQ(PartialViewCL::Lower, transpose(PartialViewCL::Upper));
  EXPECT_EQ(PartialViewCL::Upper, transpose(PartialViewCL::Lower));
  EXPECT_EQ(PartialViewCL::Diagonal, transpose(PartialViewCL::Diagonal));
  EXPECT_EQ(PartialViewCL::Entire, transpose(PartialViewCL::Entire));
}

TEST(PartialViewCL, invert) {
  using stan::math::PartialViewCL;
  using stan::math::invert;
  EXPECT_EQ(PartialViewCL::Lower, invert(PartialViewCL::Upper));
  EXPECT_EQ(PartialViewCL::Upper, invert(PartialViewCL::Lower));
  EXPECT_EQ(PartialViewCL::Entire, invert(PartialViewCL::Diagonal));
  EXPECT_EQ(PartialViewCL::Diagonal, invert(PartialViewCL::Entire));
}

TEST(PartialViewCL, from_eigen_triangular_type) {
  using stan::math::PartialViewCL;
  using stan::math::from_eigen_triangular_type;
  EXPECT_EQ(PartialViewCL::Lower, from_eigen_triangular_type(Eigen::Lower));
  EXPECT_EQ(PartialViewCL::Upper, from_eigen_triangular_type(Eigen::Upper));
  EXPECT_EQ(PartialViewCL::Entire,
            from_eigen_triangular_type(Eigen::SelfAdjoint));
  EXPECT_EQ(PartialViewCL::Entire, from_eigen_triangular_type(Eigen::UnitDiag));
}

#endif
