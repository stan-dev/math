#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>

TEST(AgradRevMatrix, LDLT_alloc_default_constructor) {
  using stan::math::LDLT_alloc;
  using stan::math::var;

  // DO NOT DELETE, allocated on the vari stack
  LDLT_alloc<-1, -1> *alloc = new LDLT_alloc<-1, -1>();
  EXPECT_EQ(0U, alloc->N_);
#ifdef EIGEN_NO_DEBUG
  EXPECT_NO_THROW(alloc->log_abs_det());
  EXPECT_NO_THROW(alloc->ldlt_.info());
#else
  // Note: If -DEIGEN_NO_DEBUG is not included in the compilation flags
  //       asserts will force these calls to die instead of the above
  //       behavior

#ifndef _WIN32
  // Google test under Windows is having trouble with these tests.
  EXPECT_DEATH(alloc->log_abs_det(),
               "m_isInitialized && \"LDLT is not initialized.\"");
  EXPECT_DEATH(alloc->ldlt_.info(),
               "m_isInitialized && \"LDLT is not initialized.\"");
#endif

#endif
}

TEST(AgradRevMatrix, LDLT_alloc_constructor) {
  using stan::math::LDLT_alloc;
  using stan::math::var;

  Eigen::Matrix<var, -1, -1> A(2, 2);
  A << 2, 1, 1, 2;

  // DO NOT DELETE, allocated on the vari stack
  LDLT_alloc<-1, -1> *alloc = new LDLT_alloc<-1, -1>(A);

  EXPECT_EQ(2U, alloc->N_);
  EXPECT_FLOAT_EQ(1.0986122886681096, alloc->log_abs_det());
  EXPECT_EQ(Eigen::Success, alloc->ldlt_.info());

  Eigen::Matrix<double, -1, -1> expectedL(2, 2);
  expectedL(0, 0) = 1.0;
  expectedL(0, 1) = 0.0;
  expectedL(1, 0) = 0.5;
  expectedL(1, 1) = 1.0;

  Eigen::Matrix<double, -1, -1> L = alloc->ldlt_.matrixL();
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      EXPECT_FLOAT_EQ(expectedL(i, j), L(i, j));
}

TEST(AgradRevMatrix, LDLT_alloc_compute) {
  using stan::math::LDLT_alloc;
  using stan::math::var;

  // DO NOT DELETE, allocated on the vari stack
  LDLT_alloc<-1, -1> *alloc = new LDLT_alloc<-1, -1>();

  Eigen::Matrix<var, -1, -1> A(2, 2);
  A << 2, 1, 1, 2;

  EXPECT_NO_THROW(alloc->compute(A));
  EXPECT_EQ(2U, alloc->N_);
  EXPECT_EQ(Eigen::Success, alloc->ldlt_.info());
  EXPECT_FLOAT_EQ(alloc->log_abs_det(), 1.0986122886681096);

  Eigen::Matrix<double, -1, -1> expectedL(2, 2);
  expectedL(0, 0) = 1.0;
  expectedL(0, 1) = 0.0;
  expectedL(1, 0) = 0.5;
  expectedL(1, 1) = 1.0;

  Eigen::Matrix<double, -1, -1> L = alloc->ldlt_.matrixL();
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      EXPECT_FLOAT_EQ(expectedL(i, j), L(i, j));
}
