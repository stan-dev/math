#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <limits>

class ErrorHandlingMatrixIsPosDef : public ::testing::Test {
 public:
  void SetUp() {}

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y;
};

TEST_F(ErrorHandlingMatrixIsPosDef, isPosDefinite) {
  using stan::math::is_pos_definite;
  y.resize(1, 1);
  y << 1;
  EXPECT_TRUE(is_pos_definite(y));

  Eigen::LLT<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > llt_1(y);
  EXPECT_TRUE(is_pos_definite(llt_1));

  Eigen::LDLT<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > ldlt_1
      = y.ldlt();
  EXPECT_TRUE(is_pos_definite(ldlt_1));

  y.resize(3, 3);
  y << 1, 0, 0, 0, 1, 0, 0, 0, 1;
  EXPECT_TRUE(is_pos_definite(y));

  Eigen::LLT<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > llt_2(y);
  EXPECT_TRUE(is_pos_definite(llt_2));

  Eigen::LDLT<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > ldlt_2
      = y.ldlt();
  EXPECT_TRUE(is_pos_definite(ldlt_2));
}

TEST_F(ErrorHandlingMatrixIsPosDef, isPosDefinite_not_square) {
  using stan::math::is_pos_definite;
  y.resize(3, 4);
  EXPECT_FALSE(is_pos_definite(y));
}

TEST_F(ErrorHandlingMatrixIsPosDef, isPosDefinite_0_size) {
  using stan::math::is_pos_definite;
  EXPECT_FALSE(is_pos_definite(y));

  Eigen::MatrixXd x;
  x.resize(0, 0);
  Eigen::LLT<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > llt(
      x.rows());
  llt.compute(x);
  EXPECT_TRUE(is_pos_definite(llt));
}

TEST_F(ErrorHandlingMatrixIsPosDef, isPosDefinite_non_symmetric) {
  using stan::math::is_pos_definite;
  y.resize(3, 3);
  y << 1, 0, 0, 0, 1, 0.5, 0, 0, 1;
  EXPECT_FALSE(is_pos_definite(y));

  Eigen::LLT<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > llt(y);
  EXPECT_TRUE(is_pos_definite(llt));

  Eigen::LDLT<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > ldlt
      = y.ldlt();
  EXPECT_TRUE(is_pos_definite(ldlt));
}

TEST_F(ErrorHandlingMatrixIsPosDef, isPosDefinite_non_pos_definite) {
  using stan::math::is_pos_definite;
  y.resize(3, 3);
  y << -1, 0, 0, 0, -1, 0, 0, 0, -1;
  EXPECT_FALSE(is_pos_definite(y));

  Eigen::LLT<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > llt_err1(
      y);
  EXPECT_FALSE(is_pos_definite(llt_err1));

  Eigen::LDLT<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > ldlt_err1
      = y.ldlt();
  EXPECT_FALSE(is_pos_definite(ldlt_err1));

  y.resize(2, 2);
  y << 1, 2, 2, 1;
  EXPECT_FALSE(is_pos_definite(y));

  Eigen::LLT<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > llt_err2(
      y);
  EXPECT_FALSE(is_pos_definite(llt_err2));

  Eigen::LDLT<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > ldlt_err2
      = y.ldlt();
  EXPECT_FALSE(is_pos_definite(ldlt_err2));

  y << 1, 1, 1, 1;
  EXPECT_FALSE(is_pos_definite(y));

  Eigen::LLT<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > llt_err3(
      y);
  EXPECT_FALSE(is_pos_definite(llt_err3));

  Eigen::LDLT<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > ldlt_err3
      = y.ldlt();
  EXPECT_FALSE(is_pos_definite(ldlt_err3));
}

TEST_F(ErrorHandlingMatrixIsPosDef, isPosDefinite_nan) {
  using stan::math::is_pos_definite;
  double nan = std::numeric_limits<double>::quiet_NaN();

  y.resize(1, 1);
  y << nan;
  EXPECT_FALSE(is_pos_definite(y));

  Eigen::LLT<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > llt_err1(
      y);
  EXPECT_FALSE(is_pos_definite(llt_err1));

  Eigen::LDLT<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > ldlt_err1
      = y.ldlt();
  EXPECT_FALSE(is_pos_definite(ldlt_err1));

  y.resize(3, 3);
  y << 2, -1, 0, -1, 2, -1, 0, -1, 2;
  EXPECT_TRUE(is_pos_definite(y));
  for (int i = 0; i < y.rows(); i++)
    for (int j = 0; j < y.cols(); j++) {
      y << 2, -1, 0, -1, 2, -1, 0, -1, 2;
      y(i, j) = nan;
      if (i >= j) {
        EXPECT_FALSE(is_pos_definite(y));
      }
    }

  y << 2, -1, nan, -1, 2, -1, nan, -1, nan;
  Eigen::LLT<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > llt_err2(
      y);
  EXPECT_FALSE(is_pos_definite(llt_err2));

  Eigen::LDLT<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > ldlt_err2
      = y.ldlt();
  EXPECT_FALSE(is_pos_definite(ldlt_err2));

  y << 0, 0, 0, 0, 0, 0, 0, 0, 0;
  EXPECT_FALSE(is_pos_definite(y));
}
