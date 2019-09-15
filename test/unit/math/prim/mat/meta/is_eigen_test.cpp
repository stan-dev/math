#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <Eigen/Sparse>
#include <gtest/gtest.h>
#include <vector>

TEST(MathMeta, primitive) {
  using stan::is_eigen;
  EXPECT_FALSE((is_eigen<bool>::value));
  EXPECT_FALSE((is_eigen<double>::value));
  EXPECT_FALSE((is_eigen<int>::value));

  EXPECT_FALSE((is_eigen<std::vector<double>>::value));

  EXPECT_TRUE((is_eigen<Eigen::EigenBase<Eigen::MatrixXd>>::value));
  EXPECT_TRUE((is_eigen<Eigen::Matrix<double, -1, -1>>::value));
  EXPECT_TRUE((is_eigen<Eigen::SparseMatrix<double>>::value));
  EXPECT_TRUE((is_eigen<Eigen::MatrixBase<Eigen::MatrixXd>>::value));

  EXPECT_TRUE((is_eigen<const Eigen::Matrix<double, -1, -1>>::value));
  EXPECT_TRUE((is_eigen<Eigen::SparseMatrix<double>&>::value));
  EXPECT_TRUE(
      (is_eigen<Eigen::MatrixBase<Eigen::Matrix<double, -1, -1>>&&>::value));

  Eigen::Matrix<double, -1, -1> a;
  Eigen::Matrix<double, -1, -1> b;

  EXPECT_TRUE((is_eigen<decltype(a * b)>::value));
  EXPECT_TRUE((is_eigen<decltype(a * b + a.transpose())>::value));
}

TEST(MathMeta, expression) {
  using stan::is_eigen_matrix;
  EXPECT_FALSE((is_eigen_matrix<bool>::value));
  EXPECT_FALSE((is_eigen_matrix<double>::value));
  EXPECT_FALSE((is_eigen_matrix<int>::value));

  EXPECT_FALSE((is_eigen_matrix<std::vector<double>>::value));

  EXPECT_FALSE((is_eigen_matrix<Eigen::EigenBase<Eigen::MatrixXd>>::value));
  EXPECT_TRUE((is_eigen_matrix<Eigen::Matrix<double, -1, -1>>::value));
  EXPECT_TRUE((is_eigen_matrix<Eigen::SparseMatrix<double>>::value));
  EXPECT_FALSE((is_eigen_matrix<Eigen::MatrixBase<Eigen::MatrixXd>>::value));

  EXPECT_TRUE((is_eigen_matrix<const Eigen::Matrix<double, -1, -1>>::value));
  EXPECT_TRUE((is_eigen_matrix<Eigen::SparseMatrix<double>&>::value));
  EXPECT_FALSE(
      (is_eigen_matrix<Eigen::MatrixBase<Eigen::Matrix<double, -1, -1>>&&>::value));

  Eigen::Matrix<double, -1, -1> a;
  Eigen::Matrix<double, -1, -1> b;

  EXPECT_FALSE((is_eigen_matrix<decltype(a * b)>::value));
  EXPECT_FALSE((is_eigen_matrix<decltype(a * b + a.transpose())>::value));
}

#include <iostream>
#include <type_traits>
#include <stan/math/opencl/opencl.hpp>


template <typename Mat>
inline Eigen::MatrixXd foo(Mat&& A) {
  if (std::is_const<std::remove_reference_t<Mat>>::value) {
    std::cout << "You Shall Not Pass!" << "\n";
    // for example purposes they do in fact Pass
    return A;
  } else {
    stan::math::matrix_cl<double> mA = stan::math::matrix_cl<double>::constant(A);
    stan::math::matrix_cl<double> B(A.rows(), A.cols());
    B.zeros();
    const stan::math::matrix_cl<double> C = mA+B;
    return stan::math::from_matrix_cl(C);
  }
}

template <typename Mat>
inline Eigen::MatrixXd bad_boye(const Mat& A) {
  return foo(A);
}
TEST(MathMatrixCL, test_cache) {
  Eigen::MatrixXd A = Eigen::MatrixXd::Random(2,2);
  const Eigen::MatrixXd AA(2, 2);

  std::cout << A(0, 0) << std::endl;
  // pass non-const
  Eigen::MatrixXd B = bad_boye(A);
  // pass const
  Eigen::MatrixXd BB = bad_boye(AA);
  A(0,0) = 2;
  Eigen::MatrixXd C = bad_boye(A);
  std::cout << C(0, 0) << std::endl;
}
