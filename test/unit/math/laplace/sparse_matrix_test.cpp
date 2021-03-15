#include <stan/math/laplace/block_matrix_sqrt.hpp>

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/LU>

#include <gtest/gtest.h>
#include <iostream>
#include <istream>
#include <fstream>
#include <vector>


TEST(sparse_matrix, eigen_example) {
  typedef Eigen::Triplet<double> trp;
  using Eigen::SparseMatrix;
  using Eigen::VectorXi;
  using Eigen::MatrixXd;

  int m = 2;  // size of each block

  std::vector<trp> triplet_list(8);
  triplet_list[0] = trp(0, 0, 4);
  triplet_list[1] = trp(0, 1, 3);
  triplet_list[2] = trp(1, 0, 3);
  triplet_list[3] = trp(1, 1, 6);
  triplet_list[4] = trp(2, 2, 4);
  triplet_list[5] = trp(2, 3, 2);
  triplet_list[6] = trp(3, 2, 7);
  triplet_list[7] = trp(3, 3, 8);

  SparseMatrix<double> A(4, 4);
  A.setFromTriplets(triplet_list.begin(), triplet_list.end());

  std::cout << "A: " << A << std::endl;

  // Alternatively, can construct a matrix without using triplets.
  SparseMatrix<double> B(4, 4);
  B.reserve(VectorXi::Constant(B.cols(), 2));
  B.insert(0, 0) = 1;
  B.insert(0, 1) = 3;
  B.insert(1, 0) = 5;
  B.insert(1, 1) = 6;
  B.insert(2, 2) = 4;
  B.insert(2, 3) = 2;
  B.insert(3, 2) = 7;
  B.insert(3, 3) = 8;
  B.makeCompressed();

  std::cout << "B: " << B << std::endl;

  // If storage order mathces, we can sparse matrices.
  std::cout << "A + B: " << A + B << std::endl;

  std::cout << "A * B: " << A * B << std::endl;

  MatrixXd C(4, 4);
  C << 1, 3, 4, 5,
       3, 4, 4, 1,
       8, 1, 0, 12,
       3, 4, 5, 1;

  std::cout << "A * C: " << A * C << std::endl;

  SparseMatrix<double> sqrt_A = stan::math::block_matrix_sqrt(A, 2);

  std::cout << "sqrt(A): " << sqrt_A << std::endl;

  std::cout << "Check we recover A: " << sqrt_A * sqrt_A << std::endl;


  SparseMatrix<double> D(4, 4);

  std::cout << "sqrt(D): " << stan::math::block_matrix_sqrt(D, 2) << std::endl;

  MatrixXd E = MatrixXd::Zero(4, 4);
  std::cout << "sqrt(E): " << E.sqrt() << std::endl;
}

TEST(LU_decomposition, eigen_example) {
  using namespace Eigen;
  // typedef Matrix<double, 5, 3> Matrix5x3;
  typedef Matrix<double, 5, 5> Matrix5x5;
  using std::cout;
  using std::endl;

  Matrix5x5 m = Matrix5x5::Random();
  cout << "Here is the matrix m:" << endl << m << endl;
  Eigen::FullPivLU<Matrix5x5> lu(m);
  cout << "Here is, up to permutations, its LU decomposition matrix:"
       << endl << lu.matrixLU() << endl;
  cout << "Here is the L part:" << endl;
  // Matrix5x5 l = Matrix5x5::Identity();
  // l.block<5,3>(0,0).triangularView<StrictlyLower>() = lu.matrixLU();
  Eigen::MatrixXd l(5, 5);
  // l.block<5, 3>(0, 0).triangularView<StrictlyLower>() = lu.matrixLU();
  l.triangularView<Lower>() = lu.matrixLU();
  cout << l << endl;
  cout << "Here is the U part:" << endl;
  Matrix5x5 u = lu.matrixLU().triangularView<Upper>();
  cout << u << endl;
  // cout << "Here are the two triangular matrices we would store: " << endl;
  // cout << lu.permutationP().inverse() * l << endl;
  // cout << u * lu.permutationQ().inverse() << endl;
  cout << "Let us now reconstruct the original matrix m:" << endl;
  cout << lu.permutationP().inverse() * l * u * lu.permutationQ().inverse()
       << endl;

  Eigen::MatrixXd L(5, 5);
  L.triangularView<StrictlyLower>() = lu.permutationP().inverse() * l;
  std::cout << "L: " << L << std::endl;


  Eigen::MatrixXd U(5, 5);
  U = u * lu.permutationQ().inverse();
  std::cout << "U: " << U << std::endl;

  cout << "Determinant through decomposition: "
       << l.determinant() * u.determinant() << std::endl
       << "Determinant through direct comp: " << m.determinant() << std::endl;
}

TEST(LU_decomposition, eigen_example_2) {
  using Eigen::MatrixXd;

  MatrixXd A = MatrixXd::Random(3, 3);
  MatrixXd B = MatrixXd::Random(3, 2);
  std::cout << "Here is matrix A:" << std::endl << A << std::endl;
  std::cout << "Here us matrix B:" << std::endl << B << std::endl;

  Eigen::PartialPivLU<MatrixXd> LU;
  LU = Eigen::PartialPivLU<MatrixXd>(A);
  std::cout << "LU determinant: " << LU.determinant() << std::endl;
  std::cout << "A determinant: " << A.determinant() << std::endl;

  std::cout << "A.solve(B): " << std::endl
            << LU.solve(B) << std::endl;

  std::cout << "Check solution: " << std::endl
            << A * LU.solve(B) << std::endl;
}
