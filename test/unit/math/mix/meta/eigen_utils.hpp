#include <stan/math/prim.hpp>
#include <Eigen/Sparse>
#include <gtest/gtest.h>


template <bool eigen_base_v, bool eigen_sparse_compressed_v, bool eigen_sparse_matrix_v, bool eigen_sparse_map_v, typename Scalar,
  template <class...> class Checker>
void test_eigen_sparse() {
  using Eigen::EigenBase;
  using Eigen::SparseMatrix;
  using Eigen::SparseMatrixBase;
  using Eigen::SparseCompressedBase;
  using Eigen::SparseMapBase;
  using SparseMatrixXd = SparseMatrix<Scalar>;
  SparseMatrixXd C(10, 10);
  SparseMatrixXd D(10, 10);
  // scalars
  EXPECT_FALSE((Checker<bool>::value));
  EXPECT_FALSE((Checker<double>::value));
  EXPECT_FALSE((Checker<int>::value));

  // Sparse
  EXPECT_TRUE((eigen_base_v == Checker<EigenBase<decltype(C * D)>>::value));
  EXPECT_TRUE((eigen_sparse_compressed_v == Checker<SparseCompressedBase<decltype(C * D)>>::value));
  EXPECT_TRUE((eigen_sparse_matrix_v == Checker<SparseMatrixBase<decltype(C * D)>>::value));
  EXPECT_TRUE((eigen_sparse_matrix_v == Checker<decltype(C * D)>::value));
  EXPECT_TRUE((eigen_sparse_matrix_v == Checker<decltype(C * D)>::value));
  EXPECT_TRUE((eigen_sparse_map_v == Checker<SparseMapBase<decltype(C * D)>>::value));
};

template <bool eigen_base_v, bool eigen_dense_v, bool eigen_matrix_v, bool eigen_map_v, bool eigen_mat_exp_v, bool eigen_dense_solver_v, typename Scalar,
  template <class...> class Checker>
void test_eigen_dense() {
  using Eigen::EigenBase;
  using Eigen::DenseBase;
  using Eigen::MatrixBase;
  using Eigen::MapBase;
  using test_mat = Eigen::Matrix<Scalar, -1, -1>;
  test_mat A(10, 10);
  test_mat B(10, 10);
  // scalars
  EXPECT_FALSE((Checker<bool>::value));
  EXPECT_FALSE((Checker<double>::value));
  EXPECT_FALSE((Checker<int>::value));

  //  Dense Matrix Hierarchy
  EXPECT_TRUE((eigen_base_v == Checker<EigenBase<test_mat>>::value));
  EXPECT_TRUE((eigen_base_v == Checker<EigenBase<decltype(A * B)>>::value));
  EXPECT_TRUE((eigen_dense_v == Checker<DenseBase<decltype(A * B)>>::value));
  EXPECT_TRUE((eigen_matrix_v == Checker<MatrixBase<decltype(A * B)>>::value));
  EXPECT_TRUE((eigen_map_v == Checker<MapBase<decltype(A * B)>>::value));
  EXPECT_TRUE((eigen_mat_exp_v == Checker<decltype(A * B)>::value));
  // Dense Solver
  EXPECT_TRUE((eigen_dense_solver_v == Checker<Eigen::LDLT<test_mat>>::value));
  EXPECT_TRUE((eigen_dense_solver_v == Checker<Eigen::LLT<test_mat>>::value));

};
