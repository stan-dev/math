#include <stan/math/mix.hpp>
#include <Eigen/Sparse>
#include <gtest/gtest.h>


template <bool eigen_base_v, bool eigen_sparse_compressed_v, bool eigen_sparse_matrix_v, bool eigen_sparse_map_v, typename Scalar,
  template <class...> class Checker>
void test_eigen_sparse_matrix() {
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

template <bool eigen_base_v, bool eigen_dense_v, bool eigen_matrix_v, bool eigen_mat_exp_v, bool eigen_map_v,
  template <class...> class Checker, typename Scalar, int... EigenDims>
void test_eigen_dense_matrix() {
  using Eigen::EigenBase;
  using Eigen::DenseBase;
  using Eigen::MatrixBase;
  using Eigen::MapBase;
  using test_mat = Eigen::Matrix<Scalar, EigenDims...>;
  test_mat A;
  Eigen::Matrix<Scalar, -1, -1> B;
  // scalars
  EXPECT_FALSE((Checker<bool>::value));
  EXPECT_FALSE((Checker<double>::value));
  EXPECT_FALSE((Checker<int>::value));

  //  Dense Matrix Hierarchy
  EXPECT_TRUE((eigen_base_v == Checker<EigenBase<test_mat>>::value));
  EXPECT_TRUE((eigen_base_v == Checker<EigenBase<decltype(A * B)>>::value));
  EXPECT_TRUE((eigen_dense_v == Checker<DenseBase<decltype(A * B)>>::value));
  EXPECT_TRUE((eigen_matrix_v == Checker<MatrixBase<decltype(A * B)>>::value));
  EXPECT_TRUE((eigen_mat_exp_v == Checker<decltype(A * B)>::value));
  EXPECT_TRUE((eigen_map_v == Checker<MapBase<test_mat>>::value));
};


template <bool eigen_base_v, bool eigen_dense_v, bool eigen_matrix_v, bool eigen_map_v, bool eigen_mat_exp_v,
  template <class...> class Checker, typename Scalar, int... EigenDims>
void test_eigen_dense_array() {
  using Eigen::EigenBase;
  using Eigen::DenseBase;
  using Eigen::MatrixBase;
  using Eigen::MapBase;
  using test_mat = Eigen::Array<Scalar, EigenDims...>;
  test_mat A;
  Eigen::Matrix<Scalar, -1, -1> B;
  // scalars
  EXPECT_FALSE((Checker<bool>::value));
  EXPECT_FALSE((Checker<double>::value));
  EXPECT_FALSE((Checker<int>::value));

  //  Dense Matrix Hierarchy
  EXPECT_TRUE((eigen_base_v == Checker<EigenBase<test_mat>>::value));
  EXPECT_TRUE((eigen_base_v == Checker<EigenBase<decltype(A * B)>>::value));
  EXPECT_TRUE((eigen_dense_v == Checker<DenseBase<decltype(A * B)>>::value));
  EXPECT_TRUE((eigen_matrix_v == Checker<MatrixBase<decltype(A * B)>>::value));
  EXPECT_TRUE((eigen_map_v == Checker<MapBase<test_mat>>::value));
  EXPECT_TRUE((eigen_mat_exp_v == Checker<decltype(A * B)>::value));
};

template <bool eigen_dense_solver_v, typename Scalar, template <class...> class Checker>
void test_eigen_dense_decomp_matrix() {
  using Eigen::EigenBase;
  using Eigen::DenseBase;
  using Eigen::MatrixBase;
  using Eigen::MapBase;
  using test_mat = Eigen::Matrix<Scalar, -1, -1>;
  // Dense Solver
  EXPECT_TRUE((eigen_dense_solver_v == Checker<Eigen::LDLT<test_mat>>::value));
  EXPECT_TRUE((eigen_dense_solver_v == Checker<Eigen::LLT<test_mat>>::value));
};

template <bool eigen_dense_solver_v, typename Scalar, template <class...> class Checker>
void test_eigen_dense_decomp_array() {
  using Eigen::EigenBase;
  using Eigen::DenseBase;
  using Eigen::MatrixBase;
  using Eigen::MapBase;
  using test_arr = Eigen::Array<Scalar, -1, -1>;
  test_arr A(10, 10);
  // Dense Solver
  EXPECT_TRUE((eigen_dense_solver_v == Checker<Eigen::LDLT<decltype(A.matrix())>>::value));
  EXPECT_TRUE((eigen_dense_solver_v == Checker<Eigen::LLT<decltype(A.matrix())>>::value));
};


namespace stan {
  namespace test {
    namespace internal {
      template <typename T>
      using eigen_mat = Eigen::Matrix<T, -1, -1>;
      template <typename T>
      using eigen_arr = Eigen::Array<T, -1, -1>;
    }
  }
}

template <bool eigen_base_v, bool eigen_dense_v, bool eigen_matrix_v, bool eigen_map_v, bool eigen_mat_exp_v,
  template <class...> class Checker, int... EigenDims>
void test_all_eigen_dense_matrix() {
  using stan::math::var;
  using stan::math::fvar;
  using stan::test::internal::eigen_mat;
  using stan::test::internal::eigen_arr;

  test_eigen_dense_matrix<eigen_base_v, eigen_dense_v, eigen_matrix_v, eigen_map_v, eigen_mat_exp_v, Checker, double, EigenDims...>();
  test_eigen_dense_matrix<eigen_base_v, eigen_dense_v, eigen_matrix_v, eigen_map_v, eigen_mat_exp_v, Checker, var, EigenDims...>();
  test_eigen_dense_matrix<eigen_base_v, eigen_dense_v, eigen_matrix_v, eigen_map_v, eigen_mat_exp_v, Checker, fvar<double>, EigenDims...>();
  test_eigen_dense_matrix<eigen_base_v, eigen_dense_v, eigen_matrix_v, eigen_map_v, eigen_mat_exp_v, Checker, fvar<var>, EigenDims...>();
  test_eigen_dense_matrix<eigen_base_v, eigen_dense_v, eigen_matrix_v, eigen_map_v, eigen_mat_exp_v, Checker, fvar<fvar<double>>, EigenDims...>();

}

template <bool eigen_base_v, bool eigen_dense_v, bool eigen_matrix_v, bool eigen_map_v, bool eigen_mat_exp_v,
  template <class...> class Checker>
void test_all_eigen_dense_array() {
  using stan::math::var;
  using stan::math::fvar;
  using stan::test::internal::eigen_mat;
  using stan::test::internal::eigen_arr;

  test_eigen_dense_array<eigen_base_v, eigen_dense_v, eigen_matrix_v, eigen_map_v, eigen_mat_exp_v, double, Checker>();
  test_eigen_dense_array<eigen_base_v, eigen_dense_v, eigen_matrix_v, eigen_map_v, eigen_mat_exp_v, var, Checker>();
  test_eigen_dense_array<eigen_base_v, eigen_dense_v, eigen_matrix_v, eigen_map_v, eigen_mat_exp_v, fvar<double>, Checker>();
  test_eigen_dense_array<eigen_base_v, eigen_dense_v, eigen_matrix_v, eigen_map_v, eigen_mat_exp_v, fvar<var>, Checker>();
  test_eigen_dense_array<eigen_base_v, eigen_dense_v, eigen_matrix_v, eigen_map_v, eigen_mat_exp_v, fvar<fvar<double>>, Checker>();

}

template <bool eigen_dense_solver_v, template <class...> class Checker>
void test_all_eigen_dense_decomp() {
  using stan::math::var;
  using stan::math::fvar;

  test_eigen_dense_decomp_matrix<eigen_dense_solver_v, double, Checker>();
  test_eigen_dense_decomp_matrix<eigen_dense_solver_v, var, Checker>();
  test_eigen_dense_decomp_matrix<eigen_dense_solver_v, fvar<double>, Checker>();
  test_eigen_dense_decomp_matrix<eigen_dense_solver_v, fvar<var>, Checker>();
  test_eigen_dense_decomp_matrix<eigen_dense_solver_v, fvar<fvar<double>>, Checker>();

  test_eigen_dense_decomp_array<eigen_dense_solver_v, double, Checker>();
  test_eigen_dense_decomp_array<eigen_dense_solver_v, var, Checker>();
  test_eigen_dense_decomp_array<eigen_dense_solver_v, fvar<double>, Checker>();
  test_eigen_dense_decomp_array<eigen_dense_solver_v, fvar<var>, Checker>();
  test_eigen_dense_decomp_array<eigen_dense_solver_v, fvar<fvar<double>>, Checker>();
}

template <bool eigen_base_v, bool eigen_sparse_compressed_v, bool eigen_sparse_matrix_v, bool eigen_sparse_map_v,
  template <class...> class Checker>
void test_all_eigen_sparse() {
  using stan::math::var;
  using stan::math::fvar;

  test_eigen_sparse_matrix<eigen_base_v, eigen_sparse_compressed_v, eigen_sparse_matrix_v, eigen_sparse_map_v, double, Checker>();
  test_eigen_sparse_matrix<eigen_base_v, eigen_sparse_compressed_v, eigen_sparse_matrix_v, eigen_sparse_map_v, var, Checker>();
  test_eigen_sparse_matrix<eigen_base_v, eigen_sparse_compressed_v, eigen_sparse_matrix_v, eigen_sparse_map_v, fvar<double>, Checker>();
  test_eigen_sparse_matrix<eigen_base_v, eigen_sparse_compressed_v, eigen_sparse_matrix_v, eigen_sparse_map_v, fvar<var>, Checker>();
  test_eigen_sparse_matrix<eigen_base_v, eigen_sparse_compressed_v, eigen_sparse_matrix_v, eigen_sparse_map_v, fvar<fvar<var>>, Checker>();

}
