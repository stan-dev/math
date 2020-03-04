#ifndef TEST_UNIT_MATH_MIX_EIGEN_UTILS_HPP
#define TEST_UNIT_MATH_MIX_EIGEN_UTILS_HPP
#include <stan/math/mix.hpp>
#include <Eigen/Sparse>
#include <gtest/gtest.h>

namespace stan {
  namespace math {
  namespace test {

template <bool eigen_base_v, bool eigen_sparse_compressed_v,
          bool eigen_sparse_matrix_v, bool eigen_sparse_map_v, typename Scalar,
          template <class...> class Checker>
void test_eigen_sparse_matrix() {
  using Eigen::EigenBase;
  using Eigen::SparseCompressedBase;
  using Eigen::SparseMapBase;
  using Eigen::SparseMatrix;
  using Eigen::SparseMatrixBase;
  using SparseMatrixXd = SparseMatrix<Scalar>;
  SparseMatrixXd C(10, 10);
  SparseMatrixXd D(10, 10);
  // scalars
  EXPECT_FALSE((Checker<bool>::value));
  EXPECT_FALSE((Checker<double>::value));
  EXPECT_FALSE((Checker<int>::value));

  // Sparse
  EXPECT_TRUE((eigen_base_v == Checker<EigenBase<decltype(C * D)>>::value));
  EXPECT_TRUE((eigen_sparse_compressed_v
               == Checker<SparseCompressedBase<decltype(C * D)>>::value));
  EXPECT_TRUE((eigen_sparse_matrix_v
               == Checker<SparseMatrixBase<decltype(C * D)>>::value));
  EXPECT_TRUE((eigen_sparse_matrix_v == Checker<decltype(C * D)>::value));
  EXPECT_TRUE((eigen_sparse_matrix_v == Checker<decltype(C * D)>::value));
  EXPECT_TRUE(
      (eigen_sparse_map_v == Checker<SparseMapBase<decltype(C * D)>>::value));
};

template <bool eigen_base_v, bool eigen_dense_v, bool eigen_matrix_v,
          bool eigen_matrix_expr_v, bool eigen_map_v,
          template <class...> class Checker, typename Scalar, int... EigenDims>
void test_eigen_dense_matrix() {
  using Eigen::DenseBase;
  using Eigen::EigenBase;
  using Eigen::Map;
  using Eigen::MatrixBase;
  using test_mat = Eigen::Matrix<Scalar, EigenDims...>;
  test_mat A;
  test_mat B;
  // scalars
  EXPECT_FALSE((Checker<bool>::value));
  EXPECT_FALSE((Checker<double>::value));
  EXPECT_FALSE((Checker<int>::value));

  //  Dense Matrix Hierarchy
  EXPECT_TRUE((eigen_base_v == Checker<EigenBase<test_mat>>::value));
  EXPECT_TRUE((eigen_base_v == Checker<EigenBase<decltype(A * B)>>::value));
  EXPECT_TRUE((eigen_dense_v == Checker<DenseBase<decltype(A * B)>>::value));
  EXPECT_TRUE((eigen_matrix_v == Checker<MatrixBase<decltype(A * B)>>::value));
  EXPECT_TRUE((eigen_matrix_expr_v == Checker<decltype(A * B)>::value));
  EXPECT_TRUE((eigen_map_v == Checker<Map<test_mat>>::value));
};

template <bool eigen_base_v, bool eigen_dense_v, bool eigen_array_v,
          bool eigen_matrix_expr_v, bool eigen_map_v,
          template <class...> class Checker, typename Scalar, int... EigenDims>
void test_eigen_dense_array() {
  using Eigen::DenseBase;
  using Eigen::EigenBase;
  using Eigen::Map;
  using Eigen::ArrayBase;
  using test_mat = Eigen::Array<Scalar, EigenDims...>;
  test_mat A;
  test_mat B;
  // scalars
  EXPECT_FALSE((Checker<bool>::value));
  EXPECT_FALSE((Checker<double>::value));
  EXPECT_FALSE((Checker<int>::value));

  //  Dense array Hierarchy
  EXPECT_TRUE((eigen_base_v == Checker<EigenBase<test_mat>>::value));
  EXPECT_TRUE((eigen_base_v == Checker<EigenBase<decltype(A * B)>>::value));
  EXPECT_TRUE((eigen_dense_v == Checker<DenseBase<decltype(A * B)>>::value));
  EXPECT_TRUE((eigen_array_v == Checker<ArrayBase<decltype(A * B)>>::value));
  // TODO(Steve): Figure out why uncommenting this causes compiler error
  //EXPECT_TRUE((eigen_matrix_expr_v == Checker<decltype(A * B)>::value));
  EXPECT_TRUE((eigen_map_v == Checker<Map<test_mat>>::value));
};

template <bool eigen_dense_solver_v, typename Scalar,
          template <class...> class Checker>
void test_eigen_dense_decomp_matrix() {
  using Eigen::DenseBase;
  using Eigen::EigenBase;
  using Eigen::MapBase;
  using Eigen::MatrixBase;
  using test_mat = Eigen::Matrix<Scalar, -1, -1>;
  // Dense Solver
  EXPECT_TRUE((eigen_dense_solver_v == Checker<Eigen::LDLT<test_mat>>::value));
  EXPECT_TRUE((eigen_dense_solver_v == Checker<Eigen::LLT<test_mat>>::value));
};

template <bool eigen_dense_solver_v, typename Scalar,
          template <class...> class Checker>
void test_eigen_dense_decomp_array() {
  using Eigen::DenseBase;
  using Eigen::EigenBase;
  using Eigen::MapBase;
  using Eigen::MatrixBase;
  using test_arr = Eigen::Array<Scalar, -1, -1>;
  test_arr A(10, 10);
  // Dense Solver
  EXPECT_TRUE((eigen_dense_solver_v
               == Checker<Eigen::LDLT<decltype(A.matrix())>>::value));
  EXPECT_TRUE((eigen_dense_solver_v
               == Checker<Eigen::LLT<decltype(A.matrix())>>::value));
};


template <typename T>
using eigen_mat = Eigen::Matrix<T, -1, -1>;
template <typename T>
using eigen_arr = Eigen::Array<T, -1, -1>;

template <bool eigen_base_v, bool eigen_dense_v, bool eigen_matrix_v,
          bool eigen_matrix_expr_v, bool eigen_map_v,
          template <class...> class Checker, int... EigenDims>
void all_eigen_dense_matrix() {
  using stan::math::fvar;
  using stan::math::var;

  test_eigen_dense_matrix<eigen_base_v, eigen_dense_v, eigen_matrix_v,
                          eigen_matrix_expr_v, eigen_map_v, Checker, double,
                          EigenDims...>();
/*  test_eigen_dense_matrix<eigen_base_v, eigen_dense_v, eigen_matrix_v,
                          eigen_map_v, eigen_matrix_expr_v, Checker, var,
                          EigenDims...>();
  test_eigen_dense_matrix<eigen_base_v, eigen_dense_v, eigen_matrix_v,
                          eigen_map_v, eigen_matrix_expr_v, Checker, fvar<double>,
                          EigenDims...>();
  test_eigen_dense_matrix<eigen_base_v, eigen_dense_v, eigen_matrix_v,
                          eigen_map_v, eigen_matrix_expr_v, Checker, fvar<var>,
                          EigenDims...>();
  test_eigen_dense_matrix<eigen_base_v, eigen_dense_v, eigen_matrix_v,
                          eigen_map_v, eigen_matrix_expr_v, Checker,
                          fvar<fvar<double>>, EigenDims...>();
*/}

template <bool eigen_base_v, bool eigen_dense_v, bool eigen_array_v,
          bool eigen_array_expr_v, bool eigen_map_v,
          template <class...> class Checker, int... EigenDims>
void all_eigen_dense_array() {
  using stan::math::fvar;
  using stan::math::var;

  test_eigen_dense_array<eigen_base_v, eigen_dense_v, eigen_array_v,
                         eigen_array_expr_v, eigen_map_v, Checker, double, EigenDims...>();
/*  test_eigen_dense_array<eigen_base_v, eigen_dense_v, eigen_array_v,
                         eigen_map_v, eigen_array_expr_v, Checker, var, EigenDims...>();
  test_eigen_dense_array<eigen_base_v, eigen_dense_v, eigen_array_v,
                         eigen_map_v, eigen_array_expr_v, Checker, fvar<double>, EigenDims...>();
  test_eigen_dense_array<eigen_base_v, eigen_dense_v, eigen_array_v,
                         eigen_map_v, eigen_array_expr_v, Checker, fvar<var>, EigenDims...>();
  test_eigen_dense_array<eigen_base_v, eigen_dense_v, eigen_array_v,
                         eigen_map_v, eigen_array_expr_v, Checker, fvar<fvar<double>>, EigenDims...>();
*/}

template <bool eigen_dense_solver_v, template <class...> class Checker>
void all_eigen_dense_decomp() {
  using stan::math::fvar;
  using stan::math::var;

  test_eigen_dense_decomp_matrix<eigen_dense_solver_v, double, Checker>();
  test_eigen_dense_decomp_matrix<eigen_dense_solver_v, var, Checker>();
  test_eigen_dense_decomp_matrix<eigen_dense_solver_v, fvar<double>, Checker>();
  test_eigen_dense_decomp_matrix<eigen_dense_solver_v, fvar<var>, Checker>();
  test_eigen_dense_decomp_matrix<eigen_dense_solver_v, fvar<fvar<double>>,
                                 Checker>();

  test_eigen_dense_decomp_array<eigen_dense_solver_v, double, Checker>();
  test_eigen_dense_decomp_array<eigen_dense_solver_v, var, Checker>();
  test_eigen_dense_decomp_array<eigen_dense_solver_v, fvar<double>, Checker>();
  test_eigen_dense_decomp_array<eigen_dense_solver_v, fvar<var>, Checker>();
  test_eigen_dense_decomp_array<eigen_dense_solver_v, fvar<fvar<double>>,
                                Checker>();
}

template <bool eigen_base_v, bool eigen_sparse_compressed_v,
          bool eigen_sparse_matrix_v, bool eigen_sparse_map_v,
          template <class...> class Checker>
void all_eigen_sparse() {
  using stan::math::fvar;
  using stan::math::var;

  test_eigen_sparse_matrix<eigen_base_v, eigen_sparse_compressed_v,
                           eigen_sparse_matrix_v, eigen_sparse_map_v, double,
                           Checker>();
  test_eigen_sparse_matrix<eigen_base_v, eigen_sparse_compressed_v,
                           eigen_sparse_matrix_v, eigen_sparse_map_v, var,
                           Checker>();
  test_eigen_sparse_matrix<eigen_base_v, eigen_sparse_compressed_v,
                           eigen_sparse_matrix_v, eigen_sparse_map_v,
                           fvar<double>, Checker>();
  test_eigen_sparse_matrix<eigen_base_v, eigen_sparse_compressed_v,
                           eigen_sparse_matrix_v, eigen_sparse_map_v, fvar<var>,
                           Checker>();
  test_eigen_sparse_matrix<eigen_base_v, eigen_sparse_compressed_v,
                           eigen_sparse_matrix_v, eigen_sparse_map_v,
                           fvar<fvar<var>>, Checker>();
}

}
}
}

#endif
