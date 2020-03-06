#ifndef TEST_UNIT_MATH_MIX_EIGEN_UTILS_HPP
#define TEST_UNIT_MATH_MIX_EIGEN_UTILS_HPP
#include <stan/math/mix.hpp>
#include <Eigen/Sparse>
#include <gtest/gtest.h>

namespace stan {
namespace math {
namespace test {

/**
 * Test type checking for sparse matrices
 * @tparam eigen_base_v logical for whether `Checker` should pass for
 *   `EigenBase` types.
 * @tparam eigen_sparse_compressed_v logical for whether `Checker` should pass
 * for `SparseCompressedBase` types.
 * @tparam eigen_sparse_matrix_v logical for whether `Checker` should pass for
 *   `SparseMatrixBase` types.
 * @tparam eigen_sparse_map_v logical for whether `Checker` should pass for
 *   `SparseMapBase` types.
 * @tparam Scalar The scalar type of the Eigen object to check.
 * @tparam Checker A struct taking in an expression template type and returning
 *  a true if check is satisfied and false otherwise.
 */
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
}

/**
 * Test type checking for dense matrices
 * @tparam eigen_base_v logical for whether `Checker` should pass for
 *   `EigenBase` types.
 * @tparam eigen_dense_v logical for whether `Checker` should pass for
 *   `DenseBase` types.
 * @tparam eigen_matrix_v logical for whether `Checker` should pass for
 *   `EigenMatrixBase` types.
 * @tparam eigen_matrix_expr_v logical for whether `Checker` should pass for
 *   `EigenMatrixBase` expression types.
 * @tparam eigen_map_v logical for whether `Checker` should pass for
 *   `Map` expression types.
 * @tparam Scalar The scalar type of the Eigen object to check.
 * @tparam Checker A struct taking in an expression template type and returning
 *  a true if check is satisfied and false otherwise.
 * @tparam EigenDims template parameter pack of `int` defining matrix
 *  compile time dimensions.
 */
template <
    bool eigen_base_v, bool eigen_dense_v, bool eigen_matrix_v,
    bool eigen_array_v, bool eigen_map_v, template <class...> class Checker,
    template <class, int...> class EigenType, typename Scalar, int... EigenDims>
void test_eigen_dense() {
  using Eigen::ArrayBase;
  using Eigen::DenseBase;
  using Eigen::EigenBase;
  using Eigen::Map;
  using Eigen::MatrixBase;
  using test_mat = EigenType<Scalar, EigenDims...>;
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
  EXPECT_TRUE((eigen_array_v == Checker<decltype(A * B)>::value));
  EXPECT_TRUE((eigen_map_v == Checker<Map<test_mat>>::value));
}

/*
 * Test type checking for Eigen dense solvers with matrics
 * @tparam eigen_dense_solver_v logical for whether `Checker` should pass for
 * `LDLT` and `LLT` Eigen types.
 * @tparam Scalar The scalar type of the Eigen type to test.
 * @tparam Checker A struct taking in a template parameter type and returning
 *  true when it's conditional is successful or false otherwise.
 */
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
}

/*
 * Test type checking for Eigen dense solvers with matrics
 * @tparam eigen_dense_solver_v logical for whether `Checker` should pass for
 * `LDLT` and `LLT` Eigen types.
 * @tparam Scalar The scalar type of the Eigen type to test.
 * @tparam Checker A struct taking in a template parameter type and returning
 *  true when it's conditional is successful or false otherwise.
 */
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
}

/**
 * Test type checking for dense matrices with all stan scalar types.
 * @tparam eigen_base_v logical for whether `Checker` should pass for
 *   `EigenBase` types.
 * @tparam eigen_dense_v logical for whether `Checker` should pass for
 *   `DenseBase` types.
 * @tparam eigen_matrix_v logical for whether `Checker` should pass for
 *   `EigenMatrixBase` types.
 * @tparam eigen_matrix_expr_v logical for whether `Checker` should pass for
 *   `EigenMatrixBase` expression types.
 * @tparam eigen_map_v logical for whether `Checker` should pass for
 *   `Map` expression types.
 * @tparam Checker A struct taking in an expression template type and returning
 *  a true if check is satisfied and false otherwise.
 * @tparam EigenDims template parameter pack of `int` defining matrix
 *  compile time dimensions.
 */
template <bool eigen_base_v, bool eigen_dense_v, bool eigen_matrix_v,
          bool eigen_matrix_expr_v, bool eigen_map_v,
          template <class...> class Checker, int... EigenDims>
void all_eigen_dense_matrix() {
  using stan::math::fvar;
  using stan::math::var;

  test_eigen_dense<eigen_base_v, eigen_dense_v, eigen_matrix_v,
                   eigen_matrix_expr_v, eigen_map_v, Checker, Eigen::Matrix,
                   double, EigenDims...>();
  test_eigen_dense<eigen_base_v, eigen_dense_v, eigen_matrix_v,
                   eigen_matrix_expr_v, eigen_map_v, Checker, Eigen::Matrix,
                   var, EigenDims...>();
  test_eigen_dense<eigen_base_v, eigen_dense_v, eigen_matrix_v,
                   eigen_matrix_expr_v, eigen_map_v, Checker, Eigen::Matrix,
                   fvar<double>, EigenDims...>();
  test_eigen_dense<eigen_base_v, eigen_dense_v, eigen_matrix_v,
                   eigen_matrix_expr_v, eigen_map_v, Checker, Eigen::Matrix,
                   fvar<var>, EigenDims...>();
  test_eigen_dense<eigen_base_v, eigen_dense_v, eigen_matrix_v,
                   eigen_matrix_expr_v, eigen_map_v, Checker, Eigen::Matrix,
                   fvar<fvar<double>>, EigenDims...>();
}

/**
 * Test type checking for dense arrays with all stan scalar types.
 * @tparam eigen_base_v logical for whether `Checker` should pass for
 *   `EigenBase` types.
 * @tparam eigen_dense_v logical for whether `Checker` should pass for
 *   `DenseBase` types.
 * @tparam eigen_matrix_v logical for whether `Checker` should pass for
 *   `EigenMatrixBase` types.
 * @tparam eigen_matrix_expr_v logical for whether `Checker` should pass for
 *   `EigenMatrixBase` expression types.
 * @tparam eigen_map_v logical for whether `Checker` should pass for
 *   `Map` expression types.
 * @tparam Checker A struct taking in an expression template type and returning
 *  a true if check is satisfied and false otherwise.
 * @tparam EigenDims template parameter pack of `int` defining matrix
 *  compile time dimensions.
 */
template <bool eigen_base_v, bool eigen_dense_v, bool eigen_array_v,
          bool eigen_array_expr_v, bool eigen_map_v,
          template <class...> class Checker, int... EigenDims>
void all_eigen_dense_array() {
  using stan::math::fvar;
  using stan::math::var;

  test_eigen_dense<eigen_base_v, eigen_dense_v, eigen_array_v,
                   eigen_array_expr_v, eigen_map_v, Checker, Eigen::Array,
                   double, EigenDims...>();

  test_eigen_dense<eigen_base_v, eigen_dense_v, eigen_array_v,
                   eigen_array_expr_v, eigen_map_v, Checker, Eigen::Array, var,
                   EigenDims...>();

  test_eigen_dense<eigen_base_v, eigen_dense_v, eigen_array_v,
                   eigen_array_expr_v, eigen_map_v, Checker, Eigen::Array,
                   fvar<double>, EigenDims...>();

  test_eigen_dense<eigen_base_v, eigen_dense_v, eigen_array_v,
                   eigen_array_expr_v, eigen_map_v, Checker, Eigen::Array,
                   fvar<var>, EigenDims...>();

  test_eigen_dense<eigen_base_v, eigen_dense_v, eigen_array_v,
                   eigen_array_expr_v, eigen_map_v, Checker, Eigen::Array,
                   fvar<fvar<double>>, EigenDims...>();
}

/*
 * Test type checking for Eigen dense solvers with matrics with all stan scalars
 * @tparam eigen_dense_solver_v logical for whether `Checker` should pass for
 * `LDLT` and `LLT` Eigen types.
 * @tparam Checker A struct taking in a template parameter type and returning
 *  true when it's conditional is successful or false otherwise.
 */
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

/**
 * Test type checking for sparse matrices with all stan scalar types.
 * @tparam eigen_base_v logical for whether `Checker` should pass for
 *   `EigenBase` types.
 * @tparam eigen_sparse_compressed_v logical for whether `Checker` should pass
 * for `SparseCompressedBase` types.
 * @tparam eigen_sparse_matrix_v logical for whether `Checker` should pass for
 *   `SparseMatrixBase` types.
 * @tparam eigen_sparse_map_v logical for whether `Checker` should pass for
 *   `SparseMapBase` types.
 * @tparam Checker A struct taking in an expression template type and returning
 *  a true if check is satisfied and false otherwise.
 */
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

}  // namespace test
}  // namespace math
}  // namespace stan

#endif
