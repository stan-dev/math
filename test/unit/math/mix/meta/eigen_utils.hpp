#ifndef TEST_UNIT_MATH_MIX_EIGEN_UTILS_HPP
#define TEST_UNIT_MATH_MIX_EIGEN_UTILS_HPP
#include <stan/math/mix.hpp>
#include <test/unit/pretty_print_types.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <gtest/gtest.h>
#include <vector>

namespace stan {
namespace math {
namespace test {

/**
 * Test type checking for sparse matrices
 * @tparam Base_v logical for whether `Checker` should pass for
 *   `EigenBase` types.
 * @tparam SparseCompressed_v logical for whether `Checker` should pass
 * for `SparseCompressedBase` types.
 * @tparam SparseMatrix_v logical for whether `Checker` should pass for
 *   `SparseMatrixBase` types.
 * @tparam SparseMap_v logical for whether `Checker` should pass for
 *   `SparseMapBase` types.
 * @tparam Scalar The scalar type of the Eigen object to check.
 * @tparam Checker A struct taking in an expression template type and returning
 *  a true if check is satisfied and false otherwise.
 */
template <bool Base_v, bool SparseCompressed_v, bool SparseMatrix_v,
          bool SparseMap_v, typename Scalar, template <class...> class Checker>
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
  EXPECT_TRUE((Base_v == Checker<EigenBase<decltype(C * D)>>::value));
  EXPECT_TRUE((SparseCompressed_v
               == Checker<SparseCompressedBase<decltype(C * D)>>::value));
  EXPECT_TRUE(
      (SparseMatrix_v == Checker<SparseMatrixBase<decltype(C * D)>>::value));
  EXPECT_TRUE((SparseMatrix_v == Checker<decltype(C * D)>::value));
  EXPECT_TRUE((SparseMatrix_v == Checker<decltype(C + D)>::value));
  EXPECT_TRUE((SparseMap_v == Checker<SparseMapBase<decltype(C * D)>>::value));
}

/**
 * Test type checking for dense matrices
 * @tparam Base_v logical for whether `Checker` should pass for
 *   `EigenBase` types.
 * @tparam Dense_v logical for whether `Checker` should pass for
 *   `DenseBase` types.
 * @tparam Matrix_v logical for whether `Checker` should pass for
 *   `MatrixBase` types.
 * @tparam Array_v logical for whether `Checker` should pass for
 *   `ArrayBase` types.
 * @tparam Map_v logical for whether `Checker` should pass for
 *   `Map` expression types.
 * @tparam Checker A struct taking in an expression template type and returning
 *  a true if check is satisfied and false otherwise.
 * @tparam EigenType The type from the Eigen namespace to validate type traits
 * with.
 */
template <bool Base_v, bool Dense_v, bool Matrix_v, bool Array_v, bool Map_v,
          template <class...> class Checker, typename EigenType>
void test_eigen_dense_hierarchy() {
  using Eigen::ArrayBase;
  using Eigen::DenseBase;
  using Eigen::EigenBase;
  using Eigen::Map;
  using Eigen::MatrixBase;
  std::remove_reference_t<EigenType> A;
  std::remove_reference_t<EigenType> B;

  //  Dense Matrix Hierarchy
  EXPECT_TRUE((Base_v == Checker<EigenBase<EigenType>>::value))
      << "Failed For Base: " << type_name<EigenType>()
      << "\nChecking Type: " << type_name<EigenBase<EigenType>>();
  EXPECT_TRUE((Base_v == Checker<EigenBase<decltype(A * B)>>::value))
      << "Failed For Base: " << type_name<EigenType>()
      << "\nChecking Type: " << type_name<EigenBase<decltype(A * B)>>();
  EXPECT_TRUE((Dense_v == Checker<DenseBase<decltype(A * B)>>::value))
      << "Failed For Base: " << type_name<EigenType>()
      << "\nChecking Type: " << type_name<DenseBase<decltype(A * B)>>();
  EXPECT_TRUE(
      (Matrix_v == Checker<MatrixBase<decltype(A * B + A.transpose())>>::value))
      << "Failed For Base: " << type_name<EigenType>()
      << "\nChecking Type: " << type_name<decltype(A * B + A.transpose())>();
  EXPECT_TRUE((Array_v == Checker<ArrayBase<decltype(A + B)>>::value))
      << "Failed For Base: " << type_name<EigenType>()
      << "\nChecking Type: " << type_name<ArrayBase<decltype(A + B)>>();
  EXPECT_TRUE((Map_v == Checker<Map<EigenType>>::value))
      << "Failed For Base: " << type_name<EigenType>()
      << "\nChecking Type: " << type_name<Map<EigenType>>();
}

/**
 * Check if Eigen type traits satsify more complex Eigen expressions.
 * @tparam Base_v `bool` for the type of `EigenType` should pass.
 * @tparam Expr_v `bool` for if a product addition and transpose expression
 * passes of `EigenType` should pass.
 * @tparam Segment_v `bool` for if calling `segment()` member from `EigenType`
 * should pass.
 * @tparam Block_v `bool` for if calling `block` member from `EigenType` should
 * pass.
 * @tparam Checker A type trait returning a `bool` value.
 * @tparam EigenType A type from the Eigen namespace to check against.
 *  Must have addition, multiplication, `segment(int, int)`, and
 * `block(int, int, int, int)` members.
 */
template <bool Base_v, bool Expr_v, bool Segment_v, bool Block_v,
          template <class...> class Checker, typename EigenType>
void test_eigen_dense_exprs() {
  std::remove_reference_t<EigenType> A;
  std::remove_reference_t<EigenType> B;

  //  Dense Ops
  EXPECT_TRUE((Base_v == Checker<decltype(A)&>::value))
      << "Failed For Base: " << type_name<EigenType&>()
      << "\nChecking Type: " << type_name<EigenType>();
  EXPECT_TRUE((Base_v == Checker<decltype(A)&&>::value))
      << "Failed For Base: " << type_name<EigenType&&>()
      << "\nChecking Type: " << type_name<EigenType>();
  EXPECT_TRUE((Base_v == Checker<const EigenType&>::value))
      << "Failed For Base: " << type_name<const EigenType&>()
      << "\nChecking Type: " << type_name<EigenType>();
  EXPECT_TRUE((Base_v == Checker<EigenType>::value))
      << "Failed For Base: " << type_name<EigenType>()
      << "\nChecking Type: " << type_name<EigenType>();
  EXPECT_TRUE((Expr_v == Checker<decltype(A * B + A.transpose())>::value))
      << "Failed For Base: " << type_name<EigenType>()
      << "\nChecking Type: " << type_name<decltype(A * B + A.transpose())>();
  EXPECT_TRUE((Segment_v == Checker<decltype(A.segment(0, 1))>::value))
      << "Failed For Base: " << type_name<EigenType>()
      << "\nChecking Type: " << type_name<decltype(A.segment(0, 1))>();
  EXPECT_TRUE((Block_v == Checker<decltype(A.block(0, 0, 1, 1))>::value))
      << "Failed For Base: " << type_name<EigenType>()
      << "\nChecking Type: " << type_name<decltype(A.block(0, 0, 1, 1))>();
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
 * Test type checking for dense arrays with all stan scalar types.
 * @tparam Base_v logical for whether `Checker` should pass for
 *   `EigenBase` types.
 * @tparam Dense_v logical for whether `Checker` should pass for
 *   `DenseBase` types.
 * @tparam Matrix_v logical for whether `Checker` should pass for
 *   `MatrixBase` types.
 * @tparam Array_v logical for whether `Checker` should pass for
 *   `ArrayBase` types.
 * @tparam Map_v logical for whether `Checker` should pass for
 *   `Map` expression types.
 * @tparam Checker A struct taking in an expression template type and returning
 *  a true if check is satisfied and false otherwise.
 * @tparam EigenType A template type from Eigen to check type traits against.
 * @tparam EigenDims template parameter pack of `int` defining `EigenType`
 * compile time dimensions. compile time dimensions.
 */
template <bool Base_v, bool Dense_v, bool Matrix_v, bool Array_v, bool Map_v,
          template <class...> class Checker,
          template <class Scalar, int...> class EigenType, int... EigenDims>
void all_eigen_dense() {
  using stan::math::fvar;
  using stan::math::var;

  // scalars
  EXPECT_FALSE((Checker<bool>::value));
  EXPECT_FALSE((Checker<double>::value));
  EXPECT_FALSE((Checker<int>::value));
  EXPECT_FALSE((Checker<std::vector<double>>::value));

  test_eigen_dense_hierarchy<Base_v, Dense_v, Matrix_v, Array_v, Map_v, Checker,
                             EigenType<double, EigenDims...>>();
  test_eigen_dense_hierarchy<Base_v, Dense_v, Matrix_v, Array_v, Map_v, Checker,
                             EigenType<var, EigenDims...>>();
  test_eigen_dense_hierarchy<Base_v, Dense_v, Matrix_v, Array_v, Map_v, Checker,
                             EigenType<fvar<double>, EigenDims...>>();
  test_eigen_dense_hierarchy<Base_v, Dense_v, Matrix_v, Array_v, Map_v, Checker,
                             EigenType<fvar<var>, EigenDims...>>();
  test_eigen_dense_hierarchy<Base_v, Dense_v, Matrix_v, Array_v, Map_v, Checker,
                             EigenType<fvar<fvar<double>>, EigenDims...>>();
}

/**
 * Check if Eigen type traits satsify more complex Eigen expressions.
 * @tparam Base_v `bool` for the type of `EigenType` should pass.
 * @tparam Expr_v `bool` for if a product addition and transpose expression
 * passes of `EigenType` should pass.
 * @tparam Segment_v `bool` for if calling `segment()` member from `EigenType`
 * should pass.
 * @tparam Block_v `bool` for if calling `block` member from `EigenType` should
 * pass.
 * @tparam Checker A type trait returning a `bool` value.
 * @tparam EigenType A type from the Eigen namespace to check against.
 *  Must have addition, multiplication, `segment(int, int)`, and
 * `block(int, int, int, int)` members.
 */
template <bool Base_v, bool Expr_v, bool Segment_v, bool Block_v,
          template <class...> class Checker,
          template <class Scalar, int...> class EigenType, int... EigenDims>
void all_eigen_dense_exprs() {
  using stan::math::fvar;
  using stan::math::var;
  test_eigen_dense_exprs<Base_v, Expr_v, Segment_v, Block_v, Checker,
                         EigenType<double, EigenDims...>>();
  test_eigen_dense_exprs<Base_v, Expr_v, Segment_v, Block_v, Checker,
                         EigenType<var, EigenDims...>>();
  test_eigen_dense_exprs<Base_v, Expr_v, Segment_v, Block_v, Checker,
                         EigenType<fvar<double>, EigenDims...>>();
  test_eigen_dense_exprs<Base_v, Expr_v, Segment_v, Block_v, Checker,
                         EigenType<fvar<var>, EigenDims...>>();
  test_eigen_dense_exprs<Base_v, Expr_v, Segment_v, Block_v, Checker,
                         EigenType<fvar<fvar<double>>, EigenDims...>>();
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
 * @tparam Base_v logical for whether `Checker` should pass for
 *   `EigenBase` types.
 * @tparam SparseCompressed_v logical for whether `Checker` should pass
 * for `SparseCompressedBase` types.
 * @tparam SparseMatrix_v logical for whether `Checker` should pass for
 *   `SparseMatrixBase` types.
 * @tparam SparseMap_v logical for whether `Checker` should pass for
 *   `SparseMapBase` types.
 * @tparam Checker A struct taking in an expression template type and returning
 *  a true if check is satisfied and false otherwise.
 */
template <bool Base_v, bool SparseCompressed_v, bool SparseMatrix_v,
          bool SparseMap_v, template <class...> class Checker>
void all_eigen_sparse() {
  using stan::math::fvar;
  using stan::math::var;

  test_eigen_sparse_matrix<Base_v, SparseCompressed_v, SparseMatrix_v,
                           SparseMap_v, double, Checker>();
  test_eigen_sparse_matrix<Base_v, SparseCompressed_v, SparseMatrix_v,
                           SparseMap_v, var, Checker>();
  test_eigen_sparse_matrix<Base_v, SparseCompressed_v, SparseMatrix_v,
                           SparseMap_v, fvar<double>, Checker>();
  test_eigen_sparse_matrix<Base_v, SparseCompressed_v, SparseMatrix_v,
                           SparseMap_v, fvar<var>, Checker>();
  test_eigen_sparse_matrix<Base_v, SparseCompressed_v, SparseMatrix_v,
                           SparseMap_v, fvar<fvar<var>>, Checker>();
}
}  // namespace test
}  // namespace math
}  // namespace stan

#endif
