#ifndef TEST_UNIT_MATH_TEST_AD_MATVAR_HPP
#define TEST_UNIT_MATH_TEST_AD_MATVAR_HPP

#include <stan/math/mix.hpp>
#include <test/unit/pretty_print_types.hpp>
#include <test/unit/math/ad_tolerances.hpp>
#include <test/unit/math/is_finite.hpp>
#include <test/unit/math/expect_near_rel.hpp>
#include <test/unit/math/serializer.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <algorithm>
#include <string>
#include <vector>

namespace stan {
namespace test {

/** @name expect_near_rel_matvar
 * Use `expect_near_rel` to check that the values and adjoints of the
 *  input arguments are equal.
 */
///@{
/**
 * Overload for two `std::vector` arguments
 *
 * @tparam T1 Type of first argument, will be `std::vector<var_value<T>>`
 * @tparam T2 Type of second argument, will be `std::vector<var_value<T>>`
 * @param message Debug message to print on fail
 * @param x First argument to compare
 * @param y Second argument to compare
 * @param tols Tolerances for comparison
 */
template <typename T1, typename T2,
          require_all_std_vector_st<is_var, T1, T2>* = nullptr>
void expect_near_rel_matvar(const std::string& message, T1&& x, T2&& y,
                            const ad_tolerances& tols) {
  stan::math::check_size_match("expect_near_rel_var", "x", x.size(), "y",
                               y.size());
  for (size_t i = 0; i < x.size(); ++i) {
    expect_near_rel(
        message + std::string(" values at i = ") + std::to_string(i),
        x[i].val(), y[i].val(), tols.gradient_val_);
    expect_near_rel(
        message + std::string(" adjoints at i = ") + std::to_string(i),
        x[i].adj(), y[i].adj(), tols.gradient_grad_varmat_matvar_);
  }
}

/**
 * Overload for non-`std::vector` arguments with scalar type `var`
 *
 * @tparam T1 Type of first argument, will have scalar type `var`
 * @tparam T2 Type of second argument, will have scalar type `var`
 * @param message Debug message to print on fail
 * @param x First argument to compare
 * @param y Second argument to compare
 * @param tols Tolerances for comparison
 */
template <typename T1, typename T2, require_all_st_var<T1, T2>* = nullptr,
          require_all_not_std_vector_t<T1, T2>* = nullptr>
void expect_near_rel_matvar(const std::string& message, T1&& x, T2&& y,
                            const ad_tolerances& tols) {
  expect_near_rel(message + std::string(" values"), x.val(), y.val(),
                  tols.gradient_val_);
  expect_near_rel(message + std::string(" adjoints"), x.adj(), y.adj(),
                  tols.gradient_grad_varmat_matvar_);
}

/**
 * Overload for arguments with arithmetic scalar type. This only checks
 * values since the adjoints don't exist
 *
 * @tparam T1 Type of first argument, will have arithmetic scalar type
 * @tparam T2 Type of second argument, will have arithmetic scalar type
 * @param message Debug message to print on fail
 * @param x First argument to compare
 * @param y Second argument to compare
 * @param tols Tolerances for comparison
 */
template <typename T1, typename T2,
          require_any_st_arithmetic<T1, T2>* = nullptr>
void expect_near_rel_matvar(const std::string& message, T1&& x, T2&& y,
                            const ad_tolerances& tols) {
  expect_near_rel(message + std::string(" doubles"), stan::math::value_of(x),
                  stan::math::value_of(y), tols.gradient_val_);
}

/**
 * Overload for tuples of arguments. This recursively calls
 * `expect_near_rel_matvar` on each input pair
 *
 * @tparam T1 Tuple type of first argument
 * @tparam T2 Tuple Type of second argument
 * @param message Debug message to print on fail
 * @param x First argument to compare
 * @param y Second argument to compare
 * @param tols Tolerances for comparison
 */
template <typename... T1, typename... T2>
void expect_near_rel_matvar(const std::string& message,
                            const std::tuple<T1...>& x,
                            const std::tuple<T2...>& y,
                            const ad_tolerances& tols) {
  if (!(sizeof...(T1) == sizeof...(T2))) {
    FAIL() << "The number of arguments in each tuple must match";
  }

  stan::math::for_each(
      [&message, &tols](const auto& x, const auto& y, const auto& i) {
        expect_near_rel_matvar(
            message + std::string(", argument ")
                + std::to_string(i + stan::error_index::value),
            x, y, tols);
      },
      x, y);
}
///@}

/** @name test_matvar_gradient
 * Given the inputs and outputs of a function called with Eigen
 * matrices of vars (matvars) and vars with inner Eigen types (varmats),
 * check that the values of the outputs are equal and that gradients with
 * respect to each output are the same
 *
 * @tparam ResultMV Output type of matvar calculation
 * @tparam ResultVM Output type of varmat calculation
 * @tparam MatVarArgs Argument types for matvar calculation
 * @tparam VarMatArgs Argument types for varmat calculation
 * @param tols Test tolerances
 * @param A_mv_ret Output for matvar calculation
 * @param A_vm_ret Output for varmat calculation
 * @param A_mv_tuple Input arguments to matvar calculation as tuple
 * @param A_vm_tuple Input arguments to varmat calculation as tuple
 */
///@{
/**
 * Overload for scalar outputs
 */
template <typename ResultMV, typename ResultVM, typename... MatVarArgs,
          typename... VarMatArgs,
          require_all_var_vt<std::is_arithmetic, ResultMV, ResultVM>* = nullptr>
inline void test_matvar_gradient(const ad_tolerances& tols, ResultMV& A_mv_ret,
                                 ResultVM& A_vm_ret,
                                 const std::tuple<MatVarArgs...>& A_mv_tuple,
                                 const std::tuple<VarMatArgs...>& A_vm_tuple) {
  expect_near_rel("var<Matrix> vs Matrix<var> values", value_of(A_vm_ret),
                  value_of(A_mv_ret), tols.gradient_val_);
  A_vm_ret.adj() = 1;
  A_mv_ret.adj() = 1;
  stan::math::grad();
  expect_near_rel_matvar("var<Matrix> vs Matrix<var>", A_vm_tuple, A_mv_tuple,
                         tols);
  stan::math::set_zero_all_adjoints();
}

/**
 * Overload for `std::vector` outputs
 */
template <typename ResultMV, typename ResultVM, typename... MatVarArgs,
          typename... VarMatArgs,
          require_std_vector_vt<is_var, ResultMV>* = nullptr,
          require_std_vector_vt<is_var, ResultVM>* = nullptr>
inline void test_matvar_gradient(const ad_tolerances& tols, ResultMV& A_mv_ret,
                                 ResultVM& A_vm_ret,
                                 const std::tuple<MatVarArgs...>& A_mv_tuple,
                                 const std::tuple<VarMatArgs...>& A_vm_tuple) {
  expect_near_rel("var<Matrix> vs Matrix<var> values", value_of(A_vm_ret),
                  value_of(A_mv_ret), tols.gradient_val_);
  for (size_t i = 0; i < A_vm_ret.size(); ++i) {
    A_vm_ret[i].adj() = 1;
    A_mv_ret[i].adj() = 1;
    stan::math::grad();
    expect_near_rel_matvar("var<Matrix> vs Matrix<var>", A_vm_ret[i],
                           A_mv_ret[i], tols);
    expect_near_rel_matvar("var<Matrix> vs Matrix<var>", A_vm_tuple, A_mv_tuple,
                           tols);
    stan::math::set_zero_all_adjoints();
  }
}

/**
 * Overload for matrix outputs
 */
template <typename ResultMV, typename ResultVM, typename... MatVarArgs,
          typename... VarMatArgs,  // matrix_dynamic
          require_eigen_st<is_var, ResultMV>* = nullptr>
inline void test_matvar_gradient(const ad_tolerances& tols, ResultMV& A_mv_ret,
                                 ResultVM& A_vm_ret,
                                 const std::tuple<MatVarArgs...>& A_mv_tuple,
                                 const std::tuple<VarMatArgs...>& A_vm_tuple) {
  expect_near_rel("var<Matrix> vs Matrix<var> values", value_of(A_vm_ret),
                  value_of(A_mv_ret), tols.gradient_val_);
  for (Eigen::Index j = 0; j < A_mv_ret.cols(); ++j) {
    for (Eigen::Index i = 0; i < A_mv_ret.rows(); ++i) {
      A_vm_ret.adj()(i, j) = 1;
      A_mv_ret.adj()(i, j) = 1;
      stan::math::grad();
      expect_near_rel_matvar("var<Matrix> vs Matrix<var>", A_vm_tuple,
                             A_mv_tuple, tols);
      stan::math::set_zero_all_adjoints();
    }
  }
}

/**
 * Overload for arithmetic outputs -- this is an error if this gets called
 */
template <typename ResultMV, typename ResultVM, typename... MatVarArgs,
          typename... VarMatArgs,  // matrix_dynamic
          require_any_st_arithmetic<ResultMV, ResultVM>* = nullptr,
          require_any_st_var<MatVarArgs..., VarMatArgs...>* = nullptr>
inline void test_matvar_gradient(const ad_tolerances& tols, ResultMV& A_mv_ret,
                                 ResultVM& A_vm_ret,
                                 const std::tuple<MatVarArgs...>& A_mv_tuple,
                                 const std::tuple<VarMatArgs...>& A_vm_tuple) {
  throw std::logic_error(
      "test_matvar_gradient should only be used when the return is a "
      "reverse mode autodiff type. That this got called it means that a "
      "function took an autodiff variable as input yet returned "
      "a non-autodiff type");
}
///@}

/**
 * Given the inputs and outputs of a function called with Eigen
 * matrices of vars (matvars) and vars with inner Eigen types (varmats),
 * check that the values of the outputs are equal and that gradients with
 * respect to each output are the same
 *
 * This overload works when the outputs are `std::vector` types.
 *
 * @tparam ResultMV Output type of matvar calculation
 * @tparam ResultVM Output type of varmat calculation
 * @tparam MatVar Argument type for matvar calculation
 * @tparam VarMat Argument type for varmat calculation
 * @param tols Test tolerances
 * @param A_mv_ret Output for matvar calculation
 * @param A_vm_ret Output for varmat calculation
 * @param A_mv Input argument to matvar calculation
 * @param A_vm Input argument to varmat calculation
 */
template <typename ResultMV, typename ResultVM, typename MatVar,
          typename VarMat,
          require_all_std_vector_t<ResultMV, ResultVM>* = nullptr>
inline void test_matvar_gradient(const ad_tolerances& tols, ResultMV& A_mv_ret,
                                 ResultVM& A_vm_ret, const MatVar& A_mv,
                                 const VarMat& A_vm) {
  for (size_t i = 0; i < A_vm_ret.size(); ++i) {
    test_matvar_gradient(tols, A_mv_ret[i], A_vm_ret[i], A_mv, A_vm);
  }
}

/**
 * Check that if at least one of the input types is a var with inner
 * Eigen type that the output is not an Eigen matrix of vars.
 *
 * Similarly check that if none of the inputs are vars with inner type
 * Eigen that the output is not a var with inner Eigen type.
 *
 * @tparam ReturnType Type of result of calculation
 * @tparam Types Types of input arguments
 * @param ret not used, only types are used.
 * @param x not used, only types are used.
 */
template <typename ReturnType, typename... Types,
          require_all_not_std_vector_t<ReturnType, Types...>* = nullptr,
          require_any_matrix_t<Types...>* = nullptr>
void check_return_type(const ReturnType& ret, const Types&... x) {
  using stan::is_eigen;
  using stan::is_var_matrix;
  using stan::math::var_value;
  using stan::math::test::type_name;
  if (stan::math::conjunction<is_eigen<Types>...>::value
      && !(is_eigen<ReturnType>::value
           || std::is_same<ReturnType, var_value<double>>::value)) {
    FAIL() << "All the arguments are `Eigen::Matrix<T, R, C>` types, but the "
              "return type is neither an `Eigen::Matrix<T, R2, C2>` or a `var`";
  } else if (stan::math::disjunction<is_var_matrix<Types>...>::value
             && !(is_var_matrix<ReturnType>::value
                  || std::is_same<ReturnType, var_value<double>>::value)) {
    FAIL() << "One of the arguments is a `var_value<Eigen::Matrix<double, R, "
              "C>>`, but the return type is not a `var_value<T>`";
  }
}

/**
 * Check that if the input type is a `std::vector` of vars with inner
 * Eigen type that the output is not a `std::vector` of Eigen matrices of vars.
 *
 * Similarly check that if the input is a `std::vector` of Eigen matrices of
 * vars that the output is not a `std::vector` of vars with inner Eigen type.
 *
 * @tparam ReturnType Element types of result of calculation
 * @tparam Types Element types of calculation input
 * @param ret not used, only types are used.
 * @param x not used, only types are used.
 */
template <typename ReturnType, typename Type, require_matrix_t<Type>* = nullptr>
void check_return_type(const std::vector<ReturnType>& ret,
                       const std::vector<Type>& x) {
  check_return_type(ReturnType(), Type());
}

/**
 * For a function that accepts a `std::vector` of matrices,
 * check that calculations with Eigen matrices of vars and vars with
 * inner Eigen type return the same values and adjoints.
 *
 * Functions are expected to either both throw, or both not throw for a given
 * input. The exception types are not checked.
 *
 * Functions must not return Eigen matrices of vars if the input arguments
 * contain vars with inner Eigen type. Similarly they must not return vars with
 * inner type Eigen if the input contains no vars with inner Eigen type.
 *
 * @tparam F Type of function to test
 * @tparam EigMats Types of input values
 * @param tols Test tolerance
 * @param f Function to test
 * @param x Input values to test function at
 */
template <typename F, typename EigMat, require_eigen_t<EigMat>* = nullptr>
void expect_ad_matvar_v(const ad_tolerances& tols, const F& f,
                        const std::vector<EigMat>& x) {
  using stan::plain_type_t;
  using stan::math::promote_scalar_t;
  using stan::math::var;
  using stan::math::var_value;
  using stan::math::test::type_name;
  using vec_mat_var = std::vector<promote_scalar_t<var, EigMat>>;
  using vec_var_mat = std::vector<var_value<plain_type_t<EigMat>>>;
  vec_mat_var A_vec_mv;
  for (auto xi : x) {
    A_vec_mv.push_back(xi);
  }
  plain_type_t<decltype(f(A_vec_mv))> A_mv_ret;
  vec_var_mat A_vec_vm;
  for (auto xi : x) {
    A_vec_vm.push_back(xi);
  }
  plain_type_t<decltype(f(A_vec_vm))> A_vm_ret;
  // Check return type is correct
  check_return_type(A_mv_ret, A_vec_mv);
  check_return_type(A_vm_ret, A_vec_vm);
  // If one throws, the other should throw as well
  try {
    A_mv_ret = f(A_vec_mv);
  } catch (...) {
    try {
      f(A_vec_vm);
      FAIL() << type_name<vec_mat_var>() << " throws, expected "
             << type_name<vec_var_mat>() << " version to throw";
    } catch (...) {
      SUCCEED();
      return;
    }
  }
  try {
    A_vm_ret = f(A_vec_vm);
  } catch (...) {
    try {
      f(A_vec_mv);
      FAIL() << type_name<vec_var_mat>() << " throws, expected "
             << type_name<vec_mat_var>() << " version to throw";
    } catch (...) {
      SUCCEED();
      return;
    }
  }

  test_matvar_gradient(tols, A_mv_ret, A_vm_ret, std::make_tuple(A_vec_mv),
                       std::make_tuple(A_vec_vm));

  stan::math::recover_memory();
}

/**
 * For a function that accepts between any number of scalar or matrix arguments,
 * check that calculations with Eigen matrices of vars and vars with
 * inner Eigen type return the same values and adjoints.
 *
 * Functions are expected to either both throw, or both not throw for a given
 * input. The exception types are not checked.
 *
 * Functions must not return Eigen matrices of vars if the input arguments
 * contain vars with inner Eigen type. Similarly they must not return vars with
 * inner type Eigen if the input contains no vars with inner Eigen type.
 *
 * @tparam Types Types of autodiff variables to use
 * @tparam F Type of function to test
 * @tparam EigMats Types of input values
 * @param tols Test tolerance
 * @param f Function to test
 * @param x Input values to test function at
 */
template <typename... Types, typename F, typename... EigMats>
void expect_ad_matvar_impl(const ad_tolerances& tols, const F& f,
                           const EigMats&... x) {
  using stan::is_var;
  using stan::plain_type_t;
  using stan::math::promote_scalar_t;
  using stan::math::var;
  using stan::math::var_value;
  using stan::math::test::type_name;

  auto A_mv_tuple
      = std::make_tuple(promote_scalar_t<scalar_type_t<Types>, EigMats>(x)...);
  auto A_vm_tuple = std::make_tuple(
      std::conditional_t<is_var<Types>::value,
                         return_var_matrix_t<EigMats, Types>, EigMats>(x)...);

  plain_type_t<decltype(stan::math::apply(f, A_mv_tuple))> A_mv_ret;
  plain_type_t<decltype(stan::math::apply(f, A_vm_tuple))> A_vm_ret;

  stan::math::apply(
      [&](auto&&... args) { check_return_type(A_mv_ret, args...); },
      A_mv_tuple);
  stan::math::apply(
      [&](auto&&... args) { check_return_type(A_vm_ret, args...); },
      A_vm_tuple);

  // If one throws, the other should throw as well
  try {
    A_mv_ret = stan::math::apply(f, A_mv_tuple);
  } catch (...) {
    try {
      stan::math::apply(f, A_vm_tuple);
      FAIL() << "`Eigen::Matrix<var, R, C>` function throws and "
                "`var_value<Eigen::Matrix<double, R, C>>` does not";
    } catch (...) {
      SUCCEED();
      return;
    }
  }
  try {
    A_vm_ret = stan::math::apply(f, A_vm_tuple);
  } catch (...) {
    try {
      stan::math::apply(f, A_mv_tuple);
      FAIL() << "`var_value<Eigen::Matrix<double, R, C>>` function throws and "
                "`Eigen::Matrix<var, R, C>` does not";
    } catch (...) {
      SUCCEED();
      return;
    }
  }

  test_matvar_gradient(tols, A_mv_ret, A_vm_ret, A_mv_tuple, A_vm_tuple);

  stan::math::recover_memory();
}

/** @name expect_ad_matvar_std_vector
 * For a function that accepts a `std::vector` of matrix types, check
 * that calculations with Eigen matrices of vars and vars with
 * inner Eigen type return the same values and adjoints
 */
///@{
/**
 * Overload with manually specified tolerances
 *
 * @tparam F Type of function to test
 * @tparam EigVec Test input type
 * @param tols Test tolerances
 * @param f Function to test
 * @param x Test input
 */
template <typename F, typename EigMat>
void expect_ad_matvar(const ad_tolerances& tols, const F& f,
                      const std::vector<EigMat>& x) {
  expect_ad_matvar_v(tols, f, x);
}

/**
 * Overload with default tolerances
 *
 * @tparam F Type of function to test
 * @tparam EigVec Test input type
 * @param tols Test tolerances
 * @param f Function to test
 * @param x Test input
 */
template <typename F, typename EigMat>
void expect_ad_matvar(const F& f, const std::vector<EigMat>& x) {
  ad_tolerances tols;
  expect_ad_matvar(tols, f, x);
}
///@}

/** @name expect_ad_matvar
 * For a function that accepts between one and three scalar or matrix arguments,
 * check that calculations with Eigen matrices of vars and vars with
 * inner Eigen type return the same values and adjoints.
 *
 * Only instantiations of functions that involve vars with inner Eigen types are
 * tested.
 *
 * For functions with more than one matrix argument, all combinations of Eigen
 * matrices of vars, vars with inner Eigen types are expected to work.
 *
 * All combinations of autodiff and non-autodiff arguments are expected to
 * work as well.
 *
 * Functions are expected to either both throw, or both not throw for a given
 * input. The exception types are not checked.
 */
///@{
/**
 * Overload for unary function with specified tolerances
 *
 * @tparam F Type of function to test
 * @tparam EigMat Type of argument to test
 * @param tols Test tolerances
 * @param f Function to test
 * @param x Value of argument
 */
template <typename F, typename EigMat>
void expect_ad_matvar(const ad_tolerances& tols, const F& f, const EigMat& x) {
  using varmat = stan::math::var_value<Eigen::MatrixXd>;

  expect_ad_matvar_impl<varmat>(tols, f, x);
}

/**
 * Overload for unary function with default tolerances
 *
 * @tparam F Type of function to test
 * @tparam EigMat Type of argument to test
 * @param f Function to test
 * @param x Value of argument
 */
template <typename F, typename EigMat>
void expect_ad_matvar(const F& f, const EigMat& x) {
  ad_tolerances tols;

  expect_ad_matvar(tols, f, x);
}

/**
 * Overload for binary functions with specified tolerances
 *
 * @tparam F Type of function to test
 * @tparam EigMat1 Type of first argument to test
 * @tparam EigMat2 Type of second argument to test
 * @param tols Test tolerances
 * @param f Function to test
 * @param x Value of first argument
 * @param y Value of second argument
 */
template <typename F, typename EigMat1, typename EigMat2>
void expect_ad_matvar(const ad_tolerances& tols, const F& f, const EigMat1& x,
                      const EigMat2& y) {
  using stan::math::var;
  using varmat = stan::math::var_value<Eigen::MatrixXd>;

  expect_ad_matvar_impl<double, varmat>(tols, f, x, y);
  expect_ad_matvar_impl<var, varmat>(tols, f, x, y);
  expect_ad_matvar_impl<varmat, double>(tols, f, x, y);
  expect_ad_matvar_impl<varmat, var>(tols, f, x, y);
  expect_ad_matvar_impl<varmat, varmat>(tols, f, x, y);
}

/**
 * Overload for binary functions with default tolerances
 *
 * @tparam F Type of function to test
 * @tparam EigMat1 Type of first argument to test
 * @tparam EigMat2 Type of second argument to test
 * @param f Function to test
 * @param x Value of first argument
 * @param y Value of second argument
 */
template <typename F, typename EigMat1, typename EigMat2>
void expect_ad_matvar(const F& f, const EigMat1& x, const EigMat2& y) {
  ad_tolerances tols;
  expect_ad_matvar(tols, f, x, y);
}

/**
 * Overload for ternary functions with specified tolerances
 *
 * @tparam F Type of function to test
 * @tparam EigMat1 Type of first argument to test
 * @tparam EigMat2 Type of second argument to test
 * @tparam EigMat3 Type of third argument to test
 * @param tols Test tolerances
 * @param f Function to test
 * @param x Value of first argument
 * @param y Value of second argument
 * @param z Value of third argument
 */
template <typename F, typename EigMat1, typename EigMat2, typename EigMat3>
void expect_ad_matvar(const ad_tolerances& tols, const F& f, const EigMat1& x,
                      const EigMat2& y, const EigMat3& z) {
  using stan::math::var;
  using varmat = stan::math::var_value<Eigen::MatrixXd>;

  expect_ad_matvar_impl<double, double, varmat>(tols, f, x, y, z);
  expect_ad_matvar_impl<double, var, varmat>(tols, f, x, y, z);
  expect_ad_matvar_impl<double, varmat, double>(tols, f, x, y, z);
  expect_ad_matvar_impl<double, varmat, var>(tols, f, x, y, z);
  expect_ad_matvar_impl<double, varmat, varmat>(tols, f, x, y, z);
  expect_ad_matvar_impl<var, double, varmat>(tols, f, x, y, z);
  expect_ad_matvar_impl<var, var, varmat>(tols, f, x, y, z);
  expect_ad_matvar_impl<var, varmat, double>(tols, f, x, y, z);
  expect_ad_matvar_impl<var, varmat, var>(tols, f, x, y, z);
  expect_ad_matvar_impl<var, varmat, varmat>(tols, f, x, y, z);
  expect_ad_matvar_impl<varmat, double, double>(tols, f, x, y, z);
  expect_ad_matvar_impl<varmat, double, var>(tols, f, x, y, z);
  expect_ad_matvar_impl<varmat, double, varmat>(tols, f, x, y, z);
  expect_ad_matvar_impl<varmat, var, double>(tols, f, x, y, z);
  expect_ad_matvar_impl<varmat, var, var>(tols, f, x, y, z);
  expect_ad_matvar_impl<varmat, var, varmat>(tols, f, x, y, z);
  expect_ad_matvar_impl<varmat, varmat, double>(tols, f, x, y, z);
  expect_ad_matvar_impl<varmat, varmat, var>(tols, f, x, y, z);
  expect_ad_matvar_impl<varmat, varmat, varmat>(tols, f, x, y, z);
}

/**
 * Overload for ternary functions with default tolerances
 *
 * @tparam F Type of function to test
 * @tparam EigMat1 Type of first argument to test
 * @tparam EigMat2 Type of second argument to test
 * @tparam EigMat3 Type of third argument to test
 * @param f Function to test
 * @param x Value of first argument
 * @param y Value of second argument
 * @param z Value of third argument
 */
template <typename F, typename EigMat1, typename EigMat2, typename EigMat3>
void expect_ad_matvar(const F& f, const EigMat1& x, const EigMat2& y,
                      const EigMat3& z) {
  ad_tolerances tols;
  expect_ad_matvar(tols, f, x, y, z);
}
///@}

/** @name expect_ad_vector_matvar
 * For a vectorized unary function, check
 * that calculations with Eigen matrices of vars and vars with
 * inner Eigen type return the same values and adjoints on
 * matrices, vectors, row vectors and `std::vector` s of those
 * types.
 */
///@{
/**
 * Overload with manually specified tolerances
 *
 * @tparam F Type of function to test
 * @tparam EigVec Test input type
 * @param tols Test tolerances
 * @param f Function to test
 * @param x Test input
 */
template <typename F, typename EigVec,
          require_eigen_vector_t<EigVec>* = nullptr>
void expect_ad_vector_matvar(const ad_tolerances& tols, const F& f,
                             const EigVec& x) {
  Eigen::VectorXd v = x;
  Eigen::RowVectorXd r = x.transpose();
  Eigen::MatrixXd m(x.size(), 2);
  m.col(0) = x;
  m.col(1) = x.reverse();  // reverse just to mix stuff up

  expect_ad_matvar(f, v);
  expect_ad_matvar(f, r);
  expect_ad_matvar(f, m);

  std::vector<Eigen::VectorXd> vv = {v, v};
  std::vector<Eigen::RowVectorXd> rv = {r, r};
  std::vector<Eigen::MatrixXd> mv = {m, m};

  expect_ad_matvar(f, vv);
  expect_ad_matvar(f, rv);
  expect_ad_matvar(f, mv);
}

/**
 * Overload with default tolerances
 *
 * @tparam F Type of function to test
 * @tparam EigVec Test input type
 * @param tols Test tolerances
 * @param f Function to test
 * @param x Test input
 */
template <typename F, typename EigVec,
          require_eigen_vector_t<EigVec>* = nullptr>
void expect_ad_vector_matvar(const F& f, const EigVec& x) {
  ad_tolerances tols;
  expect_ad_vector_matvar(tols, f, x);
}
///@}

}  // namespace test
}  // namespace stan
#endif
