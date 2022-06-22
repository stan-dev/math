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
template <
    typename T1, typename T2,
    require_all_not_std_vector_t<value_type_t<T1>, value_type_t<T2>>* = nullptr,
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

template <typename T1, typename T2,
          require_all_std_vector_vt<is_std_vector, T1, T2>* = nullptr,
          require_all_std_vector_st<is_var, T1, T2>* = nullptr>
void expect_near_rel_matvar(const std::string& message, T1&& x, T2&& y,
                            const ad_tolerances& tols) {
  stan::math::check_size_match("expect_near_rel_var", "x", x.size(), "y",
                               y.size());
  for (size_t i = 0; i < x.size(); ++i) {
    expect_near_rel_matvar(
        message + std::string(" elements at i = ") + std::to_string(i), x[i],
        y[i], tols);
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

namespace internal {
template <size_t... Idx>
inline constexpr auto make_tuple_seq(std::index_sequence<Idx...> idxs) {
  return std::make_tuple(Idx...);
}
}  // namespace internal

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
  constexpr auto idxs = internal::make_tuple_seq(
      std::make_index_sequence<std::tuple_size<std::tuple<T1...>>::value>{});
  stan::math::for_each(
      [&message, &tols](const auto& x, const auto& y, const auto& i) {
        expect_near_rel_matvar(
            message + std::string(", argument ")
                + std::to_string(i + stan::error_index::value),
            x, y, tols);
      },
      x, y, idxs);
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

/** @name make_matvar_compatible
 * Make the input be matvar compatible with the
 * same scalar type as T.
 */
///@{
/**
 * The scalar type of T is artithmetic, so just pass the input through.
 *
 * @tparam T Target type
 * @tparam S Argument type
 * @param x argument
 * @return x
 */
template <typename T, typename S, require_all_st_arithmetic<T, S>* = nullptr>
auto make_matvar_compatible(const S& x) {
  return x;
}

/**
 * The scalar type of T is var here, so make a copy of
 * x with scalar type promoted to var.
 *
 * @tparam T Target type
 * @tparam S Argument type
 * @param x argument
 * @return x with scalars promoted to var
 */
template <typename T, typename S, require_st_var<T>* = nullptr,
          require_st_arithmetic<S>* = nullptr>
auto make_matvar_compatible(const S& x) {
  return stan::math::promote_scalar_t<stan::math::var, S>(x);
}

/**
 * The scalar type of T is var here, so make a copy of
 * x with scalar type promoted to var.
 *
 * @tparam T Target type
 * @tparam S Argument type
 * @param x argument
 * @return x with scalars promoted to var
 */
template <typename T, typename S, require_st_var<T>* = nullptr,
          require_st_arithmetic<S>* = nullptr>
auto make_matvar_compatible(const std::vector<S>& x) {
  using vec_mat_var
      = std::vector<stan::math::promote_scalar_t<stan::math::var, S>>;
  vec_mat_var A_vec_mv;
  for (auto&& xi : x) {
    A_vec_mv.push_back(make_matvar_compatible<T>(xi));
  }
  return A_vec_mv;
}
///@}

/** @name make_matvar_compatible
 * Make input to varmat compatible if T is a varmat
 * type, otherwise make it be matvar compatible.
 */
///@{
/**
 * T is not a varmat, so just pass the input through to
 * convert to matvar.
 *
 * @tparam T Target type
 * @tparam S Argument type
 * @param x argument
 * @return x
 */
template <typename T, typename S, require_not_var_matrix_t<T>* = nullptr,
          require_st_arithmetic<S>* = nullptr>
auto make_varmat_compatible(const S& x) {
  return make_matvar_compatible<T>(x);
}

/**
 * T is a varmat, so convert x to be a varmat.
 *
 * @tparam T Target type
 * @tparam S Argument type
 * @param x argument
 * @return x as a varmat
 */
template <typename T, typename S, require_var_matrix_t<T>* = nullptr,
          require_st_arithmetic<S>* = nullptr>
auto make_varmat_compatible(const S& x) {
  return stan::math::var_value<plain_type_t<S>>(x);
}

/**
 * T is a varmat, so convert x to be a varmat.
 *
 * @tparam T Target type
 * @tparam S Argument type
 * @param x argument
 * @return x as a varmat
 */
template <typename T, typename S, require_var_matrix_t<T>* = nullptr,
          require_st_arithmetic<S>* = nullptr>
auto make_varmat_compatible(const std::vector<S>& x) {
  using vec_var_mat = std::vector<stan::math::var_value<plain_type_t<S>>>;
  vec_var_mat A_vec_vm;
  for (auto&& xi : x) {
    A_vec_vm.push_back(make_varmat_compatible<T>(xi));
  }
  return A_vec_vm;
}
template <typename T, typename S, require_var_matrix_t<T>* = nullptr,
          require_st_arithmetic<S>* = nullptr>
auto make_varmat_compatible(const std::vector<std::vector<S>>& x) {
  using vec_var_mat
      = std::vector<std::vector<stan::math::var_value<plain_type_t<S>>>>;
  vec_var_mat A_vec_vm;
  for (auto&& xi : x) {
    A_vec_vm.push_back(make_varmat_compatible<T>(xi));
  }
  return A_vec_vm;
}

///@}

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
  using stan::is_eigen;
  using stan::is_var;
  using stan::is_var_matrix;
  using stan::plain_type_t;
  using stan::math::promote_scalar_t;
  using stan::math::var;
  using stan::math::var_value;
  using stan::math::test::type_name;

  if (!stan::math::disjunction<is_var_matrix<Types>...>::value) {
    FAIL() << "expect_ad_matvar requires at least one varmat input!"
           << std::endl;
  }

  if (!stan::math::disjunction<is_var<scalar_type_t<Types>>...>::value) {
    FAIL() << "expect_ad_matvar requires at least one autodiff input!"
           << std::endl;
  }

  auto A_mv_tuple = std::make_tuple(make_matvar_compatible<Types>(x)...);
  auto A_vm_tuple = std::make_tuple(make_varmat_compatible<Types>(x)...);

  bool any_varmat = stan::math::apply(
      [](const auto&... args) {
        return stan::math::disjunction<is_var_matrix<decltype(args)>...>::value
               || stan::math::disjunction<stan::math::conjunction<
                   is_std_vector<decltype(args)>,
                   is_var_matrix<value_type_t<decltype(args)>>>...>::value;
      },
      A_vm_tuple);

  if (!any_varmat) {
    SUCCEED();  // If no varmats are created, skip this test
    return;
  }

  using T_mv_ret = plain_type_t<decltype(stan::math::apply(f, A_mv_tuple))>;
  using T_vm_ret = plain_type_t<decltype(stan::math::apply(f, A_vm_tuple))>;

  T_mv_ret A_mv_ret;
  T_vm_ret A_vm_ret;

  if (is_var_matrix<T_mv_ret>::value
      || (is_std_vector<T_mv_ret>::value
          && is_var_matrix<value_type_t<T_mv_ret>>::value)) {
    FAIL() << "A function with matvar inputs should not return "
           << type_name<T_mv_ret>() << std::endl;
  }

  if (is_eigen<T_vm_ret>::value
      || (is_std_vector<T_vm_ret>::value
          && is_eigen<value_type_t<T_vm_ret>>::value)) {
    FAIL() << "A function with varmat inputs should not return "
           << type_name<T_vm_ret>() << std::endl;
  }

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

/** @name expect_ad_matvar
 * For a function that accepts between one and four scalar, matrix, or
 * `std::vector` arguments (the `std::vector` can contain scalars or matrices),
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
template <typename F, typename EigMat1, typename EigMat2,
          require_all_not_st_integral<EigMat1, EigMat2>* = nullptr>
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

template <typename F, typename EigMat1, typename EigMat2,
          require_st_integral<EigMat1>* = nullptr,
          require_not_st_integral<EigMat2>* = nullptr>
void expect_ad_matvar(const ad_tolerances& tols, const F& f, const EigMat1& x,
                      const EigMat2& y) {
  using stan::math::var;
  using varmat = stan::math::var_value<Eigen::MatrixXd>;

  expect_ad_matvar_impl<value_type_t<EigMat1>, varmat>(tols, f, x, y);
}

template <typename F, typename EigMat1, typename EigMat2,
          require_not_st_integral<EigMat1>* = nullptr,
          require_st_integral<EigMat2>* = nullptr>
void expect_ad_matvar(const ad_tolerances& tols, const F& f, const EigMat1& x,
                      const EigMat2& y) {
  using stan::math::var;
  using varmat = stan::math::var_value<Eigen::MatrixXd>;

  expect_ad_matvar_impl<varmat, int>(tols, f, x, y);
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

/**
 * Overload for quarternary functions with specified tolerances
 *
 * @tparam F Type of function to test
 * @tparam EigMat1 Type of first argument to test
 * @tparam EigMat2 Type of second argument to test
 * @tparam EigMat3 Type of third argument to test
 * @tparam EigMat4 Type of fourth argument to test
 * @param tols Test tolerances
 * @param f Function to test
 * @param x Value of first argument
 * @param y Value of second argument
 * @param z Value of third argument
 * @param q Value of fourth argument
 */
template <typename F, typename EigMat1, typename EigMat2, typename EigMat3,
          typename EigMat4>
void expect_ad_matvar(const ad_tolerances& tols, const F& f, const EigMat1& x,
                      const EigMat2& y, const EigMat3& z, const EigMat4& q) {
  using stan::math::var;
  using varmat = stan::math::var_value<Eigen::MatrixXd>;

  expect_ad_matvar_impl<double, double, double, varmat>(tols, f, x, y, z, q);
  expect_ad_matvar_impl<double, double, var, varmat>(tols, f, x, y, z, q);
  expect_ad_matvar_impl<double, double, varmat, double>(tols, f, x, y, z, q);
  expect_ad_matvar_impl<double, double, varmat, var>(tols, f, x, y, z, q);
  expect_ad_matvar_impl<double, double, varmat, varmat>(tols, f, x, y, z, q);

  expect_ad_matvar_impl<double, var, double, varmat>(tols, f, x, y, z, q);
  expect_ad_matvar_impl<double, var, var, varmat>(tols, f, x, y, z, q);
  expect_ad_matvar_impl<double, var, varmat, double>(tols, f, x, y, z, q);
  expect_ad_matvar_impl<double, var, varmat, var>(tols, f, x, y, z, q);
  expect_ad_matvar_impl<double, var, varmat, varmat>(tols, f, x, y, z, q);

  expect_ad_matvar_impl<double, varmat, double, double>(tols, f, x, y, z, q);
  expect_ad_matvar_impl<double, varmat, double, var>(tols, f, x, y, z, q);
  expect_ad_matvar_impl<double, varmat, double, varmat>(tols, f, x, y, z, q);
  expect_ad_matvar_impl<double, varmat, var, double>(tols, f, x, y, z, q);
  expect_ad_matvar_impl<double, varmat, var, var>(tols, f, x, y, z, q);
  expect_ad_matvar_impl<double, varmat, var, varmat>(tols, f, x, y, z, q);
  expect_ad_matvar_impl<double, varmat, varmat, double>(tols, f, x, y, z, q);
  expect_ad_matvar_impl<double, varmat, varmat, var>(tols, f, x, y, z, q);
  expect_ad_matvar_impl<double, varmat, varmat, varmat>(tols, f, x, y, z, q);

  expect_ad_matvar_impl<var, double, double, varmat>(tols, f, x, y, z, q);
  expect_ad_matvar_impl<var, double, var, varmat>(tols, f, x, y, z, q);
  expect_ad_matvar_impl<var, double, varmat, double>(tols, f, x, y, z, q);
  expect_ad_matvar_impl<var, double, varmat, var>(tols, f, x, y, z, q);
  expect_ad_matvar_impl<var, double, varmat, varmat>(tols, f, x, y, z, q);

  expect_ad_matvar_impl<var, var, double, varmat>(tols, f, x, y, z, q);
  expect_ad_matvar_impl<var, var, var, varmat>(tols, f, x, y, z, q);
  expect_ad_matvar_impl<var, var, varmat, double>(tols, f, x, y, z, q);
  expect_ad_matvar_impl<var, var, varmat, var>(tols, f, x, y, z, q);
  expect_ad_matvar_impl<var, var, varmat, varmat>(tols, f, x, y, z, q);

  expect_ad_matvar_impl<var, varmat, double, double>(tols, f, x, y, z, q);
  expect_ad_matvar_impl<var, varmat, double, var>(tols, f, x, y, z, q);
  expect_ad_matvar_impl<var, varmat, double, varmat>(tols, f, x, y, z, q);
  expect_ad_matvar_impl<var, varmat, var, double>(tols, f, x, y, z, q);
  expect_ad_matvar_impl<var, varmat, var, var>(tols, f, x, y, z, q);
  expect_ad_matvar_impl<var, varmat, var, varmat>(tols, f, x, y, z, q);
  expect_ad_matvar_impl<var, varmat, varmat, double>(tols, f, x, y, z, q);
  expect_ad_matvar_impl<var, varmat, varmat, var>(tols, f, x, y, z, q);
  expect_ad_matvar_impl<var, varmat, varmat, varmat>(tols, f, x, y, z, q);

  expect_ad_matvar_impl<varmat, double, double, double>(tols, f, x, y, z, q);
  expect_ad_matvar_impl<varmat, double, double, var>(tols, f, x, y, z, q);
  expect_ad_matvar_impl<varmat, double, double, varmat>(tols, f, x, y, z, q);
  expect_ad_matvar_impl<varmat, double, var, double>(tols, f, x, y, z, q);
  expect_ad_matvar_impl<varmat, double, var, var>(tols, f, x, y, z, q);
  expect_ad_matvar_impl<varmat, double, var, varmat>(tols, f, x, y, z, q);
  expect_ad_matvar_impl<varmat, double, varmat, double>(tols, f, x, y, z, q);
  expect_ad_matvar_impl<varmat, double, varmat, var>(tols, f, x, y, z, q);
  expect_ad_matvar_impl<varmat, double, varmat, varmat>(tols, f, x, y, z, q);

  expect_ad_matvar_impl<varmat, var, double, double>(tols, f, x, y, z, q);
  expect_ad_matvar_impl<varmat, var, double, var>(tols, f, x, y, z, q);
  expect_ad_matvar_impl<varmat, var, double, varmat>(tols, f, x, y, z, q);
  expect_ad_matvar_impl<varmat, var, var, double>(tols, f, x, y, z, q);
  expect_ad_matvar_impl<varmat, var, var, var>(tols, f, x, y, z, q);
  expect_ad_matvar_impl<varmat, var, var, varmat>(tols, f, x, y, z, q);
  expect_ad_matvar_impl<varmat, var, varmat, double>(tols, f, x, y, z, q);
  expect_ad_matvar_impl<varmat, var, varmat, var>(tols, f, x, y, z, q);
  expect_ad_matvar_impl<varmat, var, varmat, varmat>(tols, f, x, y, z, q);

  expect_ad_matvar_impl<varmat, varmat, double, double>(tols, f, x, y, z, q);
  expect_ad_matvar_impl<varmat, varmat, double, var>(tols, f, x, y, z, q);
  expect_ad_matvar_impl<varmat, varmat, double, varmat>(tols, f, x, y, z, q);
  expect_ad_matvar_impl<varmat, varmat, var, double>(tols, f, x, y, z, q);
  expect_ad_matvar_impl<varmat, varmat, var, var>(tols, f, x, y, z, q);
  expect_ad_matvar_impl<varmat, varmat, var, varmat>(tols, f, x, y, z, q);
  expect_ad_matvar_impl<varmat, varmat, varmat, double>(tols, f, x, y, z, q);
  expect_ad_matvar_impl<varmat, varmat, varmat, var>(tols, f, x, y, z, q);
  expect_ad_matvar_impl<varmat, varmat, varmat, varmat>(tols, f, x, y, z, q);
}

/**
 * Overload for quarternary functions with default tolerances
 *
 * @tparam F Type of function to test
 * @tparam EigMat1 Type of first argument to test
 * @tparam EigMat2 Type of second argument to test
 * @tparam EigMat3 Type of third argument to test
 * @tparam EigMat4 Type of fourth argument to test
 * @param f Function to test
 * @param x Value of first argument
 * @param y Value of second argument
 * @param z Value of third argument
 * @param q Value of fourth argument
 */
template <typename F, typename EigMat1, typename EigMat2, typename EigMat3,
          typename EigMat4>
void expect_ad_matvar(const F& f, const EigMat1& x, const EigMat2& y,
                      const EigMat3& z, const EigMat4& q) {
  ad_tolerances tols;
  expect_ad_matvar(tols, f, x, y, z, q);
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

/**
 * Overload with manually specified tolerances
 *
 * @tparam F Type of function to test
 * @tparam EigVec Test input type
 * @param tols Test tolerances
 * @param f Function to test
 * @param x Test input
 */
template <typename F, typename T1, typename T2,
          require_all_eigen_t<T1, T2>* = nullptr,
          require_all_not_st_integral<T1, T2>* = nullptr>
void expect_ad_vectorized_matvar(const ad_tolerances& tols, const F& f,
                                 const T1& x, const T2& y) {
  auto x_scal = x.coeff(0, 0);
  auto y_scal = y.coeff(0, 0);
  auto x_vec = x.col(0).eval();
  auto y_vec = y.col(0).eval();
  auto x_rowvec = x.col(0).eval();
  auto y_rowvec = y.col(0).eval();

  std::vector<value_type_t<T1>> x_scal_stdvec{x_scal, x_scal};
  std::vector<value_type_t<T2>> y_scal_stdvec{y_scal, y_scal};
  std::vector<std::vector<value_type_t<T1>>> x_scal_stdvec_stdvec{
      x_scal_stdvec, x_scal_stdvec};
  std::vector<std::vector<value_type_t<T2>>> y_scal_stdvec_stdvec{
      x_scal_stdvec, x_scal_stdvec};
  std::vector<Eigen::MatrixXd> x_mat_stdvec{x, x};
  std::vector<Eigen::MatrixXd> y_mat_stdvec{y, y};
  std::vector<std::vector<Eigen::MatrixXd>> x_mat_stdvec_stdvec{x_mat_stdvec,
                                                                x_mat_stdvec};
  std::vector<std::vector<Eigen::MatrixXd>> y_mat_stdvec_stdvec{y_mat_stdvec,
                                                                y_mat_stdvec};
  expect_ad_matvar(tols, f, x_scal, y);  // scal, mat
  expect_ad_matvar(tols, f, x, y_scal);  // mat, scal
  expect_ad_matvar(tols, f, x, y);       // mat, mat
  expect_ad_matvar(tols, f, x_mat_stdvec,
                   y_mat_stdvec);                   // nest<mat>, nest<mat>
  expect_ad_matvar(tols, f, x_mat_stdvec, y_scal);  // nest<mat>, scal
  expect_ad_matvar(tols, f, x_scal, y_mat_stdvec);  // scal, nest<mat>
  expect_ad_matvar(tols, f, x_mat_stdvec_stdvec,
                   y_mat_stdvec_stdvec);  // nest<nest<mat>>, nest<nest<mat>>
  expect_ad_matvar(tols, f, x_mat_stdvec_stdvec,
                   y_scal);  // nest<nest<mat>, scal
  expect_ad_matvar(tols, f, x_scal,
                   y_mat_stdvec_stdvec);  // scal, nest<nest<mat>

  std::vector<Eigen::VectorXd> x_vec_stdvec{x_vec, x_vec};
  std::vector<Eigen::VectorXd> y_vec_stdvec{y_vec, y_vec};
  std::vector<std::vector<Eigen::VectorXd>> x_vec_stdvec_stdvec{x_vec_stdvec,
                                                                x_vec_stdvec};
  std::vector<std::vector<Eigen::VectorXd>> y_vec_stdvec_stdvec{y_vec_stdvec,
                                                                y_vec_stdvec};

  expect_ad_matvar(tols, f, x_vec, y_scal);  // vec, scal
  expect_ad_matvar(tols, f, x_scal, y_vec);  // scal, vec
  expect_ad_matvar(tols, f, x_vec, y_vec);   // vec, vec
  expect_ad_matvar(tols, f, x_vec_stdvec,
                   y_vec_stdvec);                   // nest<vec>, nest<vec>
  expect_ad_matvar(tols, f, x_vec_stdvec, y_scal);  // nest<vec>, scal
  expect_ad_matvar(tols, f, x_scal, y_vec_stdvec);  // scal, nest<vec>
  expect_ad_matvar(tols, f, x_vec_stdvec_stdvec,
                   y_vec_stdvec_stdvec);  // nest<nest<vec>>, nest<nest<vec>>
  expect_ad_matvar(tols, f, x_vec_stdvec_stdvec,
                   y_scal);  // nest<nest<vec>, scal
  expect_ad_matvar(tols, f, x_scal,
                   y_vec_stdvec_stdvec);  // scal, nest<nest<vec>

  std::vector<Eigen::RowVectorXd> x_rowvec_stdvec{x_rowvec, x_rowvec};
  std::vector<Eigen::RowVectorXd> y_rowvec_stdvec{y_rowvec, y_rowvec};
  std::vector<std::vector<Eigen::RowVectorXd>> x_rowvec_stdvec_stdvec{
      x_rowvec_stdvec, x_rowvec_stdvec};
  std::vector<std::vector<Eigen::RowVectorXd>> y_rowvec_stdvec_stdvec{
      y_rowvec_stdvec, y_rowvec_stdvec};

  expect_ad_matvar(tols, f, x_scal, y_rowvec);    // scal, rowvec
  expect_ad_matvar(tols, f, x_rowvec, y_scal);    // rowvec, scal
  expect_ad_matvar(tols, f, x_rowvec, y_rowvec);  // rowvec, rowvec
  expect_ad_matvar(tols, f, x_rowvec_stdvec,
                   y_rowvec_stdvec);                   // nest<vec>, nest<vec>
  expect_ad_matvar(tols, f, x_rowvec_stdvec, y_scal);  // nest<vec>, scal
  expect_ad_matvar(tols, f, x_scal, y_rowvec_stdvec);  // scal, nest<vec>
  expect_ad_matvar(tols, f, x_rowvec_stdvec_stdvec,
                   y_rowvec_stdvec_stdvec);  // nest<nest<vec>>, nest<nest<vec>>
  expect_ad_matvar(tols, f, x_rowvec_stdvec_stdvec,
                   y_scal);  // nest<nest<vec>, scal
  expect_ad_matvar(tols, f, x_scal,
                   y_rowvec_stdvec_stdvec);  // scal, nest<nest<vec>
}

/**
 * Implementation function for testing that binary functions with vector inputs
 * (both var_value<Eigen> and std::vector types) return the same first order
 * derivative as if we were using Eigen<var> inputs.
 *
 * @tparam F type of function
 * @tparam T1 An Eigen matrix of floating point types
 * @tparam T2 An std vector with inner integral type.
 * @param f function to test
 * @param x argument to test
 * @param y argument to test
 */
template <typename F, typename T1, typename T2,
          require_std_vector_vt<std::is_integral, T1>* = nullptr,
          require_eigen_t<T2>* = nullptr>
void expect_ad_vectorized_matvar(const ad_tolerances& tols, const F& f,
                                 const T1& x, const T2& y) {
  auto x_scal = x[0];
  auto y_vec = y.col(0).eval();

  std::vector<T1> x_stdvec{x, x};
  std::vector<T2> y_stdvec{y, y};
  std::vector<decltype(y_vec)> y_stdvec_vec{y_vec, y_vec};
  std::vector<std::vector<T1>> x_stdvec_stdvec{x_stdvec, x_stdvec};
  std::vector<std::vector<T2>> y_stdvec_stdvec{y_stdvec, y_stdvec};
  expect_ad_matvar(tols, f, x[0], y);                 // scal, mat
  expect_ad_matvar(tols, f, x[0], y_vec);             // scal, mat
  expect_ad_matvar(tols, f, x[0], y_stdvec);          // scal, nest<mat>
  expect_ad_matvar(tols, f, x, y_vec);                // stdvec, vec
  expect_ad_matvar(tols, f, x_stdvec, y_stdvec_vec);  // nest<stdvec>, nest<vec>
  expect_ad_matvar(tols, f, x_stdvec, y);             // nest<stdvec>, mat
  expect_ad_matvar(tols, f, x_stdvec_stdvec,
                   y_stdvec);  // nest<nest<stdvec>>, nest<mat>
}

/**
 * Implementation function for testing that binary functions with vector inputs
 * (both var_value<Eigen> and std::vector types) return the same first order
 * derivative as if we were using Eigen<var> inputs.
 *
 * This is a specialisation for use when the second input is an integer type.
 * We reuse the code in the (std::vector<int>, Eigen::Matrix) specialization
 * by writing a lambda that flips the inputs passed to the original lambda.
 *
 * @tparam F type of function
 * @tparam T1 An Eigen matrix of floating point types
 * @tparam T2 An std vector with inner integral type.
 * @param f function to test
 * @param x argument to test
 * @param y argument to test
 */
template <typename F, typename T1, typename T2, require_eigen_t<T1>* = nullptr,
          require_std_vector_vt<std::is_integral, T2>* = nullptr>
void expect_ad_vectorized_matvar(const ad_tolerances& tols, const F& f,
                                 const T1& x, const T2& y) {
  auto g = [&f](const auto& x, const auto& y) { return f(y, x); };
  expect_ad_vectorized_matvar(tols, g, y, x);
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
template <typename F, typename T1, typename T2>
void expect_ad_vectorized_matvar(const F& f, const T1& x, const T2& y) {
  ad_tolerances tols;
  expect_ad_vectorized_matvar(tols, f, x, y);
}
///@}

}  // namespace test
}  // namespace stan
#endif
