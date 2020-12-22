#ifndef TEST_UNIT_MATH_TEST_AD_MATVAR_HPP
#define TEST_UNIT_MATH_TEST_AD_MATVAR_HPP

#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>
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

/**
 * Runs `expect_near_rel` over the values and adjoints of
 *  `std::vector<var_value<T>>` types.
 * @tparam T1 A `var_value` with any template parameter
 * @tparam T2 A `var_value` with any template parameter
 * @param x A `var_value` whose values and adjoints are compared against `y`s
 * @param y A `var_value` whose values and adjoints are compared against `x`s
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
 * Runs `expect_near_rel` over the values and adjoints of `var_value` types.
 * @tparam T1 A `var_value` with any template parameter
 * @tparam T2 A `var_value` with any template parameter
 * @param x A `var_value` whose values and adjoints are compared against `y`s
 * @param y A `var_value` whose values and adjoints are compared against `x`s
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
 * Runs `expect_near_rel` over arithmetic input values.
 * @tparam T1 A type with an arithmetic scalar type.
 * @tparam T2 A type with an arithmetic scalar type.
 * @param x A `var_value` whose values and adjoints are compared against `y`s
 * @param y A `var_value` whose values and adjoints are compared against `x`s
 */
template <typename T1, typename T2,
          require_any_st_arithmetic<T1, T2>* = nullptr>
void expect_near_rel_matvar(const std::string& message, T1&& x, T2&& y,
                         const ad_tolerances& tols) {
  expect_near_rel(message + std::string(" doubles"),
		  stan::math::value_of(x), stan::math::value_of(y), tols.gradient_val_);
}

template <typename... T1, typename... T2, std::size_t... I>
void expect_near_rel_matvar(const std::string& message,
			    const std::tuple<T1...>& x,
			    const std::tuple<T2...>& y,
			    const ad_tolerances& tols,
			    std::index_sequence<I...> i) {
  std::vector<int> A({ (expect_near_rel_matvar(message, std::get<I>(x), std::get<I>(y), tols), 0) ...});
}

template <typename... T1, typename... T2>
void expect_near_rel_matvar(const std::string& message,
			    const std::tuple<T1...>& x,
			    const std::tuple<T2...>& y,
			    const ad_tolerances& tols) {
  if(!(sizeof...(T1) == sizeof...(T2))) {
    FAIL() << "The number of arguments in each tuple must match";
  }

  expect_near_rel_matvar(message, x, y, tols,
			 std::make_index_sequence<sizeof...(T1)>());
}

/**
 * Test that the jacobian for matrices of vars is equal for the
 *  var matrix when the result is a `var` type.
 * @tparam MatVar1 An Eigen type inheriting from `EigenBase` that has
 *  Scalar vars.
 * @tparam MatVar2 An Eigen type inheriting from `EigenBase` that has
 *  Scalar vars.
 * @tparam VarMat1 A `var_value` with an inner type inheriting from `EigenBase`.
 * @tparam VarMat2 A `var_value` with an inner type inheriting from `EigenBase`.
 * @tparam ResultMatVar The resulting type of applying a function to `A_mv1`
 * and `A_mv2`.
 * @tparam ResultVarMat The result type of applying a function to `A_vm1` and
 *  `A_vm2`.
 * @param tols The relative tolerances between the two. Uses `gradient_val_`.
 * @param A_mv1 A eigen type of vars.
 * @param A_mv2 A eigen type of vars.
 * @param A_vm1 A eigen type of vars.
 * @param A_vm2 A eigen type of vars.
 * @param A_mv_f The result of a function from applying `A_mv1` and `A_mv2`.
 * @param A_vm_f The result of a function from applying `A_vm1` and `A_vm2`
 */
template <typename ResultMatVar, typename ResultVarMat,
	  typename... MatVarArgs,
          typename... VarMatArgs,
          require_all_var_vt<std::is_arithmetic, ResultMatVar,
                             ResultVarMat>* = nullptr>
inline void test_matvar_gradient(const ad_tolerances& tols,
                                 ResultMatVar& A_mv_f, ResultVarMat& A_vm_f,
                                 const std::tuple<MatVarArgs...>& A_mv_tuple,
                                 const std::tuple<VarMatArgs...>& A_vm_tuple) {
  A_vm_f.adj() = 1;
  A_mv_f.adj() = 1;
  stan::math::grad();
  expect_near_rel_matvar("var<Matrix> vs Matrix<var>", A_vm_tuple, A_mv_tuple, tols);
  stan::math::set_zero_all_adjoints();
}

/**
 * Test that the jacobian for matrices of vars is equal for the
 *  var matrix when the result is a std::vector.
 *
 * @tparam MatVar An Eigen type inheriting from `EigenBase` that has
 * Scalar vars.
 * @param tols The relative tolerances between the two. Uses `gradient_val_`.
 * @tparam VarMat A `var_value` with an inner type inheriting from `EigenBase`.
 * @tparam ResultMatVar The resulting type of applying a function to `A_mv`.
 * @tparam ResultVarMat The result type of applying a function to `A_vm`.
 * @param A_mv A eigen type of vars.
 * @param A_vm A eigen type of vars.
 * @param A_mv_f The result of a function from applying `A_mv`.
 * @param A_vm_f The result of a function from applying `A_vm`
 */
template <typename ResultMatVar, typename ResultVarMat,
	  typename... MatVarArgs,
          typename... VarMatArgs,
          require_std_vector_vt<is_var, ResultMatVar>* = nullptr,
          require_std_vector_vt<is_var, ResultVarMat>* = nullptr>
inline void test_matvar_gradient(const ad_tolerances& tols,
                                 ResultMatVar& A_mv_f, ResultVarMat& A_vm_f,
                                 const std::tuple<MatVarArgs...>& A_mv_tuple,
                                 const std::tuple<VarMatArgs...>& A_vm_tuple) {
  for (size_t i = 0; i < A_vm_f.size(); ++i) {
    A_vm_f[i].adj() = 1;
    A_mv_f[i].adj() = 1;
    stan::math::grad();
    expect_near_rel_matvar("var<Matrix> vs Matrix<var>", A_vm_tuple, A_mv_tuple, tols);
    stan::math::set_zero_all_adjoints();
  }
}

/**
 * Test that the jacobian for matrices of vars is equal for the
 *  var matrix when the result is an either a `var_value` with inner Eigen type
 * or an Eigen matrix with a scalar `var` type.
 * @tparam MatVar1 An Eigen type inheriting from `EigenBase` that has
 *  Scalar vars.
 * @tparam MatVar2 An Eigen type inheriting from `EigenBase` that has
 *  Scalar vars.
 * @tparam VarMat1 A `var_value` with an inner type inheriting from `EigenBase`.
 * @tparam VarMat2 A `var_value` with an inner type inheriting from `EigenBase`.
 * @tparam ResultMatVar The resulting type of applying a function to `A_mv1`
 * and `A_mv2`.
 * @tparam ResultVarMat The result type of applying a function to `A_vm1` and
 *  `A_vm2`.
 * @param A_mv1 A eigen type of vars.
 * @param A_mv2 A eigen type of vars.
 * @param A_vm1 A eigen type of vars.
 * @param A_vm2 A eigen type of vars.
 * @param A_mv_f The result of a function from applying `A_mv1` and `A_mv2`.
 * @param A_vm_f The result of a function from applying `A_vm1` and `A_vm2`
 */
template <typename ResultMatVar, typename ResultVarMat,
	  typename... MatVarArgs,
          typename... VarMatArgs,//matrix_dynamic
          require_eigen_t<ResultMatVar>* = nullptr>
inline void test_matvar_gradient(const ad_tolerances& tols,
                                 ResultMatVar& A_mv_f, ResultVarMat& A_vm_f,
                                 const std::tuple<MatVarArgs...>& A_mv_tuple,
                                 const std::tuple<VarMatArgs...>& A_vm_tuple) {
  for (Eigen::Index i = 0; i < A_mv_f.rows(); ++i) {
    for (Eigen::Index j = 0; j < A_mv_f.cols(); ++j) {
      A_vm_f.adj()(i, j) = 1;
      A_mv_f.adj()(i, j) = 1;
      stan::math::grad();
      expect_near_rel_matvar("var<Matrix> vs Matrix<var>", A_vm_tuple, A_mv_tuple, tols);
      stan::math::set_zero_all_adjoints();
    }
  }
}

  /*template <typename ResultMatVar, typename ResultVarMat,
	  typename... MatVarArgs,
          typename... VarMatArgs,
          require_eigen_vector_t<ResultMatVar>* = nullptr>
inline void test_matvar_gradient(const ad_tolerances& tols,
                                 ResultMatVar& A_mv_f, ResultVarMat& A_vm_f,
                                 const std::tuple<MatVarArgs...>& A_mv_tuple,
                                 const std::tuple<VarMatArgs...>& A_vm_tuple) {
  for (Eigen::Index i = 0; i < A_mv_f.size(); ++i) {
    A_vm_f.adj()(i) = 1;
    A_mv_f.adj()(i) = 1;
    stan::math::grad();
    expect_near_rel_matvar("var<Matrix> vs Matrix<var>", A_vm_tuple, A_mv_tuple, tols);
    stan::math::set_zero_all_adjoints();
  }
  }*/

/**
 * Test that the jacobian for matrices of vars is equal for the
 *  var matrix when the result is a std::vector.
 *
 * @tparam MatVar An Eigen type inheriting from `EigenBase` that has
 * Scalar vars.
 * @param tols The relative tolerances between the two. Uses `gradient_val_`.
 * @tparam VarMat A `var_value` with an inner type inheriting from `EigenBase`.
 * @tparam ResultMatVar The resulting type of applying a function to `A_mv`.
 * @tparam ResultVarMat The result type of applying a function to `A_vm`.
 * @param A_mv A eigen type of vars.
 * @param A_vm A eigen type of vars.
 * @param A_mv_f The result of a function from applying `A_mv`.
 * @param A_vm_f The result of a function from applying `A_vm`
 */
template <typename ResultMatVar, typename ResultVarMat, typename MatVar,
          typename VarMat,
          require_std_vector_vt<is_eigen, ResultMatVar>* = nullptr,
          require_std_vector_vt<is_var_matrix, ResultVarMat>* = nullptr>
inline void test_matvar_gradient(const ad_tolerances& tols,
                                 ResultMatVar& A_mv_f, ResultVarMat& A_vm_f,
                                 const MatVar& A_mv, const VarMat& A_vm) {
  for (size_t i = 0; i < A_vm_f.size(); ++i) {
    test_matvar_gradient(tols, A_mv_f[i], A_vm_f[i], A_mv, A_vm);
  }
}

/**
 * Given at least one var<matrix> input expect a `var<matrix>` output and
 * `matrix<var>` when both inputs to a function are `matrix<var>`.
 * @tparam ReturnType The result of applying a `x` and `y` to a function `f()`
 * @tparam Type1 the type of the first argument to the function `f()`.
 * @tparam Type1 the type of the second argument to the function `f()`.
 * @param ret not used, only types are used.
 * @param x not used, only types are used.
 * @param y not used, only types are used.
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
    FAIL() << "All the arguments are `Eigen::Matrix<T, R, C>` types, but the return type is neither an `Eigen::Matrix<T, R2, C2>` or a `var`";
  } else if (stan::math::disjunction<is_var_matrix<Types>...>::value
             && !(is_var_matrix<ReturnType>::value
                  || std::is_same<ReturnType, var_value<double>>::value)) {
      FAIL() << "One of the arguments is a `var_value<Eigen::Matrix<double, R, C>>`, but the return type is not a `var_value<T>`";
  }
}

/**
 * Given the input and outputs are `std::vector` types, check that if
 * the input value type is a var<matrix> then the output value_type should be
 * a `var<matrix>` as well. Similarly if the input value_type is
 * `matrix<var>`, the output value type should be `matrix<var>`.
 *
 * @tparam ReturnType The result of applying `x` to a function `f()`
 * @tparam Type the type of the input to the function `f()`.
 * @param ret output of function
 * @param x input of function
 */
template <typename ReturnType, typename Type,
          require_std_vector_t<ReturnType>* = nullptr,
          require_std_vector_vt<is_matrix, Type>* = nullptr>
void check_return_type(const ReturnType& ret, const Type& x) {
  if (ret.size() > 0 && x.size() > 0)
    check_return_type(ret[0], x[0]);
}

/**
 * For an unary function check that an Eigen matrix of vars and a var with an
 * inner eigen type return the same values and adjoints to within
 * the given tolerances. This is done by calling the function on both types,
 * summing the results, then taking the gradient of the sum which will
 * propogate the adjoints upwards.
 *
 * @tparam F A lambda or functor type that calls the unary function.
 * @tparam T An Eigen matrix type
 * @param tols tolerances for test
 * @param f a lambda
 * @param x An Eigen matrix.
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
  vec_mat_var A_mv;
  for (auto xi : x) {
    A_mv.push_back(xi);
  }
  plain_type_t<decltype(f(A_mv))> A_mv_f;
  vec_var_mat A_vm;
  for (auto xi : x) {
    A_vm.push_back(xi);
  }
  plain_type_t<decltype(f(A_vm))> A_vm_f;
  // Check return type is correct
  check_return_type(A_mv_f, A_mv);
  check_return_type(A_vm_f, A_vm);
  // If one throws, the other should throw as well
  try {
    A_mv_f = f(A_mv);
  } catch (...) {
    try {
      f(A_vm);
      FAIL() << type_name<vec_mat_var>() << " throws, expected "
             << type_name<vec_var_mat>() << " version to throw";
    } catch (...) {
      SUCCEED();
      return;
    }
  }
  try {
    A_vm_f = f(A_vm);
  } catch (...) {
    try {
      f(A_mv);
      FAIL() << type_name<vec_var_mat>() << " throws, expected "
             << type_name<vec_mat_var>() << " version to throw";
    } catch (...) {
      SUCCEED();
      return;
    }
  }
  for (size_t i = 0; i < x.size(); ++i) {
    if (!x[i].allFinite()) {
      return;
    }
  }
  test_matvar_gradient(tols, A_mv_f, A_vm_f, std::make_tuple(A_mv), std::make_tuple(A_vm));
}

/**
 * For binary function check that an Eigen matrix of vars and a var with an
 * inner eigen type return the same values and adjoints. This is done by
 * calling the function on both types, summing the results, then taking the
 * gradient of the sum which will propogate the adjoints upwards.
 *
 * @tparam F A lambda or functor type that calls the unary function.
 * @tparam Type1 Either arithmetic or `var`. If arithmetic This will make the
 * first argument of the function invocation be arithmetic, if `var` this will
 * make the first argument be an Eigen matrix of vars and a var with an inner
 * matrix type.
 * @tparam Type2 Either arithmetic or `var`. If arithmetic This will make the
 * second argument of the function invocation be arithmetic, if `var` this will
 * make the second argument be an Eigen matrix of vars and a var with an inner
 * matrix type.
 * @tparam EigMat1 An eigen type
 * @tparam EigMat2 An eigen type
 * @param f a lambda
 * @param x An Eigen matrix.
 * @param y An Eigen matrix.
 */
template <typename... Types, typename F,
	  typename... EigMats>
void expect_ad_matvar_impl(const ad_tolerances& tols, const F& f,
			   const EigMats&... x) {
  using stan::is_var;
  using stan::plain_type_t;
  using stan::math::promote_scalar_t;
  using stan::math::var;
  using stan::math::var_value;
  using stan::math::test::type_name;

  auto A_mv_tuple =
    std::make_tuple(promote_scalar_t<scalar_type_t<Types>, EigMats>(x)...);
  auto A_vm_tuple =
    std::make_tuple(return_var_matrix_t<EigMats, Types>(x)...);

  decltype(stan::math::apply(f, A_mv_tuple)) A_mv_f;
  decltype(stan::math::apply(f, A_vm_tuple)) A_vm_f;

  stan::math::apply([&](auto&&... args) {
      check_return_type(A_mv_f, args...);
    }, A_mv_tuple);
  stan::math::apply([&](auto&&... args) {
      check_return_type(A_vm_f, args...);
    }, A_vm_tuple);

  // If one throws, the other should throw as well
  try {
    A_mv_f = stan::math::apply(f, A_mv_tuple);
  } catch (...) {
    try {
      stan::math::apply(f, A_vm_tuple);
      FAIL() << "`Eigen::Matrix<var, R, C>` function throws and `var_value<Eigen::Matrix<double, R, C>>` does not";
    } catch (...) {
      SUCCEED();
      return;
    }
  }
  try {
    A_vm_f = stan::math::apply(f, A_vm_tuple);
  } catch (...) {
    try {
      stan::math::apply(f, A_mv_tuple);
      FAIL() << "`var_value<Eigen::Matrix<double, R, C>>` function throws and `Eigen::Matrix<var, R, C>` does not";
    } catch (...) {
      SUCCEED();
      return;
    }
  }
  
  //if (!stan::math::is_scal_finite(x)) {
  //  return;
  //}
  
  test_matvar_gradient(tols, A_mv_f, A_vm_f, A_mv_tuple, A_vm_tuple);

  stan::math::recover_memory();
}

/**
 * For an unary function check that an Eigen matrix of vars and a var with an
 * inner eigen type return the same values and adjoints. This is done by
 * calling the function on both types, summing the results, then taking the
 * gradient of the sum which will propogate the adjoints upwards.
 *
 * @tparam F A lambda or functor type that calls the unary function.
 * @tparam T An Eigen matrix type
 * @param f a lambda
 * @param x An Eigen matrix.
 */
template <typename F, typename EigMat>
void expect_ad_matvar(const ad_tolerances& tols, const F& f, const EigMat& x) {
  using varmat = stan::math::var_value<Eigen::MatrixXd>;

  expect_ad_matvar_impl<varmat>(tols, f, x.eval());

  stan::math::recover_memory();
}

template <typename F, typename EigMat>
void expect_ad_matvar(const F& f, const EigMat& x) {
  ad_tolerances tols;

  expect_ad_matvar(tols, f, x);
}

/**
 * For an unary function check that a `std::vector` of Eigen matrix of vars
 * and a `std::vector` of vars with inner eigen type return the same values
 * and adjoints.
 *
 * @tparam F A lambda or functor type that calls the unary function.
 * @tparam EigMat An Eigen matrix type
 * @param f a lambda
 * @param x An Eigen matrix.
 */
template <typename F, typename EigMat>
void expect_ad_matvar(const F& f, const std::vector<EigMat>& x) {
  ad_tolerances tols;
  expect_ad_matvar_v(tols, f, x);
}

template <typename F, typename EigMat>
void expect_ad_matvar(const ad_tolerances& tols, const F& f,
                      const std::vector<EigMat>& x) {
  expect_ad_matvar_v(tols, f, x);
}

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
 * For binary function check that an Eigen matrix of vars and a var with an
 * inner eigen type return the same values and adjoints. This is done by
 * calling the function on both types, summing the results, then taking the
 * gradient of the sum which will propogate the adjoints upwards.
 *
 * @tparam F A lambda or functor type that calls the unary function.
 * @tparam EigMat1 An eigen type
 * @tparam EigMat2 An eigen type
 * @param f a lambda
 * @param x An Eigen matrix.
 * @param y An Eigen matrix.
 */
template <typename F, typename EigMat1, typename EigMat2>
void expect_ad_matvar(const ad_tolerances& tols, const F& f, const EigMat1& x,
                      const EigMat2& y) {
  using stan::math::var;
  using varmat = stan::math::var_value<Eigen::MatrixXd>;

  const auto& x_eval = x.eval();
  const auto& y_eval = y.eval();

  expect_ad_matvar_impl<double, var>(tols, f, x_eval, y_eval);
  expect_ad_matvar_impl<double, varmat>(tols, f, x_eval, y_eval);
  expect_ad_matvar_impl<var, double>(tols, f, x_eval, y_eval);
  expect_ad_matvar_impl<var, var>(tols, f, x_eval, y_eval);
  expect_ad_matvar_impl<var, varmat>(tols, f, x_eval, y_eval);
  expect_ad_matvar_impl<varmat, double>(tols, f, x_eval, y_eval);
  expect_ad_matvar_impl<varmat, var>(tols, f, x_eval, y_eval);
  expect_ad_matvar_impl<varmat, varmat>(tols, f, x_eval, y_eval);

  stan::math::recover_memory();
}

/**
 * For binary function check that an Eigen matrix of vars and a var with an
 * inner eigen type return the same values and adjoints. This is done by
 * calling the function on both types, summing the results, then taking the
 * gradient of the sum which will propogate the adjoints upwards.
 *
 * @tparam F A lambda or functor type that calls the unary function.
 * @tparam EigMat1 An eigen type
 * @tparam EigMat2 An eigen type
 * @param f a lambda
 * @param x An Eigen matrix.
 * @param y An Eigen matrix.
 */
template <typename F, typename EigMat1, typename EigMat2>
void expect_ad_matvar(const F& f, const EigMat1& x, const EigMat2& y) {
  ad_tolerances tols;
  expect_ad_matvar(tols, f, x, y);
}

}  // namespace test
}  // namespace stan
#endif
