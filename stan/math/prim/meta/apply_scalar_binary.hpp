#ifndef STAN_MATH_PRIM_META_APPLY_SCALAR_BINARY_HPP
#define STAN_MATH_PRIM_META_APPLY_SCALAR_BINARY_HPP

#include <stan/math/prim/meta/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/meta/is_stan_scalar.hpp>
#include <stan/math/prim/meta/is_container.hpp>
#include <stan/math/prim/meta/is_eigen.hpp>
#include <stan/math/prim/meta/require_generics.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Base template function for vectorization of binary scalar functions
 * defined by applying a functor to a combination of scalars,
 * containers of matching sizes, or a combination of a scalar and a container.
 * These containers can be a standard library vector, Eigen dense
 * matrix expression template, or container of these. For each specialisation,
 * the same type as the largest (dimension) input is returned.
 *
 * This base template function takes (and returns) scalars.
 *
 * @tparam T1 Type of first argument to which functor is applied.
 * @tparam T2 Type of second argument to which functor is applied.
 * @tparam F Type of functor to apply.
 * @param x First input to which operation is applied.
 * @param y Second input to which operation is applied.
 * @param f functor to apply to inputs.
 * @return Scalar with result of applying functor to input.
 */
template <typename T1, typename T2, typename F,
          require_all_stan_scalar_t<T1, T2>* = nullptr>
inline auto apply_scalar_binary(const T1& x, const T2& y, const F& f) {
  return f(x, y);
}

/**
 * Specialisation for use with two Eigen inputs. Eigen's binaryExpr framework
 * is used for more efficient indexing of both row- and column-major inputs
 * without separate loops.
 *
 * @tparam T1 Type of first argument to which functor is applied.
 * @tparam T2 Type of second argument to which functor is applied.
 * @tparam F Type of functor to apply.
 * @param x First Eigen input to which operation is applied.
 * @param y Second Eigen input to which operation is applied.
 * @param f functor to apply to Eigen input.
 * @return Eigen object with result of applying functor to inputs.
 */
template <typename T1, typename T2, typename F,
          require_all_eigen_t<T1, T2>* = nullptr>
inline auto apply_scalar_binary(const T1& x, const T2& y, const F& f) {
  return x.binaryExpr(y, f).eval();
}

/**
 * Specialisation for use when the first input is an Eigen type and the second
 * is a scalar. Eigen's unaryExpr framework is used for more efficient indexing
 * of both row- and column-major inputs. The unaryExpr framework also allows
 * for the scalar to be captured and applied to each element in the Eigen
 * object.
 *
 * @tparam T1 Type of Eigen object to which functor is applied.
 * @tparam T2 Type of scalar to which functor is applied.
 * @tparam F Type of functor to apply.
 * @param x Eigen input to which operation is applied.
 * @param y Scalar input to which operation is applied.
 * @param f functor to apply to Eigen and scalar inputs.
 * @return Eigen object with result of applying functor to inputs.
 *
 * Note: The return expresssion needs to be evaluated, otherwise the captured
 *         function and scalar fall out of scope.
 */
template <typename T1, typename T2, typename F, require_eigen_t<T1>* = nullptr,
          require_stan_scalar_t<T2>* = nullptr>
inline auto apply_scalar_binary(const T1& x, const T2& y, const F& f) {
  return x.unaryExpr([&f, &y](const auto& v) { return f(v, y); }).eval();
}

/**
 * Specialisation for use when the first input is an scalar and the second is
 * an Eigen type. Eigen's unaryExpr framework is used for more efficient
 * indexing of both row- and column-major inputs. The unaryExpr framework also
 * allows for the scalar to be captured and applied to each element in the
 * Eigen object.
 *
 * @tparam T1 Type of scalar to which functor is applied.
 * @tparam T2 Type of Eigen object to which functor is applied.
 * @tparam F Type of functor to apply.
 * @param x Scalar input to which operation is applied.
 * @param y Eigen input to which operation is applied.
 * @param f functor to apply to Eigen and scalar inputs.
 * @return Eigen object with result of applying functor to inputs.
 *
 * Note: The return expresssion needs to be evaluated, otherwise the captured
 *         function and scalar fall out of scope.
 */
template <typename T1, typename T2, typename F,
          require_stan_scalar_t<T1>* = nullptr, require_eigen_t<T2>* = nullptr>
inline auto apply_scalar_binary(const T1& x, const T2& y, const F& f) {
  return y.unaryExpr([&f, &x](const auto& v) { return f(x, v); }).eval();
}

/**
 * Specialisation for use with (non-nested) std::vectors. Inputs are mapped
 * to Eigen column vectors and then the result is evaluated directly into the
 * returned std::vector (via Eigen::Map).
 *
 * The returned scalar type is deduced to allow for cases where the input and
 * return scalar types differ (e.g., functions implicitly promoting
 * integers).
 *
 * @tparam T1 Type of first std::vector to which functor is applied.
 * @tparam T2 Type of second std::vector to which functor is applied.
 * @tparam F Type of functor to apply.
 * @param x First std::vector input to which operation is applied.
 * @param y Second std::vector input to which operation is applied.
 * @param f functor to apply to std::vector inputs.
 * @return std::vector with result of applying functor to inputs.
 */
template <typename T1, typename T2, typename F,
          require_all_std_vector_vt<is_stan_scalar, T1, T2>* = nullptr>
inline auto apply_scalar_binary(const T1& x, const T2& y, const F& f) {
  decltype(auto) x_vec = as_column_vector_or_scalar(x);
  decltype(auto) y_vec = as_column_vector_or_scalar(y);
  using T_return = value_type_t<decltype(x_vec.binaryExpr(y_vec, f))>;
  std::vector<T_return> result(x.size());
  Eigen::Map<Eigen::Matrix<T_return, -1, 1>>(result.data(), result.size())
      = x_vec.binaryExpr(y_vec, f);
  return result;
}

/**
 * Specialisation for use when the first input is a (non-nested) std::vector
 * and the second is a scalar. The std::vector input is mapped to an Eigen
 * column vector and then the result is evaluated directly into the returned
 * std::vector (via Eigen::Map).
 *
 * The returned scalar type is deduced to allow for cases where the input and
 * return scalar types differ (e.g., functions implicitly promoting
 * integers).
 *
 * @tparam T1 Type of std::vector to which functor is applied.
 * @tparam T2 Type of scalar to which functor is applied.
 * @tparam F Type of functor to apply.
 * @param x std::vector input to which operation is applied.
 * @param y Scalar input to which operation is applied.
 * @param f functor to apply to std::vector and scalar inputs.
 * @return std::vector with result of applying functor to inputs.
 */
template <typename T1, typename T2, typename F,
          require_std_vector_vt<is_stan_scalar, T1>* = nullptr,
          require_stan_scalar_t<T2>* = nullptr>
inline auto apply_scalar_binary(const T1& x, const T2& y, const F& f) {
  decltype(auto) x_vec = as_column_vector_or_scalar(x);
  using T_return = value_type_t<decltype(f(x[0], y))>;
  std::vector<T_return> result(x.size());
  Eigen::Map<Eigen::Matrix<T_return, -1, 1>>(result.data(), result.size())
      = x_vec.unaryExpr([&f, &y](const auto& v) { return f(v, y); });
  return result;
}

/**
 * Specialisation for use when the first input is a scalar and the second is a
 * (non-nested) std::vector. The std::vector input is mapped to an Eigen
 * column vector and then the result is evaluated directly into the returned
 * std::vector (via Eigen::Map).
 *
 * The returned scalar type is deduced to allow for cases where the input and
 * return scalar types differ (e.g., functions implicitly promoting
 * integers).
 *
 * @tparam T1 Type of scalar to which functor is applied.
 * @tparam T2 Type of std::vector to which functor is applied.
 * @tparam F Type of functor to apply.
 * @param x Scalar input to which operation is applied.
 * @param y std::vector input to which operation is applied.
 * @param f functor to apply to std::vector and scalar inputs.
 * @return std::vector with result of applying functor to inputs.
 */
template <typename T1, typename T2, typename F,
          require_stan_scalar_t<T1>* = nullptr,
          require_std_vector_vt<is_stan_scalar, T2>* = nullptr>
inline auto apply_scalar_binary(const T1& x, const T2& y, const F& f) {
  decltype(auto) y_vec = as_column_vector_or_scalar(y);
  using T_return = value_type_t<decltype(f(x, y[0]))>;
  std::vector<T_return> result(y.size());
  Eigen::Map<Eigen::Matrix<T_return, -1, 1>>(result.data(), result.size())
      = y_vec.unaryExpr([&f, &x](const auto& v) { return f(x, v); });
  return result;
}

/**
 * Specialisation for use with two nested containers (std::vectors).
 * The returned scalar type is deduced to allow for cases where the input and
 * return scalar types differ (e.g., functions implicitly promoting
 * integers).
 *
 * @tparam T1 Type of first std::vector to which functor is applied.
 * @tparam T2 Type of second std::vector to which functor is applied.
 * @tparam F Type of functor to apply.
 * @param x First std::vector input to which operation is applied.
 * @param y Second std::vector input to which operation is applied.
 * @param f functor to apply to std::vector inputs.
 * @return std::vector with result of applying functor to inputs.
 */
template <typename T1, typename T2, typename F,
          require_all_std_vector_vt<is_container, T1, T2>* = nullptr>
inline auto apply_scalar_binary(const T1& x, const T2& y, const F& f) {
  using T_return = decltype(apply_scalar_binary(x[0], y[0], f));
  size_t y_size = y.size();
  std::vector<T_return> result(y_size);
  for (size_t i = 0; i < y_size; ++i) {
    result[i] = apply_scalar_binary(x[i], y[i], f);
  }
  return result;
}

/**
 * Specialisation for use when the first input is a nested std::vector and the
 * second is a scalar. The returned scalar type is deduced to allow for cases
 * where the input and return scalar types differ (e.g., functions implicitly
 * promoting integers).
 *
 * @tparam T1 Type of std::vector to which functor is applied.
 * @tparam T2 Type of scalar to which functor is applied.
 * @tparam F Type of functor to apply.
 * @param x std::vector input to which operation is applied.
 * @param y Scalar input to which operation is applied.
 * @param f functor to apply to inputs.
 * @return std::vector with result of applying functor to inputs.
 */
template <typename T1, typename T2, typename F,
          require_std_vector_vt<is_container, T1>* = nullptr,
          require_stan_scalar_t<T2>* = nullptr>
inline auto apply_scalar_binary(const T1& x, const T2& y, const F& f) {
  using T_return = decltype(apply_scalar_binary(x[0], y, f));
  size_t x_size = x.size();
  std::vector<T_return> result(x_size);
  for (size_t i = 0; i < x_size; ++i) {
    result[i] = apply_scalar_binary(x[i], y, f);
  }
  return result;
}

/**
 * Specialisation for use when the first input is a scalar and the second is a
 * nested std::vector. The returned scalar type is deduced to allow for cases
 * where the input and return scalar types differ (e.g., functions implicitly
 * promoting integers).
 *
 * @tparam T1 Type of scalar to which functor is applied.
 * @tparam T2 Type of std::vector to which functor is applied.
 * @tparam F Type of functor to apply.
 * @param x Scalar input to which operation is applied.
 * @param y std::vector input to which operation is applied.
 * @param f functor to apply to inputs.
 * @return std::vector with result of applying functor to inputs.
 */
template <typename T1, typename T2, typename F,
          require_stan_scalar_t<T1>* = nullptr,
          require_std_vector_vt<is_container, T2>* = nullptr>
inline auto apply_scalar_binary(const T1& x, const T2& y, const F& f) {
  using T_return = decltype(apply_scalar_binary(x, y[0], f));
  size_t y_size = y.size();
  std::vector<T_return> result(y_size);
  for (size_t i = 0; i < y_size; ++i) {
    result[i] = apply_scalar_binary(x, y[i], f);
  }
  return result;
}

}  // namespace math
}  // namespace stan
#endif
