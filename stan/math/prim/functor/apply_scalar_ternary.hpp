#ifndef STAN_MATH_PRIM_FUNCTOR_APPLY_SCALAR_TERNARY_HPP
#define STAN_MATH_PRIM_FUNCTOR_APPLY_SCALAR_TERNARY_HPP

#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/meta/is_stan_scalar.hpp>
#include <stan/math/prim/meta/is_container.hpp>
#include <stan/math/prim/meta/is_eigen.hpp>
#include <stan/math/prim/meta/require_generics.hpp>
#include <stan/math/prim/err/check_matching_dims.hpp>
#include <stan/math/prim/err/check_matching_sizes.hpp>
#include <stan/math/prim/fun/num_elements.hpp>
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
 * @tparam T1 Type of first scalar to which functor is applied.
 * @tparam T2 Type of second scalar to which functor is applied.
 * @tparam T3 Type of third scalar to which functor is applied.
 * @tparam F Type of functor to apply.
 * @param x First input to which operation is applied.
 * @param y Second input to which operation is applied.
 * @param z Third scalar to which operation is applied.
 * @param f functor to apply to inputs.
 * @return Scalar with result of applying functor to input.
 */
template <typename T1, typename T2, typename T3, typename F,
          require_all_stan_scalar_t<T1, T2, T3>* = nullptr>
inline auto apply_scalar_ternary(const T1& x, const T2& y, const T3& z, const F& f) {
  return f(x, y, z);
}

/**
 * Specialisation for use with three Eigen inputs. Eigen's ternaryExpr framework
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
template <typename T1, typename T2, typename T3, typename F,
          require_all_eigen_t<T1, T2, T3>* = nullptr>
inline auto apply_scalar_ternary(const T1& x, const T2& y, const T3& z, const F& f) {
  check_matching_dims("Binary function", "x", x, "y", y);
  check_matching_dims("Binary function", "x", x, "z", z);
  return x.ternaryExpr(y, z, f);
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
 * @tparam T3 Type of third std::vector to which functor is applied.
 * @tparam F Type of functor to apply.
 * @param x First std::vector input to which operation is applied.
 * @param y Second std::vector input to which operation is applied.
 * @param z Third std::vector input to which operation is applied.
 * @param f functor to apply to std::vector inputs.
 * @return std::vector with result of applying functor to inputs.
 */
template <typename T1, typename T2, typename T3, typename F,
          require_all_std_vector_vt<is_stan_scalar, T1, T2, T3>* = nullptr>
inline auto apply_scalar_ternary(const T1& x, const T2& y, const T3& z, const F& f) {
  check_matching_sizes("Binary function", "x", x, "y", y);
  check_matching_sizes("Binary function", "x", x, "z", z);
  decltype(auto) x_vec = as_column_vector_or_scalar(x);
  decltype(auto) y_vec = as_column_vector_or_scalar(y);
  decltype(auto) z_vec = as_column_vector_or_scalar(z);
  using T_return = std::decay_t<decltype(f(x[0], y[0], z[0]))>;
  std::vector<T_return> result(x.size());
  Eigen::Map<Eigen::Matrix<T_return, -1, 1>>(result.data(), result.size())
      = x_vec.ternaryExpr(y_vec, z_vec, f);
  return result;
}


// T1 Container
template <typename T1, typename T2, typename T3, typename F,
          require_all_container_t<T1, T2>* = nullptr,
          require_stan_scalar_t<T3>* = nullptr>
inline auto apply_scalar_ternary(const T1& x, const T2& y, const T3& z, const F& f) {
  return apply_scalar_binary(x, y, [&f, z](const auto& a, const auto& b) { return f(a, b, z); });
}

template <typename T1, typename T2, typename T3, typename F,
          require_all_container_t<T1, T3>* = nullptr,
          require_stan_scalar_t<T2>* = nullptr>
inline auto apply_scalar_ternary(const T1& x, const T2& y, const T3& z, const F& f) {
  return apply_scalar_binary(x, z, [&f, y](const auto& a, const auto& c) { return f(a, y, c); });
}

template <typename T1, typename T2, typename T3, typename F,
          require_container_t<T1>* = nullptr,
          require_all_stan_scalar_t<T2, T3>* = nullptr>
inline auto apply_scalar_ternary(const T1& x, const T2& y, const T3& z, const F& f) {
  return apply_scalar_binary(x, y, [&](const auto& a, const auto& b) {
    return f(a, b, z);
  });
}

// T2 Container
template <typename T1, typename T2, typename T3, typename F,
          require_all_container_t<T2, T3>* = nullptr,
          require_stan_scalar_t<T1>* = nullptr>
inline auto apply_scalar_ternary(const T1& x, const T2& y, const T3& z, const F& f) {
  return apply_scalar_binary(y, z, [&f, x](const auto& b, const auto& c) { return f(x, b, c); });
}

template <typename T1, typename T2, typename T3, typename F,
          require_container_t<T2>* = nullptr,
          require_all_stan_scalar_t<T1, T3>* = nullptr>
inline auto apply_scalar_ternary(const T1& x, const T2& y, const T3& z, const F& f) {
  return apply_scalar_binary(y, z, [&](const auto& b, const auto& c) { return f(x, b, c); });
}

// T3 Container
template <typename T1, typename T2, typename T3, typename F,
          require_container_t<T3>* = nullptr,
          require_all_stan_scalar_t<T1, T2>* = nullptr>
inline auto apply_scalar_ternary(const T1& x, const T2& y, const T3& z, const F& f) {
  return apply_scalar_binary(y, z, [&](const auto& b, const auto& c) { return f(x, b, c); });
}
// Nested containers
template <typename T1, typename T2, typename T3, typename F,
          require_all_std_vector_vt<is_container_or_var_matrix, T1, T2, T3>* = nullptr>
inline auto apply_scalar_ternary(const T1& x, const T2& y, const T3& z, const F& f) {
  check_matching_sizes("Binary function", "x", x, "y", y);
  check_matching_sizes("Binary function", "x", x, "z", z);
  using T_return = plain_type_t<decltype(apply_scalar_ternary(x[0], y[0], z[0], f))>;
  size_t y_size = y.size();
  std::vector<T_return> result(y_size);
  for (size_t i = 0; i < y_size; ++i) {
    result[i] = apply_scalar_ternary(x[i], y[i], z[i], f);
  }
  return result;
}


}  // namespace math
}  // namespace stan
#endif
