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
 * @tparam T1 Type of first argument to which functor is applied.
 * @tparam T2 Type of second argument to which functor is applied.
 * @tparam F Type of functor to apply.
 * @param x First input to which operation is applied.
 * @param y Second input to which operation is applied.
 * @param f functor to apply to inputs.
 * @return Scalar with result of applying functor to input.
 */
template <typename T1, typename T2, typename T3, typename F,
          require_all_stan_scalar_t<T1, T2, T3>* = nullptr>
inline auto apply_scalar_ternary(const T1& x, const T2& y, const T3& z, const F& f) {
  return f(x, y, z);
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

// matrix, matrix, matrix
template <typename T1, typename T2, typename T3, typename F,
          require_all_eigen_t<T1, T2, T3>* = nullptr>
inline auto apply_scalar_ternary(const T1& x, const T2& y, const T3& z, const F& f) {
  return x.ternaryExpr(y, z, f);
}

// std::vector<real>, std::vector<real>, std::vector<real>
template <typename T1, typename T2, typename T3, typename F,
          require_all_std_vector_vt<is_stan_scalar, T1, T2, T3>* = nullptr>
inline auto apply_scalar_ternary(const T1& x, const T2& y, const T3& z, const F& f) {
  decltype(auto) x_vec = as_column_vector_or_scalar(x);
  decltype(auto) y_vec = as_column_vector_or_scalar(y);
  decltype(auto) z_vec = as_column_vector_or_scalar(z);
  using T_return = std::decay_t<decltype(f(x[0], y[0], z[0]))>;
  std::vector<T_return> result(x.size());
  Eigen::Map<Eigen::Matrix<T_return, -1, 1>>(result.data(), result.size())
      = x_vec.ternaryExpr(y_vec, z_vec, f);
  return result;
}

// matrix, matrix, std::vector<std::vector<int>>
template <typename T1, typename T2, typename T3, typename F,
          require_all_eigen_matrix_dynamic_vt<is_stan_scalar, T1, T2>* = nullptr,
          require_std_vector_vt<is_std_vector, T3>* = nullptr,
          require_std_vector_st<std::is_integral, T3>* = nullptr>
inline auto apply_scalar_ternary(const T1& x, const T2& y, const T3& z, const F& f) {
  if (num_elements(x) != num_elements(y) || num_elements(z) != num_elements(y)) {
    std::ostringstream msg;
    msg << "Inputs to vectorized ternary function must match in"
        << " size if one is not a scalar";
    throw std::invalid_argument(msg.str());
  }
  using return_st = decltype(f(x(0), y(0), z[0][0]));
  Eigen::Matrix<return_st, Eigen::Dynamic, Eigen::Dynamic> result(x.rows(),
                                                                  x.cols());
  for (size_t i = 0; i < z.size(); ++i) {
    result.row(i) = apply_scalar_ternary(x.row(i).transpose(),
                                         y.row(i).transpose(),
                                         as_column_vector_or_scalar(z[i]), f);
  }
  return result;
}

// matrix, std::vector<std::vector<int>>, matrix
template <typename T1, typename T2, typename T3, typename F,
          require_all_eigen_matrix_dynamic_vt<is_stan_scalar, T1, T3>* = nullptr,
          require_std_vector_vt<is_std_vector, T2>* = nullptr,
          require_std_vector_st<std::is_integral, T2>* = nullptr>
inline auto apply_scalar_ternary(const T1& x, const T2& y, const T3& z, const F& f) {
  if (num_elements(x) != num_elements(y) || num_elements(z) != num_elements(y)) {
    std::ostringstream msg;
    msg << "Inputs to vectorized ternary function must match in"
        << " size if one is not a scalar";
    throw std::invalid_argument(msg.str());
  }
  using return_st = decltype(f(x(0), y[0][0], z(0)));
  Eigen::Matrix<return_st, Eigen::Dynamic, Eigen::Dynamic> result(x.rows(),
                                                                  x.cols());
  for (size_t i = 0; i < z.size(); ++i) {
    result.row(i) = apply_scalar_ternary(x.row(i).transpose(),
                                         as_column_vector_or_scalar(y),
                                         z.row(i).transpose(), f);
  }
  return result;
}

// matrix, std::vector<std::vector<int>>, std::vector<std::vector<int>>
template <typename T1, typename T2, typename T3, typename F,
          require_eigen_matrix_dynamic_vt<is_stan_scalar, T1>* = nullptr,
          require_all_std_vector_vt<is_std_vector, T2, T3>* = nullptr,
          require_all_std_vector_st<std::is_integral, T2, T3>* = nullptr>
inline auto apply_scalar_ternary(const T1& x, const T2& y, const T3& z, const F& f) {
  if (num_elements(x) != num_elements(y) || num_elements(z) != num_elements(y)) {
    std::ostringstream msg;
    msg << "Inputs to vectorized ternary function must match in"
        << " size if one is not a scalar";
    throw std::invalid_argument(msg.str());
  }
  using return_st = decltype(f(x(0), y[0][0], z[0][0]));
  Eigen::Matrix<return_st, Eigen::Dynamic, Eigen::Dynamic> result(x.rows(),
                                                                  x.cols());
  for (size_t i = 0; i < z.size(); ++i) {
    result.row(i) = apply_scalar_ternary(x.row(i).transpose(),
                                         as_column_vector_or_scalar(y),
                                         as_column_vector_or_scalar(z), f);
  }
  return result;
}

// std::vector<std::vector<int>>, matrix, std::vector<std::vector<int>>
template <typename T1, typename T2, typename T3, typename F,
          require_eigen_matrix_dynamic_vt<is_stan_scalar, T2>* = nullptr,
          require_all_std_vector_vt<is_std_vector, T1, T3>* = nullptr,
          require_all_std_vector_st<std::is_integral, T1, T3>* = nullptr>
inline auto apply_scalar_ternary(const T1& x, const T2& y, const T3& z, const F& f) {
  if (num_elements(x) != num_elements(y) || num_elements(z) != num_elements(y)) {
    std::ostringstream msg;
    msg << "Inputs to vectorized ternary function must match in"
        << " size if one is not a scalar";
    throw std::invalid_argument(msg.str());
  }
  using return_st = decltype(f(x[0][0], y(0), z[0][0]));
  Eigen::Matrix<return_st, Eigen::Dynamic, Eigen::Dynamic> result(y.rows(),
                                                                  y.cols());
  for (size_t i = 0; i < z.size(); ++i) {
    result.row(i) = apply_scalar_ternary(as_column_vector_or_scalar(x),
                                         y.row(i).transpose(),
                                         as_column_vector_or_scalar(z), f);
  }
  return result;
}

// std::vector<std::vector<int>>, std::vector<std::vector<int>>, matrix
template <typename T1, typename T2, typename T3, typename F,
          require_eigen_matrix_dynamic_vt<is_stan_scalar, T3>* = nullptr,
          require_all_std_vector_vt<is_std_vector, T1, T2>* = nullptr,
          require_all_std_vector_st<std::is_integral, T1, T2>* = nullptr>
inline auto apply_scalar_ternary(const T1& x, const T2& y, const T3& z, const F& f) {
  if (num_elements(x) != num_elements(y) || num_elements(z) != num_elements(y)) {
    std::ostringstream msg;
    msg << "Inputs to vectorized ternary function must match in"
        << " size if one is not a scalar";
    throw std::invalid_argument(msg.str());
  }
  using return_st = decltype(f(x[0][0], y[0][0], z(0)));
  Eigen::Matrix<return_st, Eigen::Dynamic, Eigen::Dynamic> result(z.rows(),
                                                                  z.cols());
  for (size_t i = 0; i < z.size(); ++i) {
    result.row(i) = apply_scalar_ternary(as_column_vector_or_scalar(x),
                                         as_column_vector_or_scalar(y),
                                         z.row(i).transpose(), f);
  }
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
  }).eval();
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
  return apply_scalar_binary(y, z, [&](const auto& b, const auto& c) { return f(x, b, c); }).eval();
}

// T3 Container
template <typename T1, typename T2, typename T3, typename F,
          require_container_t<T3>* = nullptr,
          require_all_stan_scalar_t<T1, T2>* = nullptr>
inline auto apply_scalar_ternary(const T1& x, const T2& y, const T3& z, const F& f) {
  return apply_scalar_binary(y, z, [&](const auto& b, const auto& c) { return f(x, b, c); }).eval();
}
// Nested containers
template <typename T1, typename T2, typename T3, typename F,
          require_all_std_vector_vt<is_container_or_var_matrix, T1, T2, T3>* = nullptr>
inline auto apply_scalar_ternary(const T1& x, const T2& y, const T3& z, const F& f) {
  using T_return = plain_type_t<decltype(apply_scalar_ternary(x[0], y[0], z[0], f))>;
  size_t y_size = y.size();
  std::vector<T_return> result(y_size);
  for (size_t i = 0; i < y_size; ++i) {
    result[i] = apply_scalar_ternary(x[i], y[i], z[i], f);
  }
  return result;
}

template <typename T1, typename T2, typename T3, typename F,
          require_all_std_vector_vt<is_container_or_var_matrix, T1, T2>* = nullptr,
          require_stan_scalar_t<T3>* = nullptr>
inline auto apply_scalar_ternary(const T1& x, const T2& y, const T3& z, const F& f) {
  using T_return = plain_type_t<decltype(apply_scalar_ternary(x[0], y[0], z, f))>;
  size_t y_size = y.size();
  std::vector<T_return> result(y_size);
  for (size_t i = 0; i < y_size; ++i) {
    result[i] = apply_scalar_ternary(x[i], y[i], z, f);
  }
  return result;
}

template <typename T1, typename T2, typename T3, typename F,
          require_all_std_vector_vt<is_container_or_var_matrix, T1, T3>* = nullptr,
          require_stan_scalar_t<T2>* = nullptr>
inline auto apply_scalar_ternary(const T1& x, const T2& y, const T3& z, const F& f) {
  using T_return = plain_type_t<decltype(apply_scalar_ternary(x[0], y, z[0], f))>;
  size_t y_size = y.size();
  std::vector<T_return> result(y_size);
  for (size_t i = 0; i < y_size; ++i) {
    result[i] = apply_scalar_ternary(x[i], y, z[i], f);
  }
  return result;
}

template <typename T1, typename T2, typename T3, typename F,
          require_std_vector_vt<is_container_or_var_matrix, T1>* = nullptr,
          require_all_stan_scalar_t<T2, T3>* = nullptr>
inline auto apply_scalar_ternary(const T1& x, const T2& y, const T3& z, const F& f) {
  using T_return = plain_type_t<decltype(apply_scalar_ternary(x[0], y, z, f))>;
  size_t y_size = y.size();
  std::vector<T_return> result(y_size);
  for (size_t i = 0; i < y_size; ++i) {
    result[i] = apply_scalar_ternary(x[i], y, z, f);
  }
  return result;
}

template <typename T1, typename T2, typename T3, typename F,
          require_all_std_vector_vt<is_container_or_var_matrix, T2, T3>* = nullptr,
          require_stan_scalar_t<T1>* = nullptr>
inline auto apply_scalar_ternary(const T1& x, const T2& y, const T3& z, const F& f) {
  using T_return = plain_type_t<decltype(apply_scalar_ternary(x, y[0], z[0], f))>;
  size_t y_size = y.size();
  std::vector<T_return> result(y_size);
  for (size_t i = 0; i < y_size; ++i) {
    result[i] = apply_scalar_ternary(x, y[i], z[i], f);
  }
  return result;
}
template <typename T1, typename T2, typename T3, typename F,
          require_std_vector_vt<is_container_or_var_matrix, T2>* = nullptr,
          require_all_stan_scalar_t<T1, T3>* = nullptr>
inline auto apply_scalar_ternary(const T1& x, const T2& y, const T3& z, const F& f) {
  using T_return = plain_type_t<decltype(apply_scalar_ternary(x, y[0], z, f))>;
  size_t y_size = y.size();
  std::vector<T_return> result(y_size);
  for (size_t i = 0; i < y_size; ++i) {
    result[i] = apply_scalar_ternary(x, y[i], z, f);
  }
  return result;
}
template <typename T1, typename T2, typename T3, typename F,
          require_std_vector_vt<is_container_or_var_matrix, T3>* = nullptr,
          require_all_stan_scalar_t<T1, T2>* = nullptr>
inline auto apply_scalar_ternary(const T1& x, const T2& y, const T3& z, const F& f) {
  using T_return = plain_type_t<decltype(apply_scalar_ternary(x, y, z[0], f))>;
  size_t y_size = y.size();
  std::vector<T_return> result(y_size);
  for (size_t i = 0; i < y_size; ++i) {
    result[i] = apply_scalar_ternary(x, y, z[i], f);
  }
  return result;
}

}  // namespace math
}  // namespace stan
#endif
