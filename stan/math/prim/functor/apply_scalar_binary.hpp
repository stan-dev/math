#ifndef STAN_MATH_PRIM_FUNCTOR_APPLY_SCALAR_BINARY_HPP
#define STAN_MATH_PRIM_FUNCTOR_APPLY_SCALAR_BINARY_HPP

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
 * matrix expression template, or container of these. For each specialization,
 * the same type as the largest (dimension) input is returned.
 *
 * This base template function takes (and returns) scalars.
 *
 * @tparam F Type of functor to apply.
 * @tparam T1 Type of first argument to which functor is applied.
 * @tparam T2 Type of second argument to which functor is applied.
 * @param f functor to apply to inputs.
 * @param x First input to which operation is applied.
 * @param y Second input to which operation is applied.
 * @return Scalar with result of applying functor to input.
 */
template <typename F, typename T1, typename T2,
          require_all_stan_scalar_t<T1, T2>* = nullptr>
inline auto apply_scalar_binary(F&& f, T1&& x, T2&& y) {
  return std::forward<F>(f)(std::forward<T1>(x), std::forward<T2>(y));
}

/**
 * Specialization for use with two Eigen inputs. Eigen's binaryExpr framework
 * is used for more efficient indexing of both row- and column-major inputs
 * without separate loops.
 *
 * @tparam F Type of functor to apply.
 * @tparam T1 Type of first argument to which functor is applied.
 * @tparam T2 Type of second argument to which functor is applied.
 * @param f functor to apply to Eigen input.
 * @param x First Eigen input to which operation is applied.
 * @param y Second Eigen input to which operation is applied.
 * @return Eigen object with result of applying functor to inputs.
 */
template <typename F, typename T1, typename T2,
          require_all_eigen_t<T1, T2>* = nullptr>
inline auto apply_scalar_binary(F&& f, T1&& x, T2&& y) {
  check_matching_dims("Binary function", "x", x, "y", y);
  return make_holder(
      [](auto&& f_inner, auto&& x_inner, auto&& y_inner) {
        return std::forward<decltype(x_inner)>(x_inner).binaryExpr(
            std::forward<decltype(y_inner)>(y_inner),
            std::forward<decltype(f_inner)>(f_inner));
      },
      std::forward<F>(f), std::forward<T1>(x), std::forward<T2>(y));
}

/**
 * Specialization for use with one Eigen vector (row or column) and
 * a one-dimensional std::vector of integer types
 *
 * @tparam F Type of functor to apply.
 * @tparam T1 Type of first argument to which functor is applied.
 * @tparam T2 Type of second argument to which functor is applied.
 * @param f functor to apply to inputs.
 * @param x Eigen input to which operation is applied.
 * @param y Integer std::vector input to which operation is applied.
 * @return Eigen object with result of applying functor to inputs.
 */
template <typename F, typename T1, typename T2,
          require_eigen_vector_vt<is_stan_scalar, T1>* = nullptr,
          require_std_vector_vt<std::is_integral, T2>* = nullptr>
inline auto apply_scalar_binary(F&& f, T1&& x, T2&& y) {
  check_matching_sizes("Binary function", "x", x, "y", y);
  return make_holder(
      [](auto&& f_inner, auto&& x_inner, auto&& y_inner) {
        using int_vec_t = promote_scalar_t<value_type_t<decltype(y_inner)>,
                                           plain_type_t<decltype(x_inner)>>;
        auto y_map = make_holder(
            [](auto&& y_inner_) {
              return Eigen::Map<const int_vec_t>(y_inner_.data(),
                                                 y_inner_.size());
            },
            std::forward<decltype(y_inner)>(y_inner));
        return std::forward<decltype(x_inner)>(x_inner).binaryExpr(
            y_map, std::forward<decltype(f_inner)>(f_inner));
      },
      std::forward<F>(f), std::forward<T1>(x), std::forward<T2>(y));
}

/**
 * Specialization for use with a one-dimensional std::vector of integer types
 * and one Eigen vector (row or column).
 *
 * @tparam F Type of functor to apply.
 * @tparam T1 Type of first argument to which functor is applied.
 * @tparam T2 Type of second argument to which functor is applied.
 * @param f functor to apply to inputs.
 * @param x Integer std::vector input to which operation is applied.
 * @param y Eigen input to which operation is applied.
 * @return Eigen object with result of applying functor to inputs.
 */
template <typename F, typename T1, typename T2,
          require_std_vector_vt<std::is_integral, T1>* = nullptr,
          require_eigen_vector_vt<is_stan_scalar, T2>* = nullptr>
inline auto apply_scalar_binary(F&& f, T1&& x, T2&& y) {
  check_matching_sizes("Binary function", "x", x, "y", y);
  return make_holder(
      [](auto&& f_inner, auto&& x_inner, auto&& y_inner) {
        using int_vec_t = promote_scalar_t<value_type_t<decltype(x_inner)>,
                                           plain_type_t<decltype(y_inner)>>;
        auto x_map = make_holder(
            [](auto&& x_inner_) {
              return Eigen::Map<const int_vec_t>(x_inner_.data(),
                                                 x_inner_.size());
            },
            std::forward<decltype(x_inner)>(x_inner));
        return x_map.binaryExpr(std::forward<decltype(y_inner)>(y_inner),
                                std::forward<decltype(f_inner)>(f_inner));
      },
      std::forward<F>(f), std::forward<T1>(x), std::forward<T2>(y));
}

/**
 * Specialization for use with one Eigen matrix and
 * a two-dimensional std::vector of integer types
 *
 * @tparam F Type of functor to apply.
 * @tparam T1 Type of first argument to which functor is applied.
 * @tparam T2 Type of second argument to which functor is applied.
 * @param f functor to apply to inputs.
 * @param x Eigen matrix input to which operation is applied.
 * @param y Nested integer std::vector input to which operation is applied.
 * @return Eigen object with result of applying functor to inputs.
 */
template <typename F, typename T1, typename T2,
          require_eigen_matrix_dynamic_vt<is_stan_scalar, T1>* = nullptr,
          require_std_vector_vt<is_std_vector, T2>* = nullptr,
          require_std_vector_st<std::is_integral, T2>* = nullptr>
inline auto apply_scalar_binary(F&& f, T1&& x, T2&& y) {
  if (num_elements(x) != num_elements(y)) {
    std::ostringstream msg;
    msg << "Inputs to vectorized binary function must match in"
        << " size if one is not a scalar";
    throw std::invalid_argument(msg.str());
  }
  using return_st = decltype(f(x(0), y[0][0]));
  Eigen::Matrix<return_st, Eigen::Dynamic, Eigen::Dynamic> result(x.rows(),
                                                                  x.cols());
  for (size_t i = 0; i < y.size(); ++i) {
    result.row(i) = apply_scalar_binary(f, x.row(i).transpose(),
                                        as_column_vector_or_scalar(y[i]));
  }
  return result;
}

/**
 * Specialization for use with a two-dimensional std::vector of integer types
 * and one Eigen matrix.
 *
 * @tparam F Type of functor to apply.
 * @tparam T1 Type of first argument to which functor is applied.
 * @tparam T2 Type of second argument to which functor is applied.
 * @param f functor to apply to inputs.
 * @param x Nested integer std::vector input to which operation is applied.
 * @param y Eigen matrix input to which operation is applied.
 * @return Eigen object with result of applying functor to inputs.
 */
template <typename F, typename T1, typename T2,
          require_std_vector_vt<is_std_vector, T1>* = nullptr,
          require_std_vector_st<std::is_integral, T1>* = nullptr,
          require_eigen_matrix_dynamic_vt<is_stan_scalar, T2>* = nullptr>
inline auto apply_scalar_binary(F&& f, T1&& x, T2&& y) {
  if (num_elements(x) != num_elements(y)) {
    std::ostringstream msg;
    msg << "Inputs to vectorized binary function must match in"
        << " size if one is not a scalar";
    throw std::invalid_argument(msg.str());
  }
  using return_st = decltype(f(x[0][0], y(0)));
  Eigen::Matrix<return_st, Eigen::Dynamic, Eigen::Dynamic> result(y.rows(),
                                                                  y.cols());
  for (size_t i = 0; i < x.size(); ++i) {
    result.row(i) = apply_scalar_binary(f, as_column_vector_or_scalar(x[i]),
                                        y.row(i).transpose());
  }
  return result;
}

/**
 * Specialization for use when the first input is an Eigen type and the second
 * is a scalar. Eigen's unaryExpr framework is used for more efficient indexing
 * of both row- and column-major inputs. The unaryExpr framework also allows
 * for the scalar to be captured and applied to each element in the Eigen
 * object.
 *
 * @tparam F Type of functor to apply.
 * @tparam T1 Type of Eigen object to which functor is applied.
 * @tparam T2 Type of scalar to which functor is applied.
 * @param f functor to apply to Eigen and scalar inputs.
 * @param x Eigen input to which operation is applied.
 * @param y Scalar input to which operation is applied.
 * @return Eigen object with result of applying functor to inputs.
 */
template <typename F, typename T1, typename T2, require_eigen_t<T1>* = nullptr,
          require_stan_scalar_t<T2>* = nullptr>
inline auto apply_scalar_binary(F&& f, T1&& x, T2&& y) {
  return make_holder(
      [](auto&& f_inner, auto&& x_inner, auto&& y_inner) {
        return std::forward<decltype(x_inner)>(x_inner).unaryExpr(
            [f_inner_ = std::forward<decltype(f_inner)>(f_inner),
             y_inner](auto&& v) { return f_inner_(v, y_inner); });
      },
      std::forward<F>(f), std::forward<T1>(x), std::forward<T2>(y));
}

/**
 * Specialization for use when the first input is an scalar and the second is
 * an Eigen type. Eigen's unaryExpr framework is used for more efficient
 * indexing of both row- and column-major inputs. The unaryExpr framework also
 * allows for the scalar to be captured and applied to each element in the
 * Eigen object.
 *
 * @tparam F Type of functor to apply.
 * @tparam T1 Type of scalar to which functor is applied.
 * @tparam T2 Type of Eigen object to which functor is applied.
 * @param f functor to apply to Eigen and scalar inputs.
 * @param x Scalar input to which operation is applied.
 * @param y Eigen input to which operation is applied.
 * @return Eigen object with result of applying functor to inputs.
 */
template <typename F, typename T1, typename T2,
          require_stan_scalar_t<T1>* = nullptr, require_eigen_t<T2>* = nullptr>
inline auto apply_scalar_binary(F&& f, T1&& x, T2&& y) {
  return make_holder(
      [](auto&& f_inner, auto&& x_inner, auto&& y_inner) {
        return std::forward<decltype(y_inner)>(y_inner).unaryExpr(
            [f_inner_ = std::forward<decltype(f_inner)>(f_inner),
             x_inner](auto&& v) {
              return f_inner_(x_inner, std::forward<decltype(v)>(v));
            });
      },
      std::forward<F>(f), std::forward<T1>(x), std::forward<T2>(y));
}

/**
 * Specialization for use with (non-nested) std::vectors. Inputs are mapped
 * to Eigen column vectors and then the result is evaluated directly into the
 * returned std::vector (via Eigen::Map).
 *
 * The returned scalar type is deduced to allow for cases where the input and
 * return scalar types differ (e.g., functions implicitly promoting
 * integers).
 *
 * @tparam F Type of functor to apply.
 * @tparam T1 Type of first std::vector to which functor is applied.
 * @tparam T2 Type of second std::vector to which functor is applied.
 * @param f functor to apply to std::vector inputs.
 * @param x First std::vector input to which operation is applied.
 * @param y Second std::vector input to which operation is applied.
 * @return std::vector with result of applying functor to inputs.
 */
template <typename F, typename T1, typename T2,
          require_all_std_vector_vt<is_stan_scalar, T1, T2>* = nullptr>
inline auto apply_scalar_binary(F&& f, T1&& x, T2&& y) {
  check_matching_sizes("Binary function", "x", x, "y", y);
  using T_return = std::decay_t<decltype(f(x[0], y[0]))>;
  decltype(auto) x_vec = as_column_vector_or_scalar(std::forward<T1>(x));
  decltype(auto) y_vec = as_column_vector_or_scalar(std::forward<T2>(y));
  std::vector<T_return> result(x_vec.size());
  Eigen::Map<Eigen::Matrix<T_return, -1, 1>>(result.data(), result.size())
      = x_vec.binaryExpr(y_vec, std::forward<F>(f));
  return result;
}

/**
 * Specialization for use when the first input is a (non-nested) std::vector
 * and the second is a scalar. The std::vector input is mapped to an Eigen
 * column vector and then the result is evaluated directly into the returned
 * std::vector (via Eigen::Map).
 *
 * The returned scalar type is deduced to allow for cases where the input and
 * return scalar types differ (e.g., functions implicitly promoting
 * integers).
 *
 * @tparam F Type of functor to apply.
 * @tparam T1 Type of std::vector to which functor is applied.
 * @tparam T2 Type of scalar to which functor is applied.
 * @param f functor to apply to std::vector and scalar inputs.
 * @param x std::vector input to which operation is applied.
 * @param y Scalar input to which operation is applied.
 * @return std::vector with result of applying functor to inputs.
 */
template <typename F, typename T1, typename T2,
          require_std_vector_vt<is_stan_scalar, T1>* = nullptr,
          require_stan_scalar_t<T2>* = nullptr>
inline auto apply_scalar_binary(F&& f, T1&& x, T2&& y) {
  decltype(auto) x_vec = as_column_vector_or_scalar(std::forward<T1>(x));
  using T_return = std::decay_t<decltype(f(x[0], y))>;
  std::vector<T_return> result(x_vec.size());
  Eigen::Map<Eigen::Matrix<T_return, -1, 1>>(result.data(), result.size())
      = x_vec.unaryExpr(
          [f_ = std::forward<F>(f), y](auto&& v) { return f_(v, y); });
  return result;
}

/**
 * Specialization for use when the first input is a scalar and the second is a
 * (non-nested) std::vector. The std::vector input is mapped to an Eigen
 * column vector and then the result is evaluated directly into the returned
 * std::vector (via Eigen::Map).
 *
 * The returned scalar type is deduced to allow for cases where the input and
 * return scalar types differ (e.g., functions implicitly promoting
 * integers).
 *
 * @tparam F Type of functor to apply.
 * @tparam T1 Type of scalar to which functor is applied.
 * @tparam T2 Type of std::vector to which functor is applied.
 * @param f functor to apply to std::vector and scalar inputs.
 * @param x Scalar input to which operation is applied.
 * @param y std::vector input to which operation is applied.
 * @return std::vector with result of applying functor to inputs.
 */
template <typename F, typename T1, typename T2,
          require_stan_scalar_t<T1>* = nullptr,
          require_std_vector_vt<is_stan_scalar, T2>* = nullptr>
inline auto apply_scalar_binary(F&& f, T1&& x, T2&& y) {
  using T_return = std::decay_t<decltype(f(x, y[0]))>;
  decltype(auto) y_vec = as_column_vector_or_scalar(std::forward<T2>(y));
  std::vector<T_return> result(y_vec.size());
  Eigen::Map<Eigen::Matrix<T_return, -1, 1>>(result.data(), result.size())
      = y_vec.unaryExpr(
          [f_ = std::forward<F>(f), x](auto&& v) { return f_(x, v); });
  return result;
}

/**
 * Specialization for use with two nested containers (std::vectors).
 * The returned scalar type is deduced to allow for cases where the input and
 * return scalar types differ (e.g., functions implicitly promoting
 * integers).
 *
 * @tparam F Type of functor to apply.
 * @tparam T1 Type of first std::vector to which functor is applied.
 * @tparam T2 Type of second std::vector to which functor is applied.
 * @param f functor to apply to std::vector inputs.
 * @param x First std::vector input to which operation is applied.
 * @param y Second std::vector input to which operation is applied.
 * @return std::vector with result of applying functor to inputs.
 */
template <
    typename F, typename T1, typename T2,
    require_all_std_vector_vt<is_container_or_var_matrix, T1, T2>* = nullptr>
inline auto apply_scalar_binary(F&& f, T1&& x, T2&& y) {
  check_matching_sizes("Binary function", "x", x, "y", y);
  using T_return = plain_type_t<decltype(apply_scalar_binary(f, x[0], y[0]))>;
  size_t y_size = y.size();
  std::vector<T_return> result(y_size);
  for (size_t i = 0; i < y_size; ++i) {
    result[i] = apply_scalar_binary(f, x[i], y[i]);
  }
  return result;
}

/**
 * Specialization for use when the first input is a nested std::vector and the
 * second is a scalar. The returned scalar type is deduced to allow for cases
 * where the input and return scalar types differ (e.g., functions implicitly
 * promoting integers).
 *
 * @tparam F Type of functor to apply.
 * @tparam T1 Type of std::vector to which functor is applied.
 * @tparam T2 Type of scalar to which functor is applied.
 * @param f functor to apply to inputs.
 * @param x std::vector input to which operation is applied.
 * @param y Scalar input to which operation is applied.
 * @return std::vector with result of applying functor to inputs.
 */
template <typename F, typename T1, typename T2,
          require_std_vector_vt<is_container_or_var_matrix, T1>* = nullptr,
          require_stan_scalar_t<T2>* = nullptr>
inline auto apply_scalar_binary(F&& f, T1&& x, T2&& y) {
  using T_return = plain_type_t<decltype(apply_scalar_binary(f, x[0], y))>;
  size_t x_size = x.size();
  std::vector<T_return> result(x_size);
  for (size_t i = 0; i < x_size; ++i) {
    result[i] = apply_scalar_binary(f, x[i], y);
  }
  return result;
}

/**
 * Specialization for use when the first input is a scalar and the second is a
 * nested std::vector. The returned scalar type is deduced to allow for cases
 * where the input and return scalar types differ (e.g., functions implicitly
 * promoting integers).
 *
 * @tparam F Type of functor to apply.
 * @tparam T1 Type of scalar to which functor is applied.
 * @tparam T2 Type of std::vector to which functor is applied.
 * @param f functor to apply to inputs.
 * @param x Scalar input to which operation is applied.
 * @param y std::vector input to which operation is applied.
 * @return std::vector with result of applying functor to inputs.
 */
template <typename F, typename T1, typename T2,
          require_stan_scalar_t<T1>* = nullptr,
          require_std_vector_vt<is_container_or_var_matrix, T2>* = nullptr>
inline auto apply_scalar_binary(F&& f, T1&& x, T2&& y) {
  using T_return = plain_type_t<decltype(apply_scalar_binary(f, x, y[0]))>;
  size_t y_size = y.size();
  std::vector<T_return> result(y_size);
  for (size_t i = 0; i < y_size; ++i) {
    result[i] = apply_scalar_binary(f, x, y[i]);
  }
  return result;
}

}  // namespace math
}  // namespace stan
#endif
