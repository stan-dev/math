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

// Forward declaration to allow specialisations
template <typename T1, typename T2, typename Enable = void,
          typename Enable2 = void>
struct apply_scalar_binary {};

/**
 * Base template class for vectorization of binary scalar functions
 * defined by applying a functor to a combination of scalars,
 * containers of matching sizes, or a combination of a scclar and a container.
 * These containers can be a standard library vector, Eigen dense
 * matrix expression template, or container of these. For each specialisation,
 * the same type as the largest (dimension) input is returned.
 *
 * This base template class takes (and returns) Eigen expression templates.
 */
template <typename T1, typename T2>
struct apply_scalar_binary<T1, T2, require_all_stan_scalar_t<T1, T2>> {
  /**
   * Member function for applying a functor to two scalars and subsequently
   * returning a scalar. 
   *
   * @tparam T1 Type of first argument to which functor is applied.
   * @tparam T2 Type of second argument to which functor is applied.
   * @tparam F Type of functor to apply.
   * @param x First input to which operation is applied.
   * @param y Second input to which operation is applied.
   * @param f functor to apply to inputs.
   * @return Scalar with result of applying functor to input.
   */
  template <typename F>
  static inline auto apply(const T1& x, const T2& y, const F& f) {
    return f(x,y);
  }
};

/**
 * Specialisation for use with two Eigen inputs. Eigen's binaryExpr framework
 * is used for more efficient indexing of both row- and column-major inputs
 * without separate loops.
 */
template <typename T1, typename T2>
struct apply_scalar_binary<T1, T2, require_all_eigen_t<T1, T2>> {
  /**
   * Member function for applying a functor to two Eigen objects and
   * subsequently returning an Eigen object. 
   *
   * @tparam T1 Type of first argument to which functor is applied.
   * @tparam T2 Type of second argument to which functor is applied.
   * @tparam F Type of functor to apply.
   * @param x First Eigen input to which operation is applied.
   * @param y Second Eigen input to which operation is applied.
   * @param f functor to apply to Eigen input.
   * @return Eigen object with result of applying functor to inputs.
   */
  template <typename F>
  static inline auto apply(const T1& x, const T2& y, const F& f) {
    return x.binaryExpr(y, f).eval();
  }
};

/**
 * Specialisation for use with one Eigen input and one scalar. Eigen's
 * unaryExpr framework is used for more efficient indexing of both row-
 * and column-major inputs. The unaryExpr framework also allows for the
 * scalar to be captured and applied to each element in the Eigen object.
 */
template <typename T1, typename T2>
struct apply_scalar_binary<T1, T2, require_eigen_t<T1>,
                           require_stan_scalar_t<T2>> {
  /**
   * Member function for applying a functor to an Eigen object and
   * a scalar, and subsequently returning an Eigen object. 
   *
   * @tparam T1 Type of Eigen object to which functor is applied.
   * @tparam T2 Type of scalar to which functor is applied.
   * @tparam F Type of functor to apply.
   * @param x Eigen input to which operation is applied.
   * @param y Scalar input to which operation is applied.
   * @param f functor to apply to Eigen and scalar inputs.
   * @return Eigen object with result of applying functor to inputs.
   */
  template <typename F>
  static inline auto apply(const T1& x, const T2& y, const F& f) {
    return x.unaryExpr([&f,&y](const auto& v){ return f(v, y); }).eval();
  }
};

/**
 * Specialisation for use with one scalar and one Eigen input. Eigen's
 * unaryExpr framework is used for more efficient indexing of both row-
 * and column-major inputs. The unaryExpr framework also allows for the
 * scalar to be captured and applied to each element in the Eigen object.
 */
template <typename T1, typename T2>
struct apply_scalar_binary<T1, T2, require_stan_scalar_t<T1>,
                           require_eigen_t<T2>> {
  /**
   * Member function for applying a functor to a scalar and
   * an Eigen object, and subsequently returning an Eigen object. 
   *
   * @tparam T1 Type of scalar to which functor is applied.
   * @tparam T2 Type of Eigen object to which functor is applied.
   * @tparam F Type of functor to apply.
   * @param x Scalar input to which operation is applied.
   * @param y Eigen input to which operation is applied.
   * @param f functor to apply to scalar and Eigen inputs.
   * @return Eigen object with result of applying functor to inputs.
   */
  template <typename F>
  static inline auto apply(const T1& x, const T2& y, const F& f) {
    return y.unaryExpr([&f,&x](const auto& v){ return f(x, v); }).eval();
  }
};

/**
 * Specialisation for use with (non-nested) std::vectors. Inputs are mapped
 * to Eigen column vectors and then the result is evaluated directly into the
 * returned std::vector (via Eigen::Map).
 *
 * The returned scalar type is deduced to allow for cases where the input and
 * return scalar types differ (e.g., functions implicitly promoting
 * integers).
 */
template <typename T1, typename T2>
struct apply_scalar_binary<T1, T2,
                           require_all_std_vector_vt<is_stan_scalar, T1, T2>> {

  template <typename F>
  static inline auto apply(const T1& x, const T2& y, const F& f) {
    decltype(auto) x_vec = as_column_vector_or_scalar(x);
    decltype(auto) y_vec = as_column_vector_or_scalar(y);
    using T_return = value_type_t<decltype(x_vec.binaryExpr(y_vec, f))>;
    std::vector<T_return> result(x.size());
    Eigen::Map<Eigen::Matrix<T_return, -1, 1>>(result.data(), result.size())
        = x_vec.binaryExpr(y_vec, f);
    return result;
  }

};

template <typename T1, typename T2>
struct apply_scalar_binary<T1, T2, require_std_vector_vt<is_stan_scalar, T1>,
                           require_stan_scalar_t<T2>> {

  template <typename F>
  static inline auto apply(const T1& x, const T2& y, const F& f) {
    decltype(auto) x_vec = as_column_vector_or_scalar(x);
    using T_return = value_type_t<decltype(f(x[0], y))>;
    std::vector<T_return> result(x.size());
    Eigen::Map<Eigen::Matrix<T_return, -1, 1>>(result.data(), result.size())
        = x_vec.unaryExpr([&f,&y](const auto& v){ return f(v, y); });
    return result;
  }

};

template <typename T1, typename T2>
struct apply_scalar_binary<T1, T2, require_stan_scalar_t<T1>,
                           require_std_vector_vt<is_stan_scalar, T2>> {

  template <typename F>
  static inline auto apply(const T1& x, const T2& y, const F& f) {
    decltype(auto) y_vec = as_column_vector_or_scalar(y);
    using T_return = value_type_t<decltype(f(x, y[0]))>;
    std::vector<T_return> result(y.size());
    Eigen::Map<Eigen::Matrix<T_return, -1, 1>>(result.data(), result.size())
        = y_vec.unaryExpr([&f,&x](const auto& v){ return f(x, v); });
    return result;
  }

};

template <typename T1, typename T2>
struct apply_scalar_binary<T1, T2,
                           require_all_std_vector_vt<is_container, T1, T2>> {
  using T1_vt = value_type_t<T1>;
  using T2_vt = value_type_t<T2>;

  template <typename F>
  static inline auto apply(const T1& x, const T2& y, const F& f) {
    using T_return
           = decltype(apply_scalar_binary<T1_vt, T2_vt>::apply(x[0], y[0], f));
    size_t y_size = y.size();
    std::vector<T_return> result(y_size);
    for (size_t i = 0; i < y_size; ++i)
      result[i] = apply_scalar_binary<T1_vt, T2_vt>::apply(x[i], y[i], f);
    return result;
  }

};

template <typename T1, typename T2>
struct apply_scalar_binary<T1, T2, require_std_vector_vt<is_container, T1>,
                           require_stan_scalar_t<T2>> {
  using T1_vt = value_type_t<T1>;

  template <typename F>
  static inline auto apply(const T1& x, const T2& y, const F& f) {
    using T_return
            = decltype(apply_scalar_binary<T1_vt, T2>::apply(x[0], y, f));
    size_t x_size = x.size();
    std::vector<T_return> result(x_size);
    for (size_t i = 0; i < x_size; ++i)
      result[i] = apply_scalar_binary<T1_vt, T2>::apply(x[i], y, f);
    return result;
  }

};

template <typename T1, typename T2>
struct apply_scalar_binary<T1, T2, require_stan_scalar_t<T1>,
                           require_std_vector_vt<is_container, T2>> {
  using T2_vt = value_type_t<T2>;

  template <typename F>
  static inline auto apply(const T1& x, const T2& y, const F& f) {
    using T_return
            = decltype(apply_scalar_binary<T1, T2_vt>::apply(x, y[0], f));
    size_t y_size = y.size();
    std::vector<T_return> result(y_size);
    for (size_t i = 0; i < y_size; ++i)
      result[i] = apply_scalar_binary<T1, T2_vt>::apply(x, y[i], f);
    return result;
  }

};

}  // namespace math
}  // namespace stan
#endif
