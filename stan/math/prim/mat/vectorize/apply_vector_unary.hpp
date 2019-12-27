#ifndef STAN_MATH_PRIM_MAT_VECTORIZE_APPLY_VECTOR_UNARY_HPP
#define STAN_MATH_PRIM_MAT_VECTORIZE_APPLY_VECTOR_UNARY_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <vector>

namespace stan {
namespace math {

// Forward declaration to allow specialisations
template <typename T, typename Enable = void>
struct apply_vector_unary {};

/**
 * Base template class for vectorization of unary vector functions
 * defined by applying a functor to a standard library vector, Eigen dense
 * matrix expression template, or container of these.
 *
 * Three taxonomies of unary vector functions are implemented:
 * - f(vector) -> vector
 * - f(vector) -> scalar
 * - f(vector, scalar) -> vector
 *
 * This base template class takes (and returns) Eigen expression templates.
 */
template <typename T>
struct apply_vector_unary<T, require_eigen_t<T>> {
  /**
   * Member function for applying a functor to a vector and subsequently
   * returning a vector. The 'auto' return type is used here so that an
   * expression template is returned.
   *
   * @tparam T Type of argument to which functor is applied.
   * @tparam F Type of functor to apply.
   * @param x Eigen input to which operation is applied.
   * @param f functor to apply to Eigen input.
   * @return Eigen expression template with result of applying functor
   *         to input
   */
  template <typename F>
  static inline auto apply(const T& x, const F& f) {
    return f(x);
  }

  /**
   * Member function for applying a functor to a vector and a scalar
   * and subsequently returning a vector. The 'auto' return type is
   * used here so that an expression template is returned.
   *
   * @tparam T Type of argument to which functor is applied.
   * @tparam T2 Type of scalar to pass to functor.
   * @tparam F Type of functor to apply.
   * @param x Eigen input to which operation is applied.
   * @param y scalar passed to functor.
   * @param f functor to apply to Eigen input.
   * @return Eigen expression template with result of applying functor
   *         and scalar to input
   */
  template <typename F, typename T2>
  static inline auto apply_scalar(const T& x, const T2& y, const F& f) {
    return f(x, y);
  }

  /**
   * Member function for applying a functor to a vector and subsequently
   * returning a scalar.
   *
   * @tparam T Type of argument to which functor is applied.
   * @tparam F Type of functor to apply.
   * @param x Eigen input to which operation is applied.
   * @param f functor to apply to Eigen input.
   * @return scalar result of applying functor to input.
   */
  template <typename F>
  static inline auto reduce(const T& x, const F& f) {
    return f(x);
  }
};

/**
 * Specialisation for use with (non-nested) std::vectors. Inputs are mapped
 * to Eigen column vectors and then passed to the base (Eigen) template.
 */
template <typename T>
struct apply_vector_unary<T, require_std_vector_vt<is_stan_scalar, T>> {
  using scalar_t = scalar_type_t<T>;
  using map_t = typename Eigen::Map<const Eigen::Matrix<scalar_t, -1, 1>>;

  /**
   * Member function for applying a functor to a vector and subsequently
   * returning a vector. 
   *
   * @tparam T Type of argument to which functor is applied.
   * @tparam F Type of functor to apply.
   * @param x std::vector input to which operation is applied.
   * @param f functor to apply to vector input.
   * @return std::vector with result of applying functor to input.
   */
  template <typename F>
  static inline std::vector<scalar_t> apply(const T& x, const F& f) {
    Eigen::Matrix<scalar_t, -1, 1> result
          = apply_vector_unary<map_t>::apply(as_column_vector_or_scalar(x), f);
    return std::vector<scalar_t>(result.data(),
                                 result.data() + result.size());
  }

  /**
   * Member function for applying a functor to a vector and a scalar
   * and subsequently returning a vector. The 'auto' return type is
   * used here so that an expression template is returned.
   *
   * @tparam T Type of argument to which functor is applied.
   * @tparam T2 Type of scalar to pass to functor.
   * @tparam F Type of functor to apply.
   * @param x std::vector input to which operation is applied.
   * @param y scalar passed to functor.
   * @param f functor to apply to vector input.
   * @return std::vector with result of applying functor and scalar to input.
   */
  template <typename F, typename T2>
  static inline std::vector<scalar_t> apply_scalar(const T& x,
                                                   const T2& y,
                                                   const F& f) {
    Eigen::Matrix<scalar_t, -1, 1> result
       = apply_vector_unary<map_t>::apply_scalar(as_column_vector_or_scalar(x),
                                                 y, f);
    return std::vector<scalar_t>(result.data(),
                                 result.data() + result.size());
  }

  /**
   * Member function for applying a functor to a vector and subsequently
   * returning a scalar.
   *
   * @tparam T Type of argument to which functor is applied.
   * @tparam F Type of functor to apply.
   * @param x Eigen input to which operation is applied.
   * @param f functor to apply to std::vector input.
   * @return scalar result of applying functor to input vector.
   */
  template <typename F>
  static inline scalar_t reduce(const T& x, const F& f) {
    return apply_vector_unary<map_t>::reduce(as_column_vector_or_scalar(x), f);
  }
};

/**
 * Specialisation for use with std::vectors of Eigen types
 */
template <typename T>
struct apply_vector_unary<T, require_std_vector_vt<is_eigen, T>> {
  using eigen_t = value_type_t<T>;
  using scalar_t = typename eigen_t::Scalar;
  using return_t = std::vector<Eigen::Matrix<scalar_t,
                                              eigen_t::RowsAtCompileTime,
                                              eigen_t::ColsAtCompileTime>>;

  template <typename F>
  static inline return_t apply(const T& x, const F& f) {
    size_t x_size = x.size();
    return_t result(x_size);
    for (size_t i = 0; i < x_size; ++i)
      result[i] = apply_vector_unary<eigen_t>::apply(x[i], f);
    return result;
  }

  template <typename F, typename T2>
  static inline return_t apply_scalar(const T& x, const T2& y, const F& f) {
    scalar_seq_view<T2> y_vec(y);
    size_t x_size = x.size();
    return_t result(x_size);
    for (size_t i = 0; i < x_size; ++i)
      result[i] = apply_vector_unary<eigen_t>::apply_scalar(x[i], y_vec[i], f);
    return result;
  }

  template <typename F>
  static inline std::vector<scalar_t> reduce(const T& x, const F& f) {
    size_t x_size = x.size();
    std::vector<scalar_t> result(x_size);
    for (size_t i = 0; i < x_size; ++i)
      result[i] = apply_vector_unary<eigen_t>::reduce(x[i], f);
    return result;
  }
};

template <typename T>
struct apply_vector_unary<T, require_std_vector_vt<is_std_vector, T>> {
  using scalar_t = scalar_type_t<T>;
  using return_t = typename std::vector<std::vector<scalar_t>>;
  using map_t =
    typename Eigen::Map<const Eigen::Matrix<scalar_t, -1, 1>>;

  template <typename F>
  static inline return_t apply(const T& x, const F& f) {
    size_t x_size = x.size();
    return_t result(x_size);
    Eigen::Matrix<scalar_t, -1, 1> inter;
    for (size_t i = 0; i < x_size; ++i) {
      inter = apply_vector_unary<map_t>::apply(as_column_vector_or_scalar(x[i]),
                                               f);
      result[i] = std::vector<scalar_t>(inter.data(),
                                        inter.data() + inter.size());
    }
    return result;
  }

  template <typename F, typename T2>
  static inline return_t apply_scalar(const T& x, const T2& y, const F& f) {
    scalar_seq_view<T2> y_vec(y);
    size_t x_size = x.size();
    return_t result(x_size);
    Eigen::Matrix<scalar_t, -1, 1> inter;
    for (size_t i = 0; i < x_size; ++i) {
      inter = apply_vector_unary<map_t>::apply_scalar(
                    as_column_vector_or_scalar(x[i]), y_vec[i], f);
      result[i] = std::vector<scalar_t>(inter.data(),
                                        inter.data() + inter.size());
    }
    return result;
  }

  template <typename F>
  static inline std::vector<scalar_t> reduce(const T& x, const F& f) {
    size_t x_size = x.size();
    std::vector<scalar_t> result(x_size);
    for (size_t i = 0; i < x_size; ++i) {
      result[i] = apply_vector_unary<map_t>::reduce(
                    as_column_vector_or_scalar(x[i]), f);
    }
    return result;
  }
};

}  // namespace math
}  // namespace stan
#endif
