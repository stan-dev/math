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
 * matrix expression template, or container of these. For each specialisation,
 * the same vector type as the input is returned.
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
 * An std::vector (or scalar) is then returned as the result.
 */
template <typename T>
struct apply_vector_unary<T, require_std_vector_vt<is_stan_scalar, T>> {
  using T_scalar = scalar_type_t<T>;
  using T_map = typename Eigen::Map<const Eigen::Matrix<T_scalar, -1, 1>>;

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
  static inline std::vector<T_scalar> apply(const T& x, const F& f) {
    Eigen::Matrix<T_scalar, -1, 1> result
        = apply_vector_unary<T_map>::apply(as_column_vector_or_scalar(x), f);
    return std::vector<T_scalar>(result.data(), result.data() + result.size());
  }

  /**
   * Member function for applying a functor to a vector and a scalar
   * and subsequently returning a vector.
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
  static inline std::vector<T_scalar> apply_scalar(const T& x, const T2& y,
                                                   const F& f) {
    Eigen::Matrix<T_scalar, -1, 1> result
        = apply_vector_unary<T_map>::apply_scalar(as_column_vector_or_scalar(x),
                                                  y, f);
    return std::vector<T_scalar>(result.data(), result.data() + result.size());
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
  static inline T_scalar reduce(const T& x, const F& f) {
    return apply_vector_unary<T_map>::reduce(as_column_vector_or_scalar(x), f);
  }
};

/**
 * Specialisation for use with nested containers (std::vectors) of Eigen types.
 * For each of the member functions, an std::vector with the appropriate
 * type (vector or scalar) is returned.
 *
 */
template <typename T>
struct apply_vector_unary<T, require_std_vector_vt<is_eigen, T>> {
  using T_eigen = value_type_t<T>;
  using T_scalar = typename T_eigen::Scalar;
  using T_return
      = std::vector<Eigen::Matrix<T_scalar, T_eigen::RowsAtCompileTime,
                                  T_eigen::ColsAtCompileTime>>;

  /**
   * Member function for applying a functor to each Eigen type in an std::vector
   * and subsequently returning an std::vector of evaluated Eigen expressions.
   *
   * @tparam T Type of argument to which functor is applied.
   * @tparam F Type of functor to apply.
   * @param x std::vector of Eigen inputs to which operation is applied.
   * @param f functor to apply to vector input.
   * @return std::vector of Eigen objects with result of applying functor to
   *         input.
   */
  template <typename F>
  static inline T_return apply(const T& x, const F& f) {
    size_t x_size = x.size();
    T_return result(x_size);
    for (size_t i = 0; i < x_size; ++i)
      result[i] = apply_vector_unary<T_eigen>::apply(x[i], f);
    return result;
  }

  /**
   * Member function for applying a functor to each Eigen type in an
   * std::vector, as well as a scalar and subsequently returning an
   * std::vector of evaluated Eigen expressions. The provided scalar can either
   * be a single scalar which is applied to every Eigen object, or a vector of
   * scalars to be applied, one for each Eigen object in the std::vector
   *
   * @tparam T Type of argument to which functor is applied.
   * @tparam T2 Type of scalar (or vector) to pass to functor.
   * @tparam F Type of functor to apply.
   * @param x std::vector of Eigen inputs to which operation is applied.
   * @param y scalar (or vector) passed to functor.
   * @param f functor to apply to vector input.
   * @return std::vector of Eigen objects with result of applying functor and
   *         scalar to input.
   */
  template <typename F, typename T2>
  static inline T_return apply_scalar(const T& x, const T2& y, const F& f) {
    scalar_seq_view<T2> y_vec(y);
    size_t x_size = x.size();
    T_return result(x_size);
    for (size_t i = 0; i < x_size; ++i)
      result[i] = apply_vector_unary<T_eigen>::apply_scalar(x[i], y_vec[i], f);
    return result;
  }

  /**
   * Member function for applying a functor to each Eigen type in an
   * std::vector and subsequently returning an std::vector of scalars.
   *
   * @tparam T Type of argument to which functor is applied.
   * @tparam F Type of functor to apply.
   * @param x std::vector of Eigen inputs to which operation is applied.
   * @param f functor to apply to vector input.
   * @return std::vector of scalars with result of applying functor to input.
   */
  template <typename F>
  static inline std::vector<T_scalar> reduce(const T& x, const F& f) {
    size_t x_size = x.size();
    std::vector<T_scalar> result(x_size);
    for (size_t i = 0; i < x_size; ++i)
      result[i] = apply_vector_unary<T_eigen>::reduce(x[i], f);
    return result;
  }
};

/**
 * Specialisation for use with nested containers (std::vectors) of std::vectors.
 * For each of the member functions, an std::vector with the appropriate
 * type (vector or scalar) is returned.
 *
 */
template <typename T>
struct apply_vector_unary<T, require_std_vector_vt<is_std_vector, T>> {
  using T_scalar = scalar_type_t<T>;
  using T_return = typename std::vector<std::vector<T_scalar>>;
  using T_map = typename Eigen::Map<const Eigen::Matrix<T_scalar, -1, 1>>;

  /**
   * Member function for applying a functor to each std::vector in an
   * std::vector and subsequently returning an std::vector of std::vectors.
   *
   * @tparam T Type of argument to which functor is applied.
   * @tparam F Type of functor to apply.
   * @param x std::vector of std::vector to which operation is applied.
   * @param f functor to apply to vector input.
   * @return std::vector of std::vectors with result of applying functor to
   *         input.
   */
  template <typename F>
  static inline T_return apply(const T& x, const F& f) {
    size_t x_size = x.size();
    T_return result(x_size);
    Eigen::Matrix<T_scalar, -1, 1> inter;
    for (size_t i = 0; i < x_size; ++i) {
      inter = apply_vector_unary<T_map>::apply(as_column_vector_or_scalar(x[i]),
                                               f);
      result[i]
          = std::vector<T_scalar>(inter.data(), inter.data() + inter.size());
    }
    return result;
  }

  /**
   * Member function for applying a functor to each std::vector in an
   * std::vector, as well as a scalar and subsequently returning an
   * std::vector of evaluated std::vectors. The provided scalar can either
   * be a single scalar which is applied to every std::vector, or a vector of
   * scalars to be applied, one for each std::vector in the std::vector
   *
   * @tparam T Type of argument to which functor is applied.
   * @tparam T2 Type of scalar (or vector) to pass to functor.
   * @tparam F Type of functor to apply.
   * @param x std::vector of std::vectors to which operation is applied.
   * @param y scalar (or vector) passed to functor.
   * @param f functor to apply to vector input.
   * @return std::vector of std::vectors with result of applying functor and
   *         scalar to input.
   */
  template <typename F, typename T2>
  static inline T_return apply_scalar(const T& x, const T2& y, const F& f) {
    scalar_seq_view<T2> y_vec(y);
    size_t x_size = x.size();
    T_return result(x_size);
    Eigen::Matrix<T_scalar, -1, 1> inter;
    for (size_t i = 0; i < x_size; ++i) {
      inter = apply_vector_unary<T_map>::apply_scalar(
          as_column_vector_or_scalar(x[i]), y_vec[i], f);
      result[i]
          = std::vector<T_scalar>(inter.data(), inter.data() + inter.size());
    }
    return result;
  }

  /**
   * Member function for applying a functor to each std::vector in an
   * std::vector and subsequently returning an std::vector of scalars.
   *
   * @tparam T Type of argument to which functor is applied.
   * @tparam F Type of functor to apply.
   * @param x std::vector of std::vectors to which operation is applied.
   * @param f functor to apply to vector input.
   * @return std::vector of scalars with result of applying functor to input.
   */
  template <typename F>
  static inline std::vector<T_scalar> reduce(const T& x, const F& f) {
    size_t x_size = x.size();
    std::vector<T_scalar> result(x_size);
    for (size_t i = 0; i < x_size; ++i) {
      result[i] = apply_vector_unary<T_map>::reduce(
          as_column_vector_or_scalar(x[i]), f);
    }
    return result;
  }
};

}  // namespace math
}  // namespace stan
#endif
