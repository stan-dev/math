#ifndef STAN_MATH_PRIM_META_APPLY_VECTOR_UNARY_HPP
#define STAN_MATH_PRIM_META_APPLY_VECTOR_UNARY_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/meta/require_generics.hpp>
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
 * Two taxonomies of unary vector functions are implemented:
 * - f(vector) -> vector
 * - f(vector) -> scalar
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
   * Member function for applying a functor to a vector and subsequently
   * returning a scalar. The reduction to a scalar needs to be implemented
   * in the definition of the functor.
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
  using T_vt = value_type_t<T>;
  using T_map = typename Eigen::Map<const Eigen::Matrix<T_vt, -1, 1>>;

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
  static inline std::vector<T_vt> apply(const T& x, const F& f) {
    std::vector<T_vt> result(x.size());
    Eigen::Map<Eigen::Matrix<T_vt, -1, 1>>(result.data(), result.size())
        = apply_vector_unary<T_map>::apply(as_column_vector_or_scalar(x), f);
    return result;
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
  static inline T_vt reduce(const T& x, const F& f) {
    return apply_vector_unary<T_map>::reduce(as_column_vector_or_scalar(x), f);
  }
};

/**
 * Specialisation for use with nested containers (std::vectors).
 * For each of the member functions, an std::vector with the appropriate
 * type (vector or scalar) is returned.
 *
 */
template <typename T>
struct apply_vector_unary<T, require_std_vector_vt<is_container, T>> {
  using T_vt = value_type_t<T>;
  using T_st = value_type_t<T_vt>;

  /**
   * Member function for applying a functor to each container in an std::vector
   * and subsequently returning an std::vector of containers.
   *
   * @tparam T Type of argument to which functor is applied.
   * @tparam F Type of functor to apply.
   * @param x std::vector of containers to which operation is applied.
   * @param f functor to apply to vector input.
   * @return std::vector of containers with result of applying functor to
   *         input.
   */
  template <typename F>
  static inline std::vector<T_vt> apply(const T& x, const F& f) {
    size_t x_size = x.size();
    std::vector<T_vt> result(x_size);
    for (size_t i = 0; i < x_size; ++i)
      result[i] = apply_vector_unary<T_vt>::apply(x[i], f);
    return result;
  }

  /**
   * Member function for applying a functor to each container in an
   * std::vector and subsequently returning an std::vector of scalars.
   *
   * @tparam T Type of argument to which functor is applied.
   * @tparam F Type of functor to apply.
   * @param x std::vector of containers to which operation is applied.
   * @param f functor to apply to vector input.
   * @return std::vector of scalars with result of applying functor to input.
   */
  template <typename F>
  static inline std::vector<T_st> reduce(const T& x, const F& f) {
    size_t x_size = x.size();
    std::vector<T_st> result(x_size);
    for (size_t i = 0; i < x_size; ++i)
      result[i] = apply_vector_unary<T_vt>::reduce(x[i], f);
    return result;
  }
};

}  // namespace math
}  // namespace stan
#endif
