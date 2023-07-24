#ifndef STAN_MATH_PRIM_FUNCTOR_APPLY_VECTOR_UNARY_HPP
#define STAN_MATH_PRIM_FUNCTOR_APPLY_VECTOR_UNARY_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/meta/is_stan_scalar.hpp>
#include <stan/math/prim/meta/is_container.hpp>
#include <stan/math/prim/meta/is_eigen_matrix_base.hpp>
#include <stan/math/prim/meta/plain_type.hpp>
#include <stan/math/prim/meta/require_generics.hpp>
#include <vector>

namespace stan {
namespace math {
// Forward declaration to allow specializations
template <typename T, typename Enable = void>
struct apply_vector_unary {};

/**
 * Base template class for vectorization of unary vector functions
 * defined by applying a functor to a standard library vector, Eigen dense
 * matrix expression template, or container of these. For each specialization,
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
   * returning a vector.
   *
   * @tparam T Type of argument to which functor is applied.
   * @tparam F Type of functor to apply.
   * @param x Eigen input to which operation is applied.
   * @param f functor to apply to Eigen input.
   * @return Eigen object with result of applying functor to input
   */
  template <typename F, typename T2 = T,
            require_t<is_eigen_matrix_base<plain_type_t<T2>>>* = nullptr>
  static inline auto apply(const T& x, const F& f) {
    return make_holder([](const auto& a) { return a.matrix().derived(); },
                       f(x));
  }

  template <typename F, typename T2 = T,
            require_t<is_eigen_array<plain_type_t<T2>>>* = nullptr>
  static inline auto apply(const T& x, const F& f) {
    return make_holder([](const auto& a) { return a.array().derived(); }, f(x));
  }

  /**
   * Member function for applying a functor to a vector and subsequently
   * returning a vector. This is a variant of `apply` that does not construct
   * `holder`, so it is up to the caller to ensure the returned expression is
   * evaluated before `x` is destructed.
   *
   * @tparam T Type of argument to which functor is applied.
   * @tparam F Type of functor to apply.
   * @param x Eigen input to which operation is applied.
   * @param f functor to apply to Eigen input.
   * @return Eigen object with result of applying functor to input
   */
  template <typename F, typename T2 = T,
            require_t<is_eigen_matrix_base<plain_type_t<T2>>>* = nullptr>
  static inline auto apply_no_holder(const T& x, const F& f) {
    return f(x).matrix().derived();
  }

  template <typename F, typename T2 = T,
            require_t<is_eigen_array<plain_type_t<T2>>>* = nullptr>
  static inline auto apply_no_holder(const T& x, const F& f) {
    return f(x).array().derived();
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
 * Specialization for use with (non-nested) std::vectors. Inputs are mapped
 * to Eigen column vectors and then the result is evaluated directly into the
 * returned std::vector (via Eigen::Map).
 *
 * The returned scalar type is deduced to allow for cases where the input and
 * return scalar types differ (e.g., functions implicitly promoting
 * integers).
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
  static inline auto apply(const T& x, const F& f) {
    using T_return = value_type_t<decltype(f(as_column_vector_or_scalar(x)))>;
    std::vector<T_return> result(x.size());
    Eigen::Map<Eigen::Matrix<T_return, -1, 1>>(result.data(), result.size())
        = f(as_column_vector_or_scalar(x)).matrix();
    return result;
  }

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
  static inline auto apply_no_holder(const T& x, const F& f) {
    return apply(x, f);
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
  static inline auto reduce(const T& x, const F& f) {
    return apply_vector_unary<T_map>::reduce(as_column_vector_or_scalar(x), f);
  }
};

/**
 * Specialization for use with nested containers (std::vectors).
 * For each of the member functions, an std::vector with the appropriate
 * type (vector or scalar) is returned.
 *
 * The returned scalar type is deduced to allow for cases where the input and
 * return scalar types differ (e.g., functions implicitly promoting
 * integers).
 *
 */
template <typename T>
struct apply_vector_unary<
    T, require_std_vector_vt<is_container_or_var_matrix, T>> {
  using T_vt = value_type_t<T>;

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
  static inline auto apply(const T& x, const F& f) {
    using T_return
        = plain_type_t<decltype(apply_vector_unary<T_vt>::apply(x[0], f))>;
    std::vector<T_return> result(x.size());
    std::transform(x.begin(), x.end(), result.begin(), [&f](auto&& xx) {
      return apply_vector_unary<T_vt>::apply_no_holder(xx, f);
    });
    return result;
  }

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
  static inline auto apply_no_holder(const T& x, const F& f) {
    return apply(x, f);
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
  static inline auto reduce(const T& x, const F& f) {
    size_t x_size = x.size();
    using T_return = decltype(apply_vector_unary<T_vt>::reduce(x[0], f));
    std::vector<T_return> result(x_size);
    for (size_t i = 0; i < x_size; ++i)
      result[i] = apply_vector_unary<T_vt>::reduce(x[i], f);
    return result;
  }
};

}  // namespace math
}  // namespace stan
#endif
