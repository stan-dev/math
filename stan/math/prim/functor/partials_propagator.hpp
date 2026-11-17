#ifndef STAN_MATH_PRIM_FUNCTOR_PARTIALS_PROPOGATOR_HPP
#define STAN_MATH_PRIM_FUNCTOR_PARTIALS_PROPOGATOR_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/err/check_size_match.hpp>
#include <stan/math/prim/meta/is_autodiff.hpp>
#include <stan/math/prim/meta/is_stan_scalar.hpp>
#include <stan/math/prim/meta/is_var_matrix.hpp>
#include <stan/math/prim/meta/is_vector.hpp>
#include <stan/math/prim/meta/value_type.hpp>
#include <stan/math/prim/meta/require_generics.hpp>
#include <stan/math/prim/meta/return_type.hpp>
#include <stan/math/prim/meta/plain_type.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>
#include <vector>
#include <type_traits>
#include <tuple>

namespace stan {
namespace math {

namespace internal {

template <typename ReturnType, typename Enable, typename... Ops>
class partials_propagator;

/** \ingroup type_trait
 * \callergraph
 * This template builds partial derivatives with respect to a
 * set of
 * operands. There are two reason for the generality of this
 * class. The first is to handle vector and scalar arguments
 * without needing to write additional code. The second is to use
 * this class for writing probability distributions that handle
 * primitives, reverse mode, and forward mode variables
 * seamlessly.
 *
 * Conceptually, this class is used when we want to manually calculate
 * the derivative of a function and store this manual result on the
 * autodiff stack in a sort of "compressed" form. Think of it like an
 * easy-to-use interface to rev/core/precomputed_gradients.
 *
 * This class supports nested container ("multivariate") use-cases
 * as well by exposing a partials_vec_ member on edges of the
 * appropriate type.
 *
 * This base template is instantiated when all operands are
 * primitives and we don't want to calculate derivatives at
 * all. So all Op1 - Op5 must be arithmetic primitives
 * like int or double. This is controlled with the
 * T_return_type type parameter.
 *
 * @tparam Ops type of the operands
 * @tparam ReturnType The type returned from calling the `build()` method.
 */
template <typename ReturnType, typename... Ops>
class partials_propagator<ReturnType, require_arithmetic_t<ReturnType>,
                          Ops...> {
 public:
  template <typename... Types>
  explicit partials_propagator(Types&&... /* ops */) noexcept {}

  /** \ingroup type_trait
   * Build the node to be stored on the autodiff graph.
   * This should contain both the value and the tangent.
   *
   * For scalars (this implementation), we don't calculate any derivatives.
   * For reverse mode, we end up returning a type of var that will calculate
   * the appropriate adjoint using the stored operands and partials.
   * Forward mode just calculates the tangent on the spot and returns it in
   * a vanilla fvar.
   *
   * @param value the return value of the function we are compressing
   * @return the value with its derivative
   */
  inline static constexpr double build(double value) noexcept { return value; }
};

}  // namespace internal

/**
 * Access the edge of an `partials_propagator`
 * @tparam I The index of the edge to access
 * @tparam Types The types inside of the operands and partials.
 * @param x An `partials_propagator` type whose edge will be accessed.
 */
template <std::size_t I, class... Types>
inline constexpr auto& edge(
    internal::partials_propagator<Types...>& x) noexcept {
  return std::get<I>(x.edges_);
};

/**
 * Access the partials for an edge of an `partials_propagator`
 * @tparam I The index of the edge to access
 * @tparam Types The types inside of the operands and partials.
 * @param x An `partials_propagator` type whose edge will be accessed.
 */
template <std::size_t I, class... Types>
inline constexpr auto& partials(
    internal::partials_propagator<Types...>& x) noexcept {
  return std::get<I>(x.edges_).partials_;
};

/**
 * Access the partials_vec for an edge of a `partials_propagator`
 * @tparam I The index of the edge to access
 * @tparam Types The types inside of the operands and partials.
 * @param x An `partials_propagator` type whose edge will be accessed.
 */
template <std::size_t I, class... Types>
inline constexpr auto& partials_vec(
    internal::partials_propagator<Types...>& x) noexcept {
  return std::get<I>(x.edges_).partials_vec_;
};

namespace internal {
/**
 * True if `T` is a scalar, an Eigen vector, a `std::vector` of scalars or a
 * `var_value` of an Eigen vector.
 */
template <typename T, typename = void>
struct is_scalar_or_vector_operand
    : bool_constant<is_stan_scalar<std::decay_t<T>>::value
                    || is_eigen_vector<std::decay_t<T>>::value> {};
template <typename T>
struct is_scalar_or_vector_operand<T, require_std_vector_t<T>>
    : is_stan_scalar<value_type_t<T>> {};
template <typename T>
struct is_scalar_or_vector_operand<T, require_var_matrix_t<T>>
    : is_eigen_vector<typename std::decay_t<T>::value_type> {};
}  // namespace internal

/**
 * Write values into the partials of an edge of a `partials_propagator`. Use
 * it to set the partials of a scalar, vector or data argument with one call
 * and without new memory.
 *
 * For an edge of a vector operand, the values are written in place into the
 * partials that the edge holds, through a column array view. For an edge of a
 * scalar operand, the sum of the elements of `x` is written.
 *
 * @tparam I index of the edge; its operand is a scalar, an Eigen vector, a
 * `std::vector` of scalars or a `var_value` of an Eigen vector
 * @tparam T type of the values: an Eigen column or row vector or array, or an
 * expression of one; a scalar is accepted only for an edge of a scalar operand
 * @tparam Types types inside the `partials_propagator`
 * @param ops the `partials_propagator`
 * @param x the values; for a vector operand, of the size of the operand. `x`
 * must not read the partials of edge `I` in another element order.
 * @throw std::invalid_argument if the operand is a vector and the size of `x`
 * differs from its size
 */
template <std::size_t I, typename T, class... Types>
inline void assign_partials(internal::partials_propagator<Types...>& ops,
                            const T& x) {
  // Types are the return type, the enable type, then the operand types.
  using T_op = std::tuple_element_t<I + 2, std::tuple<Types...>>;
  static_assert(internal::is_scalar_or_vector_operand<T_op>::value,
                "assign_partials: the operand must be a scalar or a vector");
  if constexpr (is_autodiff_v<T_op>) {
    auto& p = partials<I>(ops);
    if constexpr (is_stan_scalar<std::decay_t<T_op>>::value) {
      p = x;
    } else {
      static_assert(is_eigen_vector<std::decay_t<T>>::value,
                    "assign_partials: the values for a vector operand must "
                    "be an Eigen vector or array");
      check_size_match("assign_partials", "size of the values", x.size(),
                       "size of the partials", p.size());
      as_array_or_scalar(as_column_vector_or_scalar(p))
          = as_array_or_scalar(as_column_vector_or_scalar(x));
    }
  }
}

/**
 * Construct an `partials_propagator`.
 * @tparam Ops The type of the operands used in the edges of the
 * `partials_propagator`.
 * @param ops The operands to be placed into the edges.
 */
template <typename... Ops>
inline auto make_partials_propagator(Ops&&... ops) {
  using return_type = return_type_t<Ops...>;
  return internal::partials_propagator<return_type, void,
                                       plain_type_t<std::decay_t<Ops>>...>(
      std::forward<Ops>(ops)...);
}

}  // namespace math
}  // namespace stan
#endif
