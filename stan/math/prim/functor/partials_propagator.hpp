#ifndef STAN_MATH_PRIM_FUNCTOR_PARTIALS_PROPOGATOR_HPP
#define STAN_MATH_PRIM_FUNCTOR_PARTIALS_PROPOGATOR_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/require_generics.hpp>
#include <stan/math/prim/meta/return_type.hpp>
#include <stan/math/prim/functor/broadcast_array.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>
#include <vector>
#include <type_traits>
#include <tuple>

namespace stan {
namespace math {

namespace internal {

template <typename ReturnType, typename Enable, typename... Ops>
struct partials_propagator;

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
  inline static double build(double value) noexcept { return value; }
  std::tuple<internal::ops_partials_edge<double, std::decay_t<Ops>>...> edges_;
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
