#ifndef STAN_MATH_PRIM_META_OPERANDS_AND_PARTIALS_HPP
#define STAN_MATH_PRIM_META_OPERANDS_AND_PARTIALS_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/require_generics.hpp>
#include <stan/math/prim/meta/return_type.hpp>
#include <stan/math/prim/functor/broadcast_array.hpp>
#include <vector>
#include <type_traits>
#include <tuple>

namespace stan {
namespace math {

namespace internal {

template <typename ReturnType, typename Enable, typename... Ops>
struct operands_and_partials_impl;
/** \ingroup type_trait
 * \callergraph
 * An edge holds both the operands and its associated
 * partial derivatives. They're held together in the
 * same class because then we can keep the templating logic that
 * specializes on type of operand in one place.
 *
 * This is the base template class that ends up getting instantiated
 * for arithmetic primitives (doubles and ints).
 *
 * NB: since ops_partials_edge.partials_ and ops_partials_edge.partials_vec
 * are sometimes represented internally as a broadcast_array, we need to take
 * care with assignments to them. Indeed, we can assign any right hand side
 * which allows for indexing to a broadcast_array. The resulting behaviour is
 * that the entry for the first index is what gets assigned. The most common
 * use-case should be where the rhs is some container of length 1.
 *
 * @tparam ViewElt the type we expect to be at partials_[i]
 * @tparam Op the type of the operand
 */
template <typename ViewElt, typename Op, typename Enable = void>
class ops_partials_edge;

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
class operands_and_partials_impl<ReturnType, require_arithmetic_t<ReturnType>,
                                 Ops...> {
 public:
  template <typename... Types>
  explicit operands_and_partials_impl(Types&&... /* ops */) noexcept {}

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
  inline double build(double value) const noexcept { return value; }
  std::tuple<internal::ops_partials_edge<double, std::decay_t<Ops>>...> edges_;
};

/**
 * Class representing an edge with an inner type of double. This class
 *  should never be used by the program and only exists so that
 *  developer can write functions using `operands_and_partials` that works for
 *  double, vars, and fvar types.
 * @tparam ViewElt One of `double`, `var`, `fvar`.
 * @tparam Op The type of the input operand. It's scalar type
 *  for this specialization must be a `Arithmetic`
 */
template <typename ViewElt, typename Op>
struct ops_partials_edge<ViewElt, Op, require_st_arithmetic<Op>> {
  using inner_op = std::conditional_t<is_eigen<value_type_t<Op>>::value,
                                      value_type_t<Op>, Op>;
  using partials_t = empty_broadcast_array<ViewElt, inner_op>;
  /**
   * The `partials_` are always called in `if` statements that will be
   *  removed by the dead code elimination pass of the compiler. So if we ever
   *  move up to C++17 these can be made into `constexpr if` and
   *  this can be deleted.
   */
  partials_t partials_;
  empty_broadcast_array<partials_t, inner_op> partials_vec_;
  static constexpr double operands_{0};
  ops_partials_edge() {}
  template <typename T>
  explicit ops_partials_edge(T&& /* op */) noexcept {}

  /**
   * Get the operand for the edge. For doubles this is a compile time
   * expression returning zero.
   */
  static constexpr auto operand() noexcept { return static_cast<double>(0.0); }

  /**
   * Get the partial for the edge. For doubles this is a compile time
   * expression returning zero.
   */
  static constexpr auto partial() noexcept { return static_cast<double>(0.0); }
  /**
   * Return the tangent for the edge. For doubles this is a comple time
   * expression returning zero.
   */
  static constexpr double dx() { return static_cast<double>(0); }
  /**
   * Return the size of the operand for the edge. For doubles this is a comple
   * time expression returning zero.
   */
  static constexpr int size() { return 0; }  // reverse mode
};

template <typename ViewElt, typename Op>
constexpr double
    ops_partials_edge<ViewElt, Op, require_st_arithmetic<Op>>::operands_;

}  // namespace internal

/**
 * Construct an `operands_and_partials_impl`.
 * @tparam Ops The type of the operands used in the edges of the
 * `operands_and_partials_impl`.
 * @param ops The operands to be placed into the edges.
 */
template <typename... Ops>
inline auto operands_and_partials(Ops&&... ops) {
  using return_type = return_type_t<Ops...>;
  return internal::operands_and_partials_impl<return_type, void,
                                              std::decay_t<Ops>...>(
      std::forward<Ops>(ops)...);
}

/**
 * Access the edge of an `operands_and_partials_impl`
 * @tparam I The index of the edge to access
 * @tparam Types The types inside of the operands and partials.
 * @param x An `operands_and_partials_impl` type whose edge will be accessed.
 */
template <std::size_t I, class... Types>
inline constexpr auto& edge(
    internal::operands_and_partials_impl<Types...>& x) noexcept {
  return std::get<I>(x.edges_);
};

}  // namespace math
}  // namespace stan
#endif
