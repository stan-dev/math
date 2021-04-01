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
template <typename Op1 = double, typename Op2 = double, typename Op3 = double,
          typename Op4 = double, typename Op5 = double,
          typename T_return_type = return_type_t<Op1, Op2, Op3, Op4, Op5>>
class operands_and_partials;  // Forward declaration

namespace internal {
template <typename ViewElt, typename Op, typename = void>
struct ops_partials_edge;
/**
 * Class representing an edge with an inner type of double. This class
 *  should never be used by the program and only exists so that
 *  developer can write functions using `operands_and_partials` that works for
 *  double, vars, and fvar types.
 * @tparam ViewElt One of `double`, `var`, `fvar`.
 * @tparam Op The type of the input operand. It's scalar type
 *  for this specialization must be an `Arithmetic`
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
  static constexpr double operand() noexcept { return 0.0; }

  /**
   * Get the partial for the edge. For doubles this is a compile time
   * expression returning zero.
   */
  static constexpr double partial() noexcept { return 0.0; }
  /**
   * Return the tangent for the edge. For doubles this is a compile time
   * expression returning zero.
   */
  static constexpr double dx() noexcept { return 0.0; }
  /**
   * Return the size of the operand for the edge. For doubles this is a compile
   * time expression returning zero.
   */
  static constexpr int size() noexcept { return 0; }  // reverse mode

 private:
  template <typename, typename, typename, typename, typename, typename>
  friend class stan::math::operands_and_partials;
};
template <typename ViewElt, typename Op>
constexpr double
    ops_partials_edge<ViewElt, Op, require_st_arithmetic<Op>>::operands_;
}  // namespace internal

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
 * @tparam Op1 type of the first operand
 * @tparam Op2 type of the second operand
 * @tparam Op3 type of the third operand
 * @tparam Op4 type of the fourth operand
 * @tparam Op5 type of the fifth operand
 * @tparam T_return_type return type of the expression. This defaults
 *   to calling a template metaprogram that calculates the scalar
 *   promotion of Op1..Op4
 */
template <typename Op1, typename Op2, typename Op3, typename Op4, typename Op5,
          typename T_return_type>
class operands_and_partials {
 public:
  explicit operands_and_partials(const Op1& /* op1 */) noexcept {}
  operands_and_partials(const Op1& /* op1 */, const Op2& /* op2 */) noexcept {}
  operands_and_partials(const Op1& /* op1 */, const Op2& /* op2 */,
                        const Op3& /* op3 */) noexcept {}
  operands_and_partials(const Op1& /* op1 */, const Op2& /* op2 */,
                        const Op3& /* op3 */, const Op4& /* op4 */) noexcept {}
  operands_and_partials(const Op1& /* op1 */, const Op2& /* op2 */,
                        const Op3& /* op3 */, const Op4& /* op4 */,
                        const Op5& /* op5 */) noexcept {}

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

  // These will always be 0 size base template instantiations (above).
  internal::ops_partials_edge<double, std::decay_t<Op1>> edge1_;
  internal::ops_partials_edge<double, std::decay_t<Op2>> edge2_;
  internal::ops_partials_edge<double, std::decay_t<Op3>> edge3_;
  internal::ops_partials_edge<double, std::decay_t<Op4>> edge4_;
  internal::ops_partials_edge<double, std::decay_t<Op5>> edge5_;
};

}  // namespace math
}  // namespace stan
#endif
