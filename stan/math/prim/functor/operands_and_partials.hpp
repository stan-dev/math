#ifndef STAN_MATH_PRIM_FUNCTOR_OPERANDS_AND_PARTIALS_HPP
#define STAN_MATH_PRIM_FUNCTOR_OPERANDS_AND_PARTIALS_HPP

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
          typename Op4 = double, typename Op5 = double, typename Op6 = double,
          typename Op7 = double, typename Op8 = double,
          typename T_return_type
          = return_type_t<Op1, Op2, Op3, Op4, Op5, Op6, Op7, Op8>>
class operands_and_partials;  // Forward declaration

namespace internal {

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
class ops_partials_edge<ViewElt, Op, require_st_arithmetic<Op>> {
 public:
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
  ops_partials_edge() = default;

  template <typename T,
            require_not_same_t<
                std::decay_t<T>,
                std::decay_t<ops_partials_edge<
                    ViewElt, Op, require_st_arithmetic<Op>>>>* = nullptr>
  constexpr explicit ops_partials_edge(T&& /* op */) noexcept {}

  /**
   * Get the operand for the edge. For doubles this is a compile time
   * expression returning zero.
   */
  static constexpr double operand() noexcept {
    return static_cast<double>(0.0);
  }

  /**
   * Get the partial for the edge. For doubles this is a compile time
   * expression returning zero.
   */
  static constexpr double partial() noexcept {
    return static_cast<double>(0.0);
  }
  /**
   * Return the tangent for the edge. For doubles this is a comple time
   * expression returning zero.
   */
  static constexpr double dx() noexcept { return static_cast<double>(0); }
  /**
   * Return the size of the operand for the edge. For doubles this is a comple
   * time expression returning zero.
   */
  static constexpr int size() noexcept { return 0; }  // reverse mode
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
 * all. So all Op1 - Op8 must be arithmetic primitives
 * like int or double. This is controlled with the
 * T_return_type type parameter.
 *
 * @tparam Op1 type of the first operand
 * @tparam Op2 type of the second operand
 * @tparam Op3 type of the third operand
 * @tparam Op4 type of the fourth operand
 * @tparam Op5 type of the fifth operand
 * @tparam Op6 type of the sixth operand
 * @tparam Op7 type of the seventh operand
 * @tparam Op8 type of the eighth operand
 * @tparam T_return_type return type of the expression. This defaults
 *   to calling a template metaprogram that calculates the scalar
 *   promotion of Op1..Op8
 */
template <typename Op1, typename Op2, typename Op3, typename Op4, typename Op5,
          typename Op6, typename Op7, typename Op8, typename T_return_type>
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
  operands_and_partials(const Op1& /* op1 */, const Op2& /* op2 */,
                        const Op3& /* op3 */, const Op4& /* op4 */,
                        const Op5& /* op5 */, const Op6& /* op6 */) noexcept {}
  operands_and_partials(const Op1& /* op1 */, const Op2& /* op2 */,
                        const Op3& /* op3 */, const Op4& /* op4 */,
                        const Op5& /* op5 */, const Op6& /* op6 */,
                        const Op7& /* op7 */) noexcept {}
  operands_and_partials(const Op1& /* op1 */, const Op2& /* op2 */,
                        const Op3& /* op3 */, const Op4& /* op4 */,
                        const Op5& /* op5 */, const Op6& /* op6 */,
                        const Op7& /* op7 */, const Op8& /* op8 */) noexcept {}

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
  internal::ops_partials_edge<double, std::decay_t<Op6>> edge6_;
  internal::ops_partials_edge<double, std::decay_t<Op7>> edge7_;
  internal::ops_partials_edge<double, std::decay_t<Op8>> edge8_;
};

}  // namespace math
}  // namespace stan
#endif
