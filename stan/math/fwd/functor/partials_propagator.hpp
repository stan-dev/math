#ifndef STAN_MATH_FWD_META_PARTIALS_PROPOGATOR_HPP
#define STAN_MATH_FWD_META_PARTIALS_PROPOGATOR_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/fwd/core/fvar.hpp>
#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/functor/operands_and_partials.hpp>
#include <stan/math/prim/functor/broadcast_array.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <stan/math/prim/functor/apply.hpp>
#include <stan/math/prim/functor/for_each.hpp>
#include <vector>

namespace stan {
namespace math {
namespace internal {

/**
 * This class builds partial derivatives with respect to a set of
 * operands. There are two reason for the generality of this
 * class. The first is to handle vector and scalar arguments
 * without needing to write additional code. The second is to use
 * this class for writing probability distributions that handle
 * primitives, reverse mode, and forward mode variables
 * seamlessly.
 *
 * Conceptually, this class is used when we want to calculate manually
 * the derivative of a function and store this manual result on the
 * autodiff stack in a sort of "compressed" form. Think of it like an
 * easy-to-use interface to rev/core/precomputed_gradients.
 *
 * This class now supports multivariate use-cases as well by
 * exposing edge#_.partials_vec
 *
 * This is the specialization for when the return type is fvar,
 * which should be for forward mode and all higher-order cases.
 *
 * NB: since ops_partials_edge.partials_ and ops_partials_edge.partials_vec
 * are sometimes represented internally as a broadcast_array, we need to take
 * care with assignments to them. Indeed, we can assign any right hand side
 * which allows for indexing to a broadcast_array. The resulting behaviour is
 * that the entry for the first index is what gets assigned. The most common
 * use-case should be where the rhs is some container of length 1.
 *
 * @tparam Ops The types included in the edges.
 * @tparam ReturnType return type of the expression. This defaults
 *   to a template metaprogram that calculates the scalar promotion of
 *   Ops
 */
template <typename ReturnType, typename... Ops>
class partials_propagator<ReturnType, require_fvar_t<ReturnType>, Ops...> {
 public:
  using Dx = partials_type_t<ReturnType>;
  std::tuple<
      internal::ops_partials_edge<Dx, plain_type_t<std::decay_t<Ops>>>...>
      edges_;
  using T_return_type = fvar<Dx>;

  template <typename... Types>
  explicit partials_propagator(Types&&... ops)
      : edges_(internal::ops_partials_edge<Dx, plain_type_t<std::decay_t<Ops>>>(
          std::forward<Types>(ops))...) {}

  /** \ingroup type_trait
   * Build the node to be stored on the autodiff graph.
   * This should contain both the value and the tangent.
   *
   * For scalars, we don't calculate any tangents.
   * For reverse mode, we end up returning a type of var that will calculate
   * the appropriate adjoint using the stored operands and partials.
   * Forward mode just calculates the tangent on the spot and returns it in
   * a vanilla fvar.
   *
   * @param value the return value of the function we are compressing
   * @return the value with its derivative
   */
  inline T_return_type build(Dx value) {
    auto deriv = stan::math::apply(
        [](auto&... args) { return sum_dx(args...); }, edges_);
    return T_return_type(value, deriv);
  }
};

}  // namespace internal
}  // namespace math
}  // namespace stan
#endif
