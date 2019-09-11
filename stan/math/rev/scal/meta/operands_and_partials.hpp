#ifndef STAN_MATH_REV_SCAL_META_OPERANDS_AND_PARTIALS_HPP
#define STAN_MATH_REV_SCAL_META_OPERANDS_AND_PARTIALS_HPP

#include <stan/math/rev/core/chainablestack.hpp>
#include <stan/math/rev/core/precomputed_gradients.hpp>
#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/core/vari.hpp>
#include <stan/math/rev/scal/meta/is_var.hpp>
#include <stan/math/prim/scal/meta/broadcast_array.hpp>
#include <stan/math/prim/scal/meta/operands_and_partials.hpp>
#include <stan/math/prim/scal/meta/is_vector_like.hpp>
#include <stan/math/prim/scal/meta/conjunction.hpp>
#include <type_traits>
#include <utility>

namespace stan {
namespace math {
namespace internal {
template <typename ViewElt, typename Op>
class ops_partials_edge<ViewElt, Op,
                        std::enable_if_t<std::is_floating_point<ViewElt>::value
                                         && is_var<Op>::value>> {
 public:
  double partial_;
  broadcast_array<ViewElt> partials_;
  template <typename OpC, typename = std::enable_if_t<std::is_same<std::decay_t<Op>, std::decay_t<OpC>>::value>>
  explicit ops_partials_edge(OpC&& op)
      : partial_(0), partials_(partial_), operand_(std::forward<OpC>(op)) {}

 private:
  template <typename, typename, typename, typename, typename, typename,
            typename>
  friend class stan::math::operands_and_partials;
  Op operand_;

  void dump_partials(double* partials) { *partials = this->partial_; }
  void dump_operands(vari** varis) { *varis = this->operand_.vi_; }
  int size() const { return 1; }
};
}  // namespace internal

/**
 * This class builds partial derivatives with respect to a set of
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
 * This class now supports multivariate use-cases as well by
 * exposing edge#_.partials_vec
 *
 * This is the specialization for when the return type is var,
 * which should be for all of the reverse mode cases.
 *
 * NB: since ops_partials_edge.partials_ and ops_partials_edge.partials_vec
 * are sometimes represented internally as a broadcast_array, we need to take
 * care with assignments to them. Indeed, we can assign any right hand side
 * which allows for indexing to a broadcast_array. The resulting behaviour is
 * that the entry for the first index is what gets assigned. The most common
 * use-case should be where the rhs is some container of length 1.
 *
 * @tparam Op1 type of the first operand
 * @tparam Op2 type of the second operand
 * @tparam Op3 type of the third operand
 * @tparam Op4 type of the fourth operand
 * @tparam Op5 type of the fifth operand
 */
template <typename Op1, typename Op2, typename Op3, typename Op4, typename Op5,
          typename T_return_type>
class operands_and_partials<Op1, Op2, Op3, Op4, Op5, T_return_type,
                            std::enable_if_t<is_var<T_return_type>::value>> {
 public:
  internal::ops_partials_edge<double, Op1> edge1_;
  internal::ops_partials_edge<double, Op2> edge2_;
  internal::ops_partials_edge<double, Op3> edge3_;
  internal::ops_partials_edge<double, Op4> edge4_;
  internal::ops_partials_edge<double, Op5> edge5_;

  template <typename T1, typename T2>
  using is_same_op = std::is_same<std::decay_t<T1>, std::decay_t<T2>>;

  template <typename OpC1,
   std::enable_if_t<is_same_op<Op1, OpC1>::value>...>
  explicit operands_and_partials(OpC1&& o1) : edge1_(std::forward<OpC1>(o1)) {}

  template <typename OpC1, typename OpC2,
   std::enable_if_t<conjunction<is_same_op<Op1, OpC1>, is_same_op<Op2, OpC2>>::value>...>
  operands_and_partials(OpC1&& o1, OpC2&& o2)
      : edge1_(std::forward<OpC1>(o1)), edge2_(std::forward<OpC2>(o2)) {}

  template <typename OpC1, typename OpC2, typename OpC3,
   std::enable_if_t<conjunction<is_same_op<Op1, OpC1>, is_same_op<Op2, OpC2>,
    is_same_op<Op3, OpC3>>::value>...>
  operands_and_partials(OpC1&& o1, OpC2&& o2, OpC3&& o3)
      : edge1_(std::forward<OpC1>(o1)), edge2_(std::forward<OpC2>(o2)), edge3_(std::forward<OpC3>(o3)) {}

  template <typename OpC1, typename OpC2, typename OpC3, typename OpC4,
   std::enable_if_t<conjunction<is_same_op<Op1, OpC1>, is_same_op<Op2, OpC2>,
    is_same_op<Op3, OpC3>, is_same_op<Op4, OpC4>>::value>...>
  operands_and_partials(OpC1&& o1, OpC2&& o2, OpC3&& o3, OpC4&& o4)
      : edge1_(std::forward<OpC1>(o1)), edge2_(std::forward<OpC2>(o2)),
        edge3_(std::forward<OpC3>(o3)), edge4_(std::forward<OpC4>(o4)) {}

  template <typename OpC1, typename OpC2, typename OpC3, typename OpC4, typename OpC5,
   std::enable_if_t<conjunction<is_same_op<Op1, OpC1>, is_same_op<Op2, OpC2>,
    is_same_op<Op3, OpC3>, is_same_op<Op4, OpC4>, is_same_op<Op5, OpC5>>::value>...>
  operands_and_partials(OpC1&& o1, OpC2&& o2, OpC3&& o3, OpC4&& o4, OpC5&& o5)
      : edge1_(std::forward<OpC1>(o1)), edge2_(std::forward<OpC2>(o2)),
        edge3_(std::forward<OpC3>(o3)), edge4_(std::forward<OpC4>(o4)),
        edge5_(std::forward<OpC5>(o5)) {}

  /**
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
   * @return the node to be stored in the expression graph for autodiff
   */
  T_return_type build(double value) {
    size_t size = edge1_.size() + edge2_.size() + edge3_.size() + edge4_.size()
                  + edge5_.size();
    vari** varis
        = ChainableStack::instance_->memalloc_.alloc_array<vari*>(size);
    double* partials
        = ChainableStack::instance_->memalloc_.alloc_array<double>(size);
    int idx = 0;
    edge1_.dump_operands(&varis[idx]);
    edge1_.dump_partials(&partials[idx]);
    edge2_.dump_operands(&varis[idx += edge1_.size()]);
    edge2_.dump_partials(&partials[idx]);
    edge3_.dump_operands(&varis[idx += edge2_.size()]);
    edge3_.dump_partials(&partials[idx]);
    edge4_.dump_operands(&varis[idx += edge3_.size()]);
    edge4_.dump_partials(&partials[idx]);
    edge5_.dump_operands(&varis[idx += edge4_.size()]);
    edge5_.dump_partials(&partials[idx]);

    return {new precomputed_gradients_vari(value, size, varis, partials)};
  }
};
}  // namespace math
}  // namespace stan
#endif
