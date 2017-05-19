#ifndef STAN_MATH_PRIM_SCAL_META_OPERANDS_AND_PARTIALS_HPP
#define STAN_MATH_PRIM_SCAL_META_OPERANDS_AND_PARTIALS_HPP

#include <stan/math/prim/scal/meta/broadcast_array.hpp>
#include <stan/math/prim/scal/meta/return_type.hpp>

namespace stan {
  namespace math {
    template <typename Op1 = double, typename Op2 = double,
              typename Op3 = double, typename Op4 = double,
              typename T_return_type
              = typename return_type<Op1, Op2, Op3, Op4>::type>
    class operands_and_partials;

    namespace internal {
      /**
       * Here an edge is intended to hold both the operands and its associated
       * partial derivatives. They're held together in the
       * same class because then we can keep the templating logic that
       * specializes on type of operand in one place.
       *
       * This is the base template class that ends up getting instantiated
       * for arithmetic primitives (doubles and ints).
       *
       * @tparam ViewElt the type we expect to be inside each partials_[i]
       * @tparam Op the type of the operand(s)
       */
      template <typename ViewElt, typename Op>
      class ops_partials_edge {
      public:
        empty_broadcast_array<ViewElt> partials_;

        ops_partials_edge() {}
        explicit ops_partials_edge(const Op& /* op */) {}
        template<typename, typename, typename, typename, typename>
        friend class stan::math::operands_and_partials;
      private:
        void dump_partials(ViewElt* /* partials */) const {}
        void dump_operands(void* /* operands */) const {}
        ViewElt dx() const { return 0; }  // used for fvars

        /* size is used just for reverse mode to calculate the amount of memory
         * slots required to store the operands and partials. */
        int size() const { return 0; }
      };

      /**
       * Base class shared between fvar and var specializations of uni- and
       * multivariate edges will be in the prim/scal|arr|mat files.
       * Single param version - underlying Var is the same as Op, so there's
       * just one.
       *
       * @tparam ViewElt the type we expect to be inside each partials_[i]
       * @tparam Op the type of the operand(s)
       */
      template <typename ViewElt, typename Op>
      class ops_partials_edge_single {
      protected:
        template<typename, typename, typename, typename, typename>
        friend class stan::math::operands_and_partials;
        ViewElt partial_;
        const Op& operand_;
        int size() const { return 1; }
      public:
        broadcast_array<ViewElt> partials_;
        explicit ops_partials_edge_single(const Op& op)
          : partial_(0), operand_(op), partials_(partial_) {}
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
     * This base template is used when all operands are primitives
     * and we don't want to calculate derivatives at all. So all
     * Op1 - Op4 must be arithmetic primitives like int or double.
     *
     * @tparam Op1 type of the first operand
     * @tparam Op2 type of the second operand
     * @tparam Op3 type of the third operand
     * @tparam Op4 type of the fourth operand
     * @tparam T_return_type return type of the expression. This defaults
     *   to a template metaprogram that calculates the scalar promotion of
     *   Op1 -- Op4
     */
    template <typename Op1, typename Op2,
              typename Op3, typename Op4,
              typename T_return_type>
    class operands_and_partials {
    public:
      explicit operands_and_partials(const Op1& op1) {}
      operands_and_partials(const Op1& op1, const Op2& op2) {}
      operands_and_partials(const Op1& op1, const Op2& op2, const Op3& op3) {}
      operands_and_partials(const Op1& op1, const Op2& op2, const Op3& op3,
                            const Op4& op4) {}
      typedef double ViewElt;

      /**
       * Build the node to be stored on the autodiff graph.
       * This should contain both the value and the tangent.
       *
       * For scalars (this implementation), we don't calculate any tangents.
       * For reverse mode, we end up returning a type of var that will calculate
       * the appropriate adjoint using the stored operands and partials.
       * Forward mode just calculates the tangent on the spot and returns it in
       * a vanilla fvar.
       *
       * @param value the return value of the function we are compressing
       * @return the value with its derivative
       */
      T_return_type build(ViewElt value) {
        return value;
      }

      internal::ops_partials_edge<ViewElt, Op1> edge1_;
      internal::ops_partials_edge<ViewElt, Op2> edge2_;
      internal::ops_partials_edge<ViewElt, Op3> edge3_;
      internal::ops_partials_edge<ViewElt, Op4> edge4_;
    };
  }  // namespace math
}  // namespace stan
#endif
