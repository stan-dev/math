#ifndef STAN_MATH_PRIM_SCAL_META_OPERANDS_AND_PARTIALS_HPP
#define STAN_MATH_PRIM_SCAL_META_OPERANDS_AND_PARTIALS_HPP

#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/prim/scal/meta/broadcast_array.hpp>

namespace stan {
  namespace math {
    namespace detail {
      template <typename ViewElt, typename Op>
      struct ops_partials_edge {
        empty_broadcast_array<ViewElt> partials_;

        ops_partials_edge() {}
        explicit ops_partials_edge(const Op& /* a */) {}
        //
        void dump_partials(ViewElt* /* partials */) const {}  // used for vars
        void dump_operands(void* /* operands */) const {}  // used for vars
        ViewElt dx() const { return 0; }  // used for fvars
        int size() const { return 0; }
      };

      // Base class shared between fvar and var specializations of uni- and
      // multivariate edges will be in the prim/scal|arr|mat files.
      // Single param version - underlying Var is the same as Op, so there's
      // just one.
      template <typename ViewElt, typename Op>
      class ops_partials_edge_single {
      private:
        ViewElt partial_;
      public:
        const Op& operand_;
        broadcast_array<ViewElt> partials_;
        explicit ops_partials_edge_single(const Op& a)
          : partial_(0), operand_(a), partials_(partial_) {}
        // dump_operands implemented in specialization
        int size() const { return 1; }
      };
    }  // end namespace detail

    /**
     * This class builds partial derivatives with respect to a set of
     * operands. There are two reason for the generality of this
     * class. The first is to handle vector and scalar arguments
     * without needing to write additional code. The second is to use
     * this class for writing probability distributions that handle
     * primitives, reverse mode, and forward mode variables
     * seamlessly.
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
    template <typename Op1 = double, typename Op2 = double,
              typename Op3 = double, typename Op4 = double,
              typename T_return_type
              = typename return_type<Op1, Op2, Op3, Op4>::type>
    class operands_and_partials {
    public:
      explicit operands_and_partials(const Op1& op1) {}
      operands_and_partials(const Op1& op1, const Op2& op2) {}
      operands_and_partials(const Op1& op1, const Op2& op2, const Op3& op3) {}
      operands_and_partials(const Op1& op1, const Op2& op2, const Op3& op3,
                            const Op4& op4) {}
      typedef double ViewElt;

      T_return_type build(ViewElt value) {
        return value;
      }

      detail::ops_partials_edge<ViewElt, Op1> edge1_;
      detail::ops_partials_edge<ViewElt, Op2> edge2_;
      detail::ops_partials_edge<ViewElt, Op3> edge3_;
      detail::ops_partials_edge<ViewElt, Op4> edge4_;
    };
  }  // end namespace math
}  // end namespace stan
#endif
