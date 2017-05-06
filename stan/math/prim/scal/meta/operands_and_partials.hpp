#ifndef STAN_MATH_PRIM_SCAL_META_OPERANDS_AND_PARTIALS_HPP
#define STAN_MATH_PRIM_SCAL_META_OPERANDS_AND_PARTIALS_HPP

#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/prim/scal/meta/broadcast_array.hpp>

namespace stan {
  namespace math {
    namespace detail {
      template <typename ViewElt, typename Op>
      class ops_partials_edge {
      public:
        ops_partials_edge() {}
        explicit ops_partials_edge(const Op& /* a */) {}
        empty_broadcast_array<ViewElt> partials;
        void dump_partials(ViewElt* /* partials */) {}  // used for vars
        void dump_operands(void* /* operands */) {}  // also used for reverse
        double dx() const { return 0; }  // used for fvars
        int size() { return 0; }
      };

      // Base class shared between fvar and var specializations of uni- and
      // multivariate edges will be in the prim/scal|arr|mat files.
      // Singular version - underlying Var is the same as Op, so there's
      // just one.
      template <typename ViewElt, typename Op>
      class ops_partials_edge_singular {
      public:
        const Op& operand;
        ViewElt partial;
        broadcast_array<ViewElt> partials;
        explicit ops_partials_edge_singular(const Op& a)
          : operand(a), partial(0), partials(partial) {}
        // dump_operands implemented in specialization
        int size() { return 1; }
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
     * @tparam Op1 First set of operands.
     * @tparam Op2 Second set of operands.
     * @tparam Op3 Third set of operands.
     * @tparam Op4 Fourth set of operands.
     * @tparam T_return_type Return type of the expression. This defaults
     *   to a template metaprogram that calculates the scalar promotion of
     *   Op1 -- Op4.
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

      double build(const double value) {
        return value;
      }

      detail::ops_partials_edge<double, Op1> edge1_;
      detail::ops_partials_edge<double, Op2> edge2_;
      detail::ops_partials_edge<double, Op3> edge3_;
      detail::ops_partials_edge<double, Op4> edge4_;
    };
  }  // end namespace math
}  // end namespace stan
#endif
