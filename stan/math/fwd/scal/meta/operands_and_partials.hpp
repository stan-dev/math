#ifndef STAN_MATH_FWD_SCAL_META_OPERANDS_AND_PARTIALS_HPP
#define STAN_MATH_FWD_SCAL_META_OPERANDS_AND_PARTIALS_HPP

#include <stan/math/prim/scal/meta/operands_and_partials.hpp>
#include <stan/math/fwd/core/fvar.hpp>

namespace stan {
  namespace math {
    namespace detail {
      // For fvars, each of these must implement a dx() method that calculates
      // the contribution to the derivative of the o&p node for the edge.
      template <typename ViewElt, typename Dx>
      class ops_partials_edge<ViewElt, fvar<Dx> >
        : public ops_partials_edge_single<ViewElt, fvar<Dx> > {
      public:
        explicit ops_partials_edge(const fvar<Dx>& op)
          : ops_partials_edge_single<ViewElt, fvar<Dx> >(op) {}
        Dx dx() {
          return this->partials_[0] * this->operand_.d_;
        }
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
     * This is the specialization for when the return type is fvar,
     * which should be for forward mode and all higher-order cases.
     *
     * @tparam Op1 type of the first operand
     * @tparam Op2 type of the second operand
     * @tparam Op3 type of the third operand
     * @tparam Op4 type of the fourth operand
     * @tparam T_return_type return type of the expression. This defaults
     *   to a template metaprogram that calculates the scalar promotion of
     *   Op1 -- Op4
     */
    template <typename Op1, typename Op2, typename Op3, typename Op4,
              typename Dx>
    class operands_and_partials<Op1, Op2, Op3, Op4, fvar<Dx> > {
    public:
      detail::ops_partials_edge<Dx, Op1> edge1_;
      detail::ops_partials_edge<Dx, Op2> edge2_;
      detail::ops_partials_edge<Dx, Op3> edge3_;
      detail::ops_partials_edge<Dx, Op4> edge4_;
      typedef fvar<Dx> T_return_type;
      explicit operands_and_partials(const Op1& o1)
        : edge1_(o1) { }
      operands_and_partials(const Op1& o1, const Op2& o2)
        : edge1_(o1), edge2_(o2) { }
      operands_and_partials(const Op1& o1, const Op2& o2, const Op3& o3)
        : edge1_(o1), edge2_(o2), edge3_(o3) { }
      operands_and_partials(const Op1& o1, const Op2& o2, const Op3& o3,
                            const Op4& o4)
        : edge1_(o1), edge2_(o2), edge3_(o3), edge4_(o4) { }

      T_return_type build(Dx value) {
        Dx deriv = edge1_.dx() + edge2_.dx() + edge3_.dx() + edge4_.dx();
        return T_return_type(value, deriv);
      }
    };
  }
}
#endif
