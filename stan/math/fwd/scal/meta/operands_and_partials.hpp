#ifndef STAN_MATH_FWD_SCAL_META_OPERANDS_AND_PARTIALS_HPP
#define STAN_MATH_FWD_SCAL_META_OPERANDS_AND_PARTIALS_HPP

#include <stan/math/prim/scal/meta/operands_and_partials.hpp>

namespace stan {
  namespace math {
    namespace detail {
      // For fvars, each of these must implement a dx() method that calculates the
      // contribution to the derivative of the o&p node for the edge.
      template <typename ViewElt, typename Dx>
      class ops_partials_edge<ViewElt, fvar<Dx> >
        : public ops_partials_edge_singular<ViewElt, fvar<Dx> > {
      public:
        ops_partials_edge(const fvar<Dx>& op)
          : ops_partials_edge_singular<ViewElt, fvar<Dx> >(op) {}
        Dx dx() {
          return this->partials[0] * this->operand.d_;
        }
      };
    } // end namespace detail
    template <typename Op1, typename Op2, typename Op3, typename Op4,
              typename Dx>
    class operands_and_partials<Op1, Op2, Op3, Op4, fvar<Dx> > {
    public:
      // these are going to be stack local and get collected
      detail::ops_partials_edge<Dx, Op1> edge1_;
      detail::ops_partials_edge<Dx, Op2> edge2_;
      detail::ops_partials_edge<Dx, Op3> edge3_;
      detail::ops_partials_edge<Dx, Op4> edge4_;
      typedef fvar<Dx> T_return_type;
      operands_and_partials(const Op1& o1)
        : edge1_(o1) { }
      operands_and_partials(const Op1& o1, const Op2& o2)
        : edge1_(o1), edge2_(o2) { }
      operands_and_partials(const Op1& o1, const Op2& o2, const Op3& o3)
        : edge1_(o1), edge2_(o2), edge3_(o3) { }
      operands_and_partials(const Op1& o1, const Op2& o2, const Op3& o3, const Op4& o4)
        : edge1_(o1), edge2_(o2), edge3_(o3), edge4_(o4) { }

      // this is what matters in terms of going on autodiff stack
      T_return_type build(Dx value) {
        // Just call chain here and end up with an fvar(val, adj)
        Dx deriv = 0;
        deriv += edge1_.dx();
        deriv += edge2_.dx();
        deriv += edge3_.dx();
        deriv += edge4_.dx();
        return T_return_type(value, deriv);
      };
    };
  }
}

#endif
