#ifndef STAN_MATH_REV_SCAL_META_OPERANDS_AND_PARTIALS_HPP
#define STAN_MATH_REV_SCAL_META_OPERANDS_AND_PARTIALS_HPP

#include <stan/math/rev/core/chainablestack.hpp>
#include <stan/math/rev/core/precomputed_gradients.hpp>
#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/core/vari.hpp>
#include <stan/math/prim/scal/meta/operands_and_partials.hpp>
#include <stan/math/prim/scal/meta/is_vector_like.hpp>

namespace stan {
  namespace math {
    namespace detail {
      // ViewElt = double, Arg = var d_ holds double
      template <typename ViewElt>
      class ops_partials_edge<ViewElt, var>
        : public ops_partials_edge_singular<ViewElt, var> {
      public:
        explicit ops_partials_edge(const var& op)
          : ops_partials_edge_singular<ViewElt, var>(op) {}
        void dump_partials(ViewElt* partials) {
          *partials = this->partials[0];
        }
        void dump_operands(vari** varis) {
          *varis = this->operand.vi_;
        }
      };
    }  // end namespace detail
    template <typename Op1, typename Op2, typename Op3, typename Op4>
    class operands_and_partials<Op1, Op2, Op3, Op4, var> {
    public:
      // these are going to be stack local and get collected
      detail::ops_partials_edge<double, Op1> edge1_;
      detail::ops_partials_edge<double, Op2> edge2_;
      detail::ops_partials_edge<double, Op3> edge3_;
      detail::ops_partials_edge<double, Op4> edge4_;
      explicit operands_and_partials(const Op1& o1)
        : edge1_(o1) { }
      operands_and_partials(const Op1& o1, const Op2& o2)
        : edge1_(o1), edge2_(o2) { }
      operands_and_partials(const Op1& o1, const Op2& o2, const Op3& o3)
        : edge1_(o1), edge2_(o2), edge3_(o3) { }
      operands_and_partials(const Op1& o1, const Op2& o2, const Op3& o3,
                            const Op4& o4)
        : edge1_(o1), edge2_(o2), edge3_(o3), edge4_(o4) { }

      // this is what matters in terms of going on autodiff stack
      var build(double value) {
        size_t size = edge1_.size() + edge2_.size() + edge3_.size()
          + edge4_.size();
        vari** varis = ChainableStack::memalloc_.alloc_array<vari*>(size);
        int idx = 0;
        edge1_.dump_operands(&varis[idx]);
        edge2_.dump_operands(&varis[idx += edge1_.size()]);
        edge3_.dump_operands(&varis[idx += edge2_.size()]);
        edge4_.dump_operands(&varis[idx += edge3_.size()]);

        double* partials = ChainableStack::memalloc_.alloc_array<double>(size);
        idx = 0;
        edge1_.dump_partials(&partials[idx]);
        edge2_.dump_partials(&partials[idx += edge1_.size()]);
        edge3_.dump_partials(&partials[idx += edge2_.size()]);
        edge4_.dump_partials(&partials[idx += edge3_.size()]);

        return var(new
                   precomputed_gradients_vari(value, size, varis, partials));
      }
    };
  }
}
#endif
