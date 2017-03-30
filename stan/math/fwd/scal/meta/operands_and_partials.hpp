#ifndef STAN_MATH_FWD_SCAL_META_OPERANDS_AND_PARTIALS_HPP
#define STAN_MATH_FWD_SCAL_META_OPERANDS_AND_PARTIALS_HPP

#include <stan/math/prim/scal/meta/operands_and_partials.hpp>

namespace stan {
  namespace math {
    namespace detail {
      template <typename Op1, typename Op2, typename Op3, typename Op4>
      class operands_and_partials<Op1, Op2, Op3, Op4, var> {
      public:
        operands_and_partials(const Op1& o1)
          : edge1_(o1) { }
        operands_and_partials(const Op1& o1, const Op2& o2)
          : edge1_(o1), edge2_(o2) { }
        operands_and_partials(const Op1& o1, const Op2& o2, const Op3& o3)
          : edge1_(o1), edge2_(o2), edge3_(o3) { }
        operands_and_partials(const Op1& o1, const Op2& o2, const Op3& o3, const Op4& o4)
          : edge1_(o1), edge2_(o2), edge3_(o3), edge4_(o4) { }

        void increment_dx1(const int n, const double adj) { edge1_.increment_dx(n, adj); }
        void increment_dx2(const int n, const double adj) { edge2_.increment_dx(n, adj); }
        void increment_dx3(const int n, const double adj) { edge3_.increment_dx(n, adj); }
        void increment_dx4(const int n, const double adj) { edge4_.increment_dx(n, adj); }

        template<typename T>
        typename boost::enable_if_c<is_vector_like<Op1>::value
                                    && is_vector_like<T>::value, void>::type
        increment_dx1_vector(const int n, const T& adj) {
          edge1_.increment_dx_vector(n, adj);
        }
        template<typename T>
        typename boost::enable_if_c<is_vector_like<Op2>::value
                                    && is_vector_like<T>::value, void>::type
        increment_dx2_vector(const int n, const T& adj) {
          edge2_.increment_dx_vector(n, adj);
        }
        template<typename T>
        typename boost::enable_if_c<is_vector_like<Op3>::value
                                    && is_vector_like<T>::value, void>::type
        increment_dx3_vector(const int n, const T& adj) {
          edge3_.increment_dx_vector(n, adj);
        }
        template<typename T>
        typename boost::enable_if_c<is_vector_like<Op4>::value
                                    && is_vector_like<T>::value, void>::type
        increment_dx4_vector(const int n, const T& adj) {
          edge4_.increment_dx_vector(n, adj);
        }

        // this is what matters in terms of going on autodiff stack
        fvar build(double value) {
          size_t size = edge1_.size() + edge2_.size() + edge3_.size() + edge4_.size();
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

          return var(new precomputed_gradients_vari(value, size, varis, partials));
        };
      private:
        // these are going to be stack local and get collected
        // TODO(sean): should we pass in ViewElt somehow? shape of it?
        ops_partials_edge<double, Op1> edge1_;
        ops_partials_edge<double, Op2> edge2_;
        ops_partials_edge<double, Op3> edge3_;
        ops_partials_edge<double, Op4> edge4_;
      };
    } // end namespace detail
  }
}

#endif
