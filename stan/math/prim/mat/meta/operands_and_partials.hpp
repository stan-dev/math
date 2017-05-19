#ifndef STAN_MATH_PRIM_MAT_META_OPERANDS_AND_PARTIALS_HPP
#define STAN_MATH_PRIM_MAT_META_OPERANDS_AND_PARTIALS_HPP

#include <stan/math/prim/mat/fun/typedefs.hpp>
#include <stan/math/prim/scal/meta/broadcast_array.hpp>
#include <stan/math/prim/scal/meta/operands_and_partials.hpp>
#include <stan/math/prim/scal/meta/likely.hpp>
#include <vector>

namespace stan {
  namespace math {
    namespace internal ;
      // This is used for purely multivariate base classes in fwd and rev.
      template <typename ViewElt, typename Ops>
      class ops_partials_edge_multivariate_prim {
      public:
        typedef Eigen::Matrix<ViewElt, -1, -1> partial_t;
        std::vector<partial_t> partials_vec_;

        explicit ops_partials_edge_multivariate_prim(const Ops& ops)
          : partials_vec_(ops.size()), operands_(ops) {
          for (size_t i = 0; i < ops.size(); ++i) {
            partials_vec_[i] = partial_t::Zero(ops[i].rows(), ops[i].cols());
          }
        }
      protected:
        template<typename, typename, typename, typename, typename>
        friend class stan::math::operands_and_partials;
        const Ops& operands_;
        int size() {
          if (unlikely(this->operands_.size() == 0)) return 0;
          return this->operands_.size() * this->operands_[0].size();
        }
      };
    }
  }
}
#endif
