#ifndef STAN_MATH_PRIM_MAT_META_OPERANDS_AND_PARTIALS_HPP
#define STAN_MATH_PRIM_MAT_META_OPERANDS_AND_PARTIALS_HPP

#include <stan/math/prim/mat/fun/typedefs.hpp>
#include <stan/math/prim/scal/meta/operands_and_partials.hpp>
#include <vector>

namespace stan {
  namespace math {
    namespace detail {
      // Handles all the vectorized cases for "views of scalars":
      // ViewElt = double, Arg = std::vector<var> d_ holds vector<double>
      // ViewElt = double, Arg = std::vector<double> d_ holds dummy
      // Op will always be some container of vars
      template <typename ViewElt, typename Op, int R, int C>
      class ops_partials_edge_mat_prim {
      protected:
        typedef Eigen::Matrix<ViewElt, R, C> partials_t;
        const Op& operands;
        partials_t partials;
      public:
        ops_partials_edge_mat_prim(const Op& ops)
          : operands(ops), partials(partials_t::Zero(ops.size())) {}
        void increment_dx(int n, const ViewElt& adj) {
          partials[n] += adj;
        }
        void increment_dx_vector(int /*n*/, const partials_t& adj) {
          partials += adj;
        }
        void dump_partials(double* partials) {
          for (int i = 0; i < this->partials.size(); ++i) {
            partials[i] = this->partials(i);
          }
        }
        int size() {
          return this->operands.size();
        }
      };

      template <typename ViewElt, typename Op>
      class ops_partials_edge_multivariate_prim {
      public:
        typedef Eigen::Matrix<ViewElt, Eigen::Dynamic, Eigen::Dynamic> partial_t;
        ops_partials_edge_multivariate_prim(const std::vector<Op>& ops)
          : partials(ops.size()), operands(ops) {
          for (size_t i = 0; i < ops.size(); ++i) {
            partials[i] = partial_t::Zero(ops[i].rows(), ops[i].cols());
          }
        }
        void increment_dx_vector(int n, const partial_t& adj) {
          partials[n] += adj;
        }
        int size() {
          return this->operands.size() * this->operands[0].size();
        }
      protected:
        std::vector<partial_t> partials;
        const std::vector<Op>& operands;
      };
    }
  }
}
#endif
