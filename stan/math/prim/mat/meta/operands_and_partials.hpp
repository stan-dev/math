#ifndef STAN_MATH_PRIM_MAT_META_OPERANDS_AND_PARTIALS_HPP
#define STAN_MATH_PRIM_MAT_META_OPERANDS_AND_PARTIALS_HPP

#include <stan/math/rev/mat/fun/typedefs.hpp>
#include <stan/math/rev/scal/meta/operands_and_partials.hpp>
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
      public:
        typedef Eigen::Matrix<ViewElt, R, C> partials_t;
        ops_partials_edge_mat_prim(const operands_t& ops)
          : operands(ops), partials(partials_t::Zero(ops.rows(), ops.cols())) {

        }
        void increment_dx_vector(int /*n*/, const partials_t& adj) {
          partials += adj;
        }
        void dump_partials(double* partials) {
          for (int i = 0; i < this->partials.size(); ++i) {
            partials[i] = this->partials(i);
          }
        }
        void dump_operands(vari** varis) {
          for (int i = 0; i < this->operands.size(); ++i) {
            varis[i] = this->operands(i).vi_;
          }
        }
        int size() {
          return this->operands.size();
        }
      protected:
        const Op& operands;
        partials_t partials;
      };

      template <typename ViewElt, typename Op>
      class ops_partials_edge_multivariate_prim {
      public:
        typedef Eigen::Matrix<ViewElt, Eigen::Dynamic, Eigen::Dynamic> partial_t;
        explicit ops_partials_edge_multivariate_prim(const std::vector<Op>& ops)
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
