#ifndef STAN_MATH_REV_MAT_META_OPERANDS_AND_PARTIALS_HPP
#define STAN_MATH_REV_MAT_META_OPERANDS_AND_PARTIALS_HPP

#include <stan/math/rev/mat/fun/typedefs.hpp>
#include <stan/math/rev/scal/meta/operands_and_partials.hpp>
#include <stan/math/prim/mat/meta/operands_and_partials.hpp>
#include <vector>

namespace stan {
  namespace math {
    namespace detail {
      // Vectorized Univariate
      template <typename ViewElt>
      class ops_partials_edge<ViewElt, std::vector<var> >
        : public ops_partials_edge_mat<ViewElt, std::vector<var>, -1, 1> {
      public:
        explicit ops_partials_edge(const std::vector<var>& ops)
          : ops_partials_edge_mat<ViewElt, std::vector<var>, -1, 1>(ops) {}
        void dump_operands(vari** varis) {
          for (size_t i = 0; i < this->operands.size(); ++i) {
            varis[i] = this->operands[i].vi_;
          }
        }
      };
      template <typename ViewElt, int R, int C>
      class ops_partials_edge<ViewElt, Eigen::Matrix<var, R, C> >
        : public ops_partials_edge_mat<ViewElt,
                                       Eigen::Matrix<var, R, C>, R, C> {
      public:
        explicit ops_partials_edge(const Eigen::Matrix<var, R, C>& ops)
          : ops_partials_edge_mat<ViewElt, Eigen::Matrix<var, R, C>, R, C>(ops)
        {}
        void dump_operands(vari** varis) {
          for (int i = 0; i < this->operands.size(); ++i) {
            varis[i] = this->operands(i).vi_;
          }
        }
      };

      // MULTIVARIATE
      template <typename ViewElt, typename Op>
      class ops_partials_edge_multi_rev
        : public ops_partials_edge_multivariate_prim<ViewElt, Op> {
      public:
        explicit ops_partials_edge_multi_rev(const std::vector<Op>& ops)
          : ops_partials_edge_multivariate_prim<ViewElt, Op>(ops) {}
        void dump_partials(double* partials) {
          int p_i = 0;
          for (size_t r = 0; r < this->partials_vec.size(); ++r) {
            for (int c = 0; c < this->partials_vec[r].size(); ++c, ++p_i) {
              partials[p_i] = this->partials_vec[r](c);
            }
          }
        }
        void dump_operands(vari** varis) {
          int p_i = 0;
          for (size_t r = 0; r < this->operands.size(); ++r) {
            for (int c = 0; c < this->operands[r].size(); ++c, ++p_i) {
              varis[p_i] = this->operands[r](c).vi_;
            }
          }
        }
      };

      template <typename ViewElt>
      class ops_partials_edge<ViewElt, std::vector<row_vector_v> >
        : public ops_partials_edge_multi_rev<ViewElt, row_vector_v> {
      public:
        explicit ops_partials_edge(const std::vector<row_vector_v>& ops)
          : ops_partials_edge_multi_rev<ViewElt, row_vector_v>(ops) {}
      };

      template <typename ViewElt>
      class ops_partials_edge<ViewElt, std::vector<vector_v> >
        : public ops_partials_edge_multi_rev<ViewElt, vector_v> {
      public:
        explicit ops_partials_edge(const std::vector<vector_v>& ops)
          : ops_partials_edge_multi_rev<ViewElt, vector_v>(ops) {}
      };

      template <typename ViewElt>
      class ops_partials_edge<ViewElt, std::vector<matrix_v> >
        : public ops_partials_edge_multi_rev<ViewElt, matrix_v> {
      public:
        explicit ops_partials_edge(const std::vector<matrix_v>& ops)
          : ops_partials_edge_multi_rev<ViewElt, matrix_v>(ops) {}
      };

    }  // end namespace detail
  }  // end namespace math
}  // end namespace stan
#endif
