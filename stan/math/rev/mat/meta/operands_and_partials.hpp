#ifndef STAN_MATH_REV_MAT_META_OPERANDS_AND_PARTIALS_HPP
#define STAN_MATH_REV_MAT_META_OPERANDS_AND_PARTIALS_HPP

#include <stan/math/rev/mat/fun/typedefs.hpp>
#include <stan/math/rev/scal/meta/operands_and_partials.hpp>
#include <stan/math/prim/mat/meta/operands_and_partials.hpp> // XXX are both of these includes okay
#include <vector>

namespace stan {
  namespace math {
    namespace detail {
      template <typename ViewElt, typename Op>
      class ops_partials_edge_vec_rev
        : public ops_partials_edge_mat_prim<ViewElt, Op, -1, 1> {
      public:
        ops_partials_edge_vec_rev(const Op& ops)
          : ops_partials_edge_mat_prim<ViewElt, Op, -1, 1>(ops) {}
        void dump_partials(double* partials) {
          for (int i = 0; i < this->size(); ++i) {
            partials[i] = this->partials[i];
          }
        }
        void dump_operands(vari** varis) {
          for (int i = 0; i < this->size(); ++i) {
            varis[i] = this->operands[i].vi_;
          }
        }
      };

      // ViewElt = double, Arg = std::vector<var> d_ holds vector<double>
      template <typename ViewElt>
      class ops_partials_edge<ViewElt, std::vector<var> >
        : public ops_partials_edge_vec_rev<ViewElt, std::vector<var> > {
      public:
        explicit ops_partials_edge(const std::vector<var>& ops)
          : ops_partials_edge_vec_rev<ViewElt, std::vector<var> >(ops) {}
      };

      // ViewElt = double, Arg = vector_v d_ holds vector<double>
      // ViewElt = double, Arg = VectorXd d_ holds dummy
      template <typename ViewElt>
      class ops_partials_edge<ViewElt, vector_v >
        : public ops_partials_edge_vec_rev<ViewElt, vector_v > {
      public:
        explicit ops_partials_edge(const vector_v& ops)
          : ops_partials_edge_vec_rev<ViewElt, vector_v>(ops) {}
      };

      // ViewElt = double, Arg = row_vector_v d_ holds = vector<double>
      // ViewElt = double, Arg = RowVectorXd, d_ holds dummy
      template <typename ViewElt>
      class ops_partials_edge<ViewElt, row_vector_v >
        : public ops_partials_edge_vec_rev<ViewElt, row_vector_v > {
      public:
        explicit ops_partials_edge(const row_vector_v& ops)
          : ops_partials_edge_vec_rev<ViewElt, row_vector_v>(ops) {}
      };

      // MULTIVARIATE
      // VIEWS of ROW VECTORS (R = 1, C = -1)
      // =====================
      // ViewElt = RowVectorXd, Arg = row_vector_v, d_ holds RowVectorXd
      // ViewElt = RowVectorXd, Arg = vector<row_vector_v> d_ holds vector<RowVectorXd>

      template <typename ViewElt, typename Op>
      class ops_partials_edge_multi_rev
        : public ops_partials_edge_multivariate_prim<ViewElt, Op> {
      public:
        ops_partials_edge_multi_rev(const std::vector<Op>& ops)
          : ops_partials_edge_multivariate_prim<ViewElt, Op>(ops) {}
        void dump_partials(double* partials) {
          int p_i = 0;
          for (size_t r = 0; r < this->partials.size(); ++r) {
            for (int c = 0; c < this->partials[r].size(); ++c, ++p_i) {
              partials[p_i] = this->partials[r](c);
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

      // VIEWS of MATRICES (R = -1, C = -1)
      // =====================
      // ViewElt = MatrixXd, Arg = matrix_v, d_ holds MatrixXd
      // ViewElt = MatrixXd, Arg = MatrixXd d_ holds dummy
      template <typename ViewElt, int R, int C>
      class ops_partials_edge<ViewElt, Eigen::Matrix<var, R, C> > {
      public:
        typedef Eigen::Matrix<ViewElt, R, C> partials_t;
        typedef Eigen::Matrix<var, R, C> operands_t;
        ops_partials_edge(const operands_t& ops)
          : operands(ops), partials(partials_t::Zero(ops.rows(), ops.cols())) {}
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
        const operands_t& operands;
        partials_t partials;
      };

      // ViewElt = MatrixXd, Arg = vector<matrix_v> d_ holds vector<MatriXd>
      // ViewElt = MatrixXd, Arg = vector<MatrixXd> d_ holds dummy
      template <typename ViewElt>
      class ops_partials_edge<ViewElt, std::vector<matrix_v> >
        : public ops_partials_edge_multi_rev<ViewElt, matrix_v> {
      public:
        ops_partials_edge(const std::vector<matrix_v>& ops)
          : ops_partials_edge_multi_rev<ViewElt, matrix_v>(ops) {}
      };

    } // end namespace detail
  } // end namespace math
} // end namespace stan
#endif
