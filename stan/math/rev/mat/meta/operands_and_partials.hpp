#ifndef STAN_MATH_REV_MAT_META_OPERANDS_AND_PARTIALS_HPP
#define STAN_MATH_REV_MAT_META_OPERANDS_AND_PARTIALS_HPP

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
      template <typename ViewElt, typename Op>
      class ops_partials_edge_vec {
      protected:
        const Op& operands;
        typedef Eigen::Matrix<ViewElt, 1, Eigen::Dynamic> partials_t;
        partials_t partials;

      public:
        explicit ops_partials_edge_vec(const Op& ops)
          : operands(ops), partials(partials_t::Zero(ops.size())) {}

        void increment_dx(int n, const ViewElt& adj) {
          partials[n] += adj;
        }
        void increment_dx_vector(int /* n */, const partials_t& adj) {
          partials += adj;
        }
        void dump_partials(double* partials) {
          for (int i = 0; i < size(); ++i) {
            partials[i] = this->partials[i];
          }
        }
        void dump_operands(vari** varis) {
          for (int i = 0; i < size(); ++i) {
            varis[i] = operands[i].vi_;
          }
        }
        int size() {
          return partials.size();
        }
      };

      // ViewElt = double, Arg = std::vector<var> d_ holds vector<double>
      template <typename ViewElt>
      class ops_partials_edge<ViewElt, std::vector<var> >
        : public ops_partials_edge_vec<ViewElt, std::vector<var> > {
      public:
        explicit ops_partials_edge(const std::vector<var>& ops)
          : ops_partials_edge_vec<ViewElt, std::vector<var> >(ops) {}
      };

      // ViewElt = double, Arg = vector_v d_ holds vector<double>
      // ViewElt = double, Arg = VectorXd d_ holds dummy
      template <typename ViewElt>
      class ops_partials_edge<ViewElt, vector_v >
        : public ops_partials_edge_vec<ViewElt, vector_v > {
      public:
        explicit ops_partials_edge(const vector_v& ops)
          : ops_partials_edge_vec<ViewElt, vector_v>(ops) {}
      };

      // ViewElt = double, Arg = row_vector_v d_ holds = vector<double>
      // ViewElt = double, Arg = RowVectorXd, d_ holds dummy
      template <typename ViewElt>
      class ops_partials_edge<ViewElt, row_vector_v >
        : public ops_partials_edge_vec<ViewElt, row_vector_v > {
      public:
        explicit ops_partials_edge(const row_vector_v& ops)
          : ops_partials_edge_vec<ViewElt, row_vector_v>(ops) {}
      };

      // TODO(Sean): multivariate cases below
      // MULTIVARIATE
      // VIEWS of ROW VECTORS (R = 1, C = -1)
      // =====================
      // ViewElt = RowVectorXd, Arg = row_vector_v, d_ holds RowVectorXd
      template <typename ViewElt, typename Op>
      class ops_partials_edge_multivariate {
      private:
        typedef Eigen::Matrix<ViewElt, Eigen::Dynamic, Eigen::Dynamic> partial_t;
        std::vector<partial_t> partials;
      public:
        const std::vector<Op>& operands;
        explicit ops_partials_edge_multivariate(const std::vector<Op>& ops)
          : partials(ops.size()), operands(ops) {
          for (size_t i = 0; i < ops.size(); ++i) {
            partials[i] = partial_t::Zero(ops[i].rows(), ops[i].cols());
          }
        }
        void increment_dx_vector(int n, const partial_t& adj) {
          partials[n] += adj;
        }
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
        int size() {
          return this->operands.size() * this->operands[0].size();
        }
      };

      // ViewElt = RowVectorXd, Arg = vector<row_vector_v> d_ holds vector<RowVectorXd>
      template <typename ViewElt>
      class ops_partials_edge<ViewElt, std::vector<row_vector_v> >
        : public ops_partials_edge_multivariate<ViewElt, row_vector_v> {
      public:
        explicit ops_partials_edge(const std::vector<row_vector_v>& ops)
          : ops_partials_edge_multivariate<ViewElt, row_vector_v>(ops) {}
      };

      template <typename ViewElt>
      class ops_partials_edge<ViewElt, std::vector<vector_v> >
        : public ops_partials_edge_multivariate<ViewElt, vector_v> {
      public:
        explicit ops_partials_edge(const std::vector<vector_v>& ops)
          : ops_partials_edge_multivariate<ViewElt, vector_v>(ops) {}
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
        : public ops_partials_edge_multivariate<ViewElt, matrix_v> {
      public:
        explicit ops_partials_edge(const std::vector<matrix_v>& ops)
          : ops_partials_edge_multivariate<ViewElt, matrix_v>(ops) {}
      };

    } // end namespace detail
  } // end namespace math
} // end namespace stan
#endif
