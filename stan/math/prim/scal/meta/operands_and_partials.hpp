#ifndef STAN_MATH_PRIM_SCAL_META_OPERANDS_AND_PARTIALS_HPP
#define STAN_MATH_PRIM_SCAL_META_OPERANDS_AND_PARTIALS_HPP

#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/prim/scal/meta/broadcast_array.hpp>

namespace stan {
  namespace math {
    namespace detail {
      // VIEWS of SCALARS
      // ====================================
      // ViewElt = double, Arg = var d_ holds double
      // ViewElt = double, Arg = std::vector<var> d_ holds vector<double>
      // ViewElt = double, Arg = vector_v d_ holds vector<double>
      // ViewElt = double, Arg = row_vector_v d_ holds = vector<double>

      // dummies
      // ViewElt = double, Arg = double, d_ holds dummy
      // ViewElt = double, Arg = std::vector<double> d_ holds dummy
      // ViewElt = double, Arg = VectorXd d_ holds dummy
      // ViewElt = double, Arg = RowVectorXd, d_ holds dummy

      // VIEWS of VECTORS (R = -1, C = 1)
      // =====================
      // ViewElt = VectorXd, Arg = vector_v, d_ holds VectorXd
      // ViewElt = VectorXd, Arg = vector<vector_v> d_ holds vector<VectorXd>
      // vector<VectrorXd>

      // dummies
      // ViewElt = VectorXd, Arg = VectorXd d_ holds dummy
      // ViewElt = VectorXd, Arg = vector<VectorXd> d_ holds dummy

      // VIEWS of ROW VECTORS (R = 1, C = -1)
      // =====================
      // ViewElt = RowVectorXd, Arg = row_vector_v, d_ holds RowVectorXd
      // ViewElt = RowVectorXd, Arg = vector<row_vector_v> d_ holds vector<RowVectorXd>

      // dummies
      // ViewElt = RowVectorXd, Arg = RowVectorXd d_ holds dummy
      // ViewElt = RowVectorXd, Arg = vector<RowVectorXd> d_ holds dummy

      // VIEWS of MATRICES (R = -1, C = -1)
      // =====================
      // ViewElt = MatrixXd, Arg = matrix_v, d_ holds MatrixXd
      // ViewElt = MatrixXd, Arg = vector<matrix_v> d_ holds vector<MatriXd>

      // dummies
      // ViewElt = MatrixXd, Arg = MatrixXd d_ holds dummy
      // ViewElt = MatrixXd, Arg = vector<MatrixXd> d_ holds dummy

      template <typename ViewElt, typename Op>
      class ops_partials_edge {
      public:
        ops_partials_edge() { }
        ops_partials_edge(const Op& /* a */){ }
        void dump_partials(ViewElt* /* partials */) {} // used for vars
        void dump_operands(void* /* operands */) {} // also used for reverse mode
        double dx() const { return 0; } //used for fvars
        int size() { return 0; }
      };

      // Base class shared between fvar and var specializations of uni- and
      // multivariate edges will be in the prim/scal|arr|mat files.
      // Singular version - underlying Var is the same as Op, so there's just one.
      template <typename ViewElt, typename Op>
      class ops_partials_edge_singular {
      public:
        const Op& operand;
        ViewElt partial;
        broadcast_array<ViewElt> partials;
        ops_partials_edge_singular(const Op& a) : operand(a), partial(0), partials(partial) { }
        // dump_operands implemented in specialization
        int size() { return 1; }
      };

      template<>
      struct ops_partials_edge<double, double> {
        ops_partials_edge() {}
        ops_partials_edge(const double& /* a */) {}
        broadcast_array<void> partials;
        void dump_operands(void*) {}
        void dump_partials(void*) {}
        double dx() { return 0; }
        int size() { return 0; }
      };
    } // end namespace detail
    template <typename Op1 = double, typename Op2 = double,
              typename Op3 = double, typename Op4 = double,
              typename T_return_type = typename return_type<Op1, Op2, Op3, Op4>::type>
    class operands_and_partials {
    public:
      operands_and_partials(const Op1& op1) {}
      operands_and_partials(const Op1& op1, const Op2& op2) {}
      operands_and_partials(const Op1& op1, const Op2& op2, const Op3& op3) {}
      operands_and_partials(const Op1& op1, const Op2& op2, const Op3& op3, const Op4& op4) {}

      double build(const double value) {
        return value;
      }

      detail::ops_partials_edge<double, Op1> edge1_;
      detail::ops_partials_edge<double, Op2> edge2_;
      detail::ops_partials_edge<double, Op3> edge3_;
      detail::ops_partials_edge<double, Op4> edge4_;
    };
  } // end namespace math
} // end namespace stan
#endif
