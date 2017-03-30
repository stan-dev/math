#ifndef STAN_MATH_PRIM_SCAL_META_OPERANDS_AND_PARTIALS_HPP
#define STAN_MATH_PRIM_SCAL_META_OPERANDS_AND_PARTIALS_HPP

#include <stan/math/prim/scal/meta/return_type.hpp>

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

      template <typename ViewElt, typename Op, typename Var>
      class ops_partials_edge {
      public:
        ops_partials_edge() { }
        ops_partials_edge(const Op& /* a */){ }
        void increment_dx(int /* n */, const ViewElt& /* adj */) {}
        void dump_partials(ViewElt* /* partials */) {}
        void dump_operands(Var** /* partials */) {}
        int size() {return 0; }
      };

      template <typename Op1 = double, typename Op2 = double,
                typename Op3 = double, typename Op4 = double,
                typename T_return_type = typename return_type<Op1, Op2, Op3, Op4>::type>
      class operands_and_partials {
      public:
        operands_and_partials(const Op1& op1) {}
        operands_and_partials(const Op1& op1, const Op2& op2) {}
        operands_and_partials(const Op1& op1, const Op2& op2, const Op3& op3) {}
        operands_and_partials(const Op1& op1, const Op2& op2, const Op3& op3, const Op4& op4) {}

        void increment_dx1(int, double) {}
        void increment_dx2(int, double) {}
        void increment_dx3(int, double) {}
        void increment_dx4(int, double) {}

        template <typename PartialVec>
        void increment_dx1_vector(int, PartialVec) {}
        template <typename PartialVec>
        void increment_dx2_vector(int, PartialVec) {}
        template <typename PartialVec>
        void increment_dx3_vector(int, PartialVec) {}
        template <typename PartialVec>
        void increment_dx4_vector(int, PartialVec) {}

        double build(const double value) {
          return value;
        }
      };

      // Base class shared between fvar and var specializations of uni- and
      // multivariate edges will be in the prim/scal|arr|mat files.
      // Singular version - underlying Var is the same as Op, so there's just one.
      template <typename ViewElt, typename Op, typename Var>
      class ops_partials_edge_singular<ViewElt, Op, Var> {
      protected:
        const Op& operand;
        ViewElt partial;

      public:
        ops_partials_edge(const Op& a) : operand(a), partial() { }
        void increment_dx(int /*n*/, const ViewElt& adj) {
          partial += adj;
        }
        void dump_partials(ViewElt* partials) {
          *partials = partial;
        }
        // dump_operands implemented in specialization
        int size() {
          return 1;
        }
      };

      // One of these for both vars and fvars
      template <typename Op1, typename Op2, typename Op3, typename Op4>
      class operands_and_partials_base<Op1, Op2, Op3, Op4, var> {
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
        var build(double value) {
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
          edge4_.dump_partkals(&partials[idx += edge3_.size()]);

          return var(new precomputed_gradients_vari(value, size, varis, partials));
        };
      private:
        // these are going to be stack local and get collected
        // TODO(sean): should we pass in ViewElt somehow? shape of it?
        ops_partials_edge<double, Op1, vari> edge1_;
        ops_partials_edge<double, Op2, vari> edge2_;
        ops_partials_edge<double, Op3, vari> edge3_;
        ops_partials_edge<double, Op4, vari> edge4_;
      };
    } // end namespace detail
  } // end namespace math
} // end namespace stan
#endif
