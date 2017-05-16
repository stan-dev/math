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
      private:
        template<typename T1, typename T2, typename T3, typename T4, typename T5>
        friend class stan::math::operands_and_partials;
        void dump_operands(vari** varis) {
          for (size_t i = 0; i < this->operands_.size(); ++i) {
            varis[i] = this->operands_[i].vi_;
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
      private:
        template<typename T1, typename T2, typename T3, typename T4, typename T5>
        friend class stan::math::operands_and_partials;
        void dump_operands(vari** varis) {
          for (int i = 0; i < this->operands_.size(); ++i) {
            varis[i] = this->operands_(i).vi_;
          }
        }
      };

      // MULTIVARIATE
      template <typename ViewElt, int R, int C>
      class ops_partials_edge<ViewElt, std::vector<Eigen::Matrix<var, R, C> > >
        : public ops_partials_edge_multivariate_prim
      <ViewElt, std::vector<Eigen::Matrix<var, R, C> > > {
      public:
        typedef std::vector<Eigen::Matrix<var, R, C> > Op;
        explicit ops_partials_edge(const Op& ops)
          : ops_partials_edge_multivariate_prim<ViewElt, Op>(ops) {}
      private:
        template<typename T1, typename T2, typename T3, typename T4, typename T5>
        friend class stan::math::operands_and_partials;
        void dump_partials(double* partials) {
          int p_i = 0;
          for (size_t i = 0; i < this->partials_vec_.size(); ++i) {
            for (int j = 0; j < this->partials_vec_[i].size(); ++j, ++p_i) {
              partials[p_i] = this->partials_vec_[i](j);
            }
          }
        }
        void dump_operands(vari** varis) {
          int p_i = 0;
          for (size_t i = 0; i < this->operands_.size(); ++i) {
            for (int j = 0; j < this->operands_[i].size(); ++j, ++p_i) {
              varis[p_i] = this->operands_[i](j).vi_;
            }
          }
        }
      };
    }  // end namespace detail
  }  // end namespace math
}  // end namespace stan
#endif
