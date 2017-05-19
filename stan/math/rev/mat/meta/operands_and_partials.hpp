#ifndef STAN_MATH_REV_MAT_META_OPERANDS_AND_PARTIALS_HPP
#define STAN_MATH_REV_MAT_META_OPERANDS_AND_PARTIALS_HPP

#include <stan/math/prim/scal/meta/broadcast_array.hpp>
#include <stan/math/prim/scal/meta/likely.hpp>
#include <stan/math/rev/mat/fun/typedefs.hpp>
#include <stan/math/rev/scal/meta/operands_and_partials.hpp>
#include <vector>

namespace stan {
  namespace math {
    namespace internal {
      // Vectorized Univariate
      template <typename ViewElt>
      class ops_partials_edge<ViewElt, std::vector<var> > {
      public:
        explicit ops_partials_edge(const std::vector<var>& op)
          : partials_(op.size()), partials_vec_(partials_),
            operands_(ops) {}

        typedef std::vector<ViewElt> partials_t;
        partials_t partials_;  // For univariate use-cases
        broadcast_array<partials_t> partials_vec_;  // For multivariate

      private:
        template<typename, typename, typename, typename, typename>
        friend class stan::math::operands_and_partials;
        const std::vector<var>& operands_;

        void dump_partials(double* partials) {
          for (int i = 0; i < this->partials_.size(); ++i) {
            partials[i] = this->partials_(i);
          }
        }
        void dump_operands(vari** varis) {
          for (size_t i = 0; i < this->operands_.size(); ++i) {
            varis[i] = this->operands_[i].vi_;
          }
        }
        int size() { return this->operands_.size(); }
      };

      template <typename ViewElt, int R, int C>
      class ops_partials_edge<ViewElt, Eigen::Matrix<var, R, C> > {
      public:
        explicit ops_partials_edge(const Eigen::Matrix<var, R, C>& ops)
          : partials_(ops.rows(), ops.cols()), partials_vec_(partials_),
            operands_(ops) {}

        typedef Eigen::Matrix<ViewElt, R, C> partials_t;
        partials_t partials_;  // For univariate use-cases
        broadcast_array<partials_t> partials_vec_;  // For multivariate

      private:
        template<typename, typename, typename, typename, typename>
        friend class stan::math::operands_and_partials;
        const Op& operands_;

        void dump_operands(vari** varis) {
          for (int i = 0; i < this->operands_.size(); ++i) {
            varis[i] = this->operands_(i).vi_;
          }
        }
        void dump_partials(double* partials) {
          for (int i = 0; i < this->partials_.size(); ++i) {
            partials[i] = this->partials_(i);
          }
        }
        int size() { return this->operands_.size(); }
      };

      // SPECIALIZATIONS FOR MULTIVARIATE VECTORIZATIONS
      template <typename ViewElt, int R, int C>
      class ops_partials_edge<ViewElt, std::vector<Eigen::Matrix<var, R, C> > >
        : public ops_partials_edge_multivariate_prim
      <ViewElt, std::vector<Eigen::Matrix<var, R, C> > > {
      public:
        explicit ops_partials_edge(const Op& ops)
          : partials_vec_(ops.size()), operands_(ops) {
          for (size_t i = 0; i < ops.size(); ++i) {
            partials_vec_[i] = partial_t::Zero(ops[i].rows(), ops[i].cols());
          }
        }

        typedef std::vector<Eigen::Matrix<var, R, C> > Op;
        typedef Eigen::Matrix<ViewElt, -1, -1> partial_t;
        std::vector<partial_t> partials_vec_;

      private:
        template<typename, typename, typename, typename, typename>
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
        int size() {
          if (unlikely(this->operands_.size() == 0)) return 0;
          return this->operands_.size() * this->operands_[0].size();
        }
      };
    }  // namespace internal
  }  // namespace math
}  // namespace stan
#endif
