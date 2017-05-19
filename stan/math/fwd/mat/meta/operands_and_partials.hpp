#ifndef STAN_MATH_FWD_MAT_META_OPERANDS_AND_PARTIALS_HPP
#define STAN_MATH_FWD_MAT_META_OPERANDS_AND_PARTIALS_HPP

#include <stan/math/fwd/mat/fun/typedefs.hpp>
#include <stan/math/fwd/scal/meta/operands_and_partials.hpp>
#include <stan/math/prim/scal/meta/broadcast_array.hpp>
#include <vector>

namespace stan {
  namespace math {
    namespace internal {
      // Vectorized Univariate
      template <typename ViewElt, typename Dx>
      class ops_partials_edge<ViewElt, std::vector<fvar<Dx> > > {
      public:
        explicit ops_partials_edge(const std::vector<fvar<Dx> >& op)
          : partials_(op.size()), partials_vec_(partials_),
            operands_(ops) {}

        typedef std::vector<ViewElt> partials_t;
        partials_t partials_;  // For univariate use-cases
        broadcast_array<partials_t> partials_vec_;  // For multivariate

      private:
        template<typename, typename, typename, typename, typename>
        friend class stan::math::operands_and_partials;
        Dx dx() {
          Dx derivative(0);
          for (int i = 0; i < this->size(); ++i) {
            derivative += this->partials_(i) * this->operands_[i].d_;
          }
          return derivative;
        }
      };

      template <typename ViewElt, typename Dx, int R, int C>
      class ops_partials_edge<ViewElt, Eigen::Matrix<fvar<Dx>, R, C> > {
      public:
        explicit ops_partials_edge(const Eigen::Matrix<fvar<Dx>, R, C>& ops)
          : partials_(ops.rows(), ops.cols()), partials_vec_(partials_),
            operands_(ops) {}

        typedef Eigen::Matrix<ViewElt, R, C> partials_t;
        partials_t partials_;  // For univariate use-cases
        broadcast_array<partials_t> partials_vec_;  // For multivariate
      private:
        template<typename, typename, typename, typename, typename>
        friend class stan::math::operands_and_partials;
        const Op& operands_;

        Dx dx() {
          Dx derivative(0);
          for (int i = 0; i < this->size(); ++i) {
            derivative += this->partials_(i) * this->operands_(i).d_;
          }
          return derivative;
        }
      };

      // Multivariate; vectors of eigen types
      template <typename ViewElt, typename Dx, int R, int C>
      class ops_partials_edge
      <ViewElt, std::vector<Eigen::Matrix<fvar<Dx>, R, C> > > {
      public:
        explicit ops_partials_edge(const Op& ops)
          : partials_vec_(ops.size()), operands_(ops) {
          for (size_t i = 0; i < ops.size(); ++i) {
            partials_vec_[i] = partial_t::Zero(ops[i].rows(), ops[i].cols());
          }
        }

        typedef std::vector<Eigen::Matrix<fvar<Dx>, R, C> > Op;
        typedef Eigen::Matrix<ViewElt, -1, -1> partial_t;
        std::vector<partial_t> partials_vec_;
      private:
        template<typename, typename, typename, typename, typename>
        friend class stan::math::operands_and_partials;
        Dx dx() {
          Dx derivative(0);
          for (size_t i = 0; i < this->operands_.size(); ++i) {
            for (int j = 0; j < this->operands_[i].size(); ++j) {
              derivative
                += this->partials_vec_[i](j) * this->operands_[i](j).d_;
            }
          }
          return derivative;
        }
      };
    }  // namespace internal
  }  // namespace math
}  // namespace stan
#endif
