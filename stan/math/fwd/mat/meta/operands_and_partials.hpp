#ifndef STAN_MATH_FWD_MAT_META_OPERANDS_AND_PARTIALS_HPP
#define STAN_MATH_FWD_MAT_META_OPERANDS_AND_PARTIALS_HPP

#include <stan/math/fwd/mat/fun/typedefs.hpp>
#include <stan/math/fwd/scal/meta/operands_and_partials.hpp>
#include <stan/math/prim/mat/meta/operands_and_partials.hpp>
#include <vector>

namespace stan {
  namespace math {
    namespace detail {
      // Vectorized Univariate
      template <typename ViewElt, typename Dx>
      class ops_partials_edge<ViewElt, std::vector<fvar<Dx> > >
        : public ops_partials_edge_mat<ViewElt,
                                       std::vector<fvar<Dx> >, -1, 1> {
      public:
        explicit ops_partials_edge(const std::vector<fvar<Dx> >& ops)
          : ops_partials_edge_mat<ViewElt,
                                  std::vector<fvar<Dx> >, -1, 1>(ops) {}
        Dx dx() {
          Dx derivative(0);
          for (int i = 0; i < this->size(); ++i) {
            derivative += this->partials_(i) * this->operands_[i].d_;
          }
          return derivative;
        }
      };

      template <typename ViewElt, typename Dx, int R, int C>
      class ops_partials_edge<ViewElt, Eigen::Matrix<fvar<Dx>, R, C> >
        : public ops_partials_edge_mat<ViewElt,
                                       Eigen::Matrix<fvar<Dx>, R, C>, R, C> {
      public:
        explicit ops_partials_edge(const Eigen::Matrix<fvar<Dx>, R, C>& ops)
          : ops_partials_edge_mat<ViewElt,
                                  Eigen::Matrix<fvar<Dx>, R, C>, R, C>(ops) {}
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
      class ops_partials_edge<ViewElt,
                              std::vector<Eigen::Matrix<fvar<Dx>, R, C> > >
        : public ops_partials_edge_multivariate_prim
      <ViewElt, std::vector<Eigen::Matrix<fvar<Dx>, R, C> > > {
      public:
        typedef std::vector<Eigen::Matrix<fvar<Dx>, R, C> > Op;
        explicit ops_partials_edge(const Op& ops)
          : ops_partials_edge_multivariate_prim<ViewElt, Op>(ops) {}
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
    }  // namespace detail
  }  // namespace math
}  // namespace stan
#endif
