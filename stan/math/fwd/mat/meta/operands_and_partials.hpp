#ifndef STAN_MATH_FWD_MAT_META_OPERANDS_AND_PARTIALS_HPP
#define STAN_MATH_FWD_MAT_META_OPERANDS_AND_PARTIALS_HPP

#include <stan/math/fwd/mat/fun/typedefs.hpp>
#include <stan/math/fwd/scal/meta/operands_and_partials.hpp>
#include <stan/math/prim/mat/meta/operands_and_partials.hpp>
#include <vector>

namespace stan {
  namespace math {
    namespace detail {
      template <typename ViewElt, typename Op, typename Dx>
      class ops_partials_edge_vec
        : public ops_partials_edge_mat_prim<ViewElt, Op, -1, 1> {
      public:
        ops_partials_edge_vec(const Op& ops)
          : ops_partials_edge_mat_prim<ViewElt, Op, -1, 1>(ops) {}
        Dx dx() {
          Dx derivative(0);
          for (int i = 0; i < this->size(); ++i) {
            derivative += this->partials[i] * this->operands[i].d_;
          }
          return derivative;
        }
      };

      // Vectorized Univariate
      template <typename ViewElt, typename Dx>
      class ops_partials_edge<ViewElt, std::vector<fvar<Dx> > >
        : public ops_partials_edge_vec<ViewElt, std::vector<fvar<Dx> >, Dx> {
      public:
        explicit ops_partials_edge(const std::vector<fvar<Dx> >& ops)
          : ops_partials_edge_vec<ViewElt, std::vector<fvar<Dx> >, Dx>(ops) {}
      };
      template <typename ViewElt, typename Dx, int R, int C>
      class ops_partials_edge<ViewElt, Eigen::Matrix<fvar<Dx>, R, C> >
        : public ops_partials_edge_vec<ViewElt, Eigen::Matrix<fvar<Dx>, R, C>, Dx> {
      public:
        explicit ops_partials_edge(const Eigen::Matrix<fvar<Dx>, R, C>& ops)
          : ops_partials_edge_vec<ViewElt, Eigen::Matrix<fvar<Dx>, R, C>, Dx>(ops) {}
      };

      // Multivariate; vectors of eigen types
      template <typename ViewElt, typename Op, typename Dx>
      class ops_partials_edge_multivariate
        : public ops_partials_edge_multivariate_prim<ViewElt, Op> {
      public:
        explicit ops_partials_edge_multivariate(const std::vector<Op>& ops)
          : ops_partials_edge_multivariate_prim<ViewElt, Op>(ops) {}
        Dx dx() {
          Dx derivative(0);
          for (int i = 0; i < this->operands.size(); ++i) {
            for (int j = 0; j < this->operands[i].size(); ++j) {
              derivative += this->partials[i](j) * this->operands[i](j);
            }
          }
        }
      };
      template <typename ViewElt, typename Dx, int R, int C>
      class ops_partials_edge<ViewElt, std::vector<Eigen::Matrix<fvar<Dx>, R, C> > >
        : public ops_partials_edge_multivariate<ViewElt, Eigen::Matrix<fvar<Dx>, R, C>, Dx> {
      public:
        explicit ops_partials_edge(const std::vector<Eigen::Matrix<fvar<Dx>, R, C> >& ops)
          : ops_partials_edge_multivariate<ViewElt, Eigen::Matrix<fvar<Dx>, R, C>, Dx>(ops) {}
      };
    }
  }
}
#endif
