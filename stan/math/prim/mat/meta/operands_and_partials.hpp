#ifndef STAN_MATH_PRIM_MAT_META_OPERANDS_AND_PARTIALS_HPP
#define STAN_MATH_PRIM_MAT_META_OPERANDS_AND_PARTIALS_HPP

#include <stan/math/prim/mat/fun/typedefs.hpp>
#include <stan/math/prim/scal/meta/operands_and_partials.hpp>
#include <stan/math/prim/scal/meta/broadcast_array.hpp>
#include <vector>

namespace stan {
  namespace math {
    namespace detail {
      template <typename T, typename Orig, int R, int C>
      struct zero_vec_or_mat {
        typedef Eigen::Matrix<T, R, C> ret;
        static ret zero(const Orig& in) {
          return ret::Zero(in.rows(), in.cols());
        }
      };
      template <typename T, typename Orig>
      struct zero_vec_or_mat<T, Orig, -1, 1> {
        typedef Eigen::Matrix<T, -1, 1> ret;
        static ret zero(const Orig& in) {
          return ret::Zero(in.size(), 1);
        }
      };
      template <typename T, typename Orig>
      struct zero_vec_or_mat<T, Orig, 1, -1> {
        typedef Eigen::Matrix<T, 1, -1> ret;
        static ret zero(const Orig& in) {
          return ret::Zero(1, in.size());
        }
      };

      // Handles all the vectorized cases for "views of scalars":
      // ViewElt = double, Arg = std::vector<var> d_ holds vector<double>
      // ViewElt = double, Arg = std::vector<double> d_ holds dummy
      // Op will always be some container of vars
      template <typename ViewElt, typename Op, int R, int C>
      class ops_partials_edge_mat {
      public:
        typedef Eigen::Matrix<ViewElt, R, C> partials_t;
        const Op& operands;
        partials_t partials;
        broadcast_array<partials_t> partials_vec;
        explicit ops_partials_edge_mat(const Op& ops)
          : operands(ops),
            partials(zero_vec_or_mat<ViewElt, Op, R, C>::zero(ops)),
            partials_vec(partials) {}
        void dump_partials(double* partials) {
          for (int i = 0; i < this->partials.size(); ++i) {
            partials[i] = this->partials(i);
          }
        }
        int size() {
          return this->operands.size();
        }
      };

      template <typename ViewElt, typename Op>
      class ops_partials_edge_multivariate_prim {
      public:
        typedef Eigen::Matrix<ViewElt, -1, -1> partial_t;
        std::vector<partial_t> partials_vec;
        const std::vector<Op>& operands;
        explicit ops_partials_edge_multivariate_prim(const std::vector<Op>& ops)
          : partials_vec(ops.size()), operands(ops) {
          for (size_t i = 0; i < ops.size(); ++i) {
            partials_vec[i] = partial_t::Zero(ops[i].rows(), ops[i].cols());
          }
        }
        int size() {
          return this->operands.size() * this->operands[0].size();
        }
      };
    }
  }
}
#endif
