#ifndef STAN_MATH_PRIM_MAT_META_OPERANDS_AND_PARTIALS_HPP
#define STAN_MATH_PRIM_MAT_META_OPERANDS_AND_PARTIALS_HPP

#include <stan/math/prim/mat/fun/typedefs.hpp>
#include <stan/math/prim/scal/meta/operands_and_partials.hpp>
#include <stan/math/prim/scal/meta/broadcast_array.hpp>
#include <vector>

namespace stan {
  namespace math {
    namespace internal {
      // The zero_vec_or_mat metaprogram provides the logic for initializing
      // an Eigen row or column vector or matrix dynamically to zeroes of the
      // appropriate size.
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

      // This class serves as the base class for the ops_partials_edge template
      // specializations in fwd/mat and rev/mat where the operands are any
      // sort of collection, and can be used both in univariate and multivariate
      // code.
      // So this is farily overloaded, though the behavior is shared.
      // Notably, you can instantiate this with Op = std::vector<T> and either R
      // or C set to -1 and 1 - this will give us an Eigen row or column vector
      // for the partial derivatives, which will be exposed as `partials`
      // (for univariate) and as a broadcasted view in
      // `partials_vec` (multivariate).
      // This can also be instantiated more normally with Op = Eigen::Matrix
      // of some kind;
      // row or column vector or matrix. The partials (and partials_vec[0]) will
      // be initialized to the same shape as the passed in matrix, both
      // statically by sharing the R and C as well as dynamically through the
      // use of the above zero_vec_or_mat metaprogram.
      template <typename ViewElt, typename Op, int R, int C>
      class ops_partials_edge_mat {
      public:
        typedef Eigen::Matrix<ViewElt, R, C> partials_t;
        partials_t partials_;  // For univariate use-cases
        broadcast_array<partials_t> partials_vec_;  // For multivariate

        explicit ops_partials_edge_mat(const Op& ops)
          : partials_(zero_vec_or_mat<ViewElt, Op, R, C>::zero(ops)),
            partials_vec_(partials_),
            operands_(ops) {}
      protected:
        template<typename, typename, typename, typename, typename>
        friend class stan::math::operands_and_partials;
        const Op& operands_;
        void dump_partials(double* partials) {
          for (int i = 0; i < this->partials_.size(); ++i) {
            partials[i] = this->partials_(i);
          }
        }
        int size() { return this->operands_.size(); }
      };

      // This is used for purely multivariate base classes in fwd and rev.
      template <typename ViewElt, typename Ops>
      class ops_partials_edge_multivariate_prim {
      public:
        typedef Eigen::Matrix<ViewElt, -1, -1> partial_t;
        std::vector<partial_t> partials_vec_;

        explicit ops_partials_edge_multivariate_prim(const Ops& ops)
          : partials_vec_(ops.size()), operands_(ops) {
          for (size_t i = 0; i < ops.size(); ++i) {
            partials_vec_[i] = partial_t::Zero(ops[i].rows(), ops[i].cols());
          }
        }
      protected:
        template<typename, typename, typename, typename, typename>
        friend class stan::math::operands_and_partials;
        const Ops& operands_;
        int size() {
          if (this->operands_.size() == 0) return 0;
          return this->operands_.size() * this->operands_[0].size();
        }
      };
    }
  }
}
#endif
