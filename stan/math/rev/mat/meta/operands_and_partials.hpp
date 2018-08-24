#ifndef STAN_MATH_REV_MAT_META_OPERANDS_AND_PARTIALS_HPP
#define STAN_MATH_REV_MAT_META_OPERANDS_AND_PARTIALS_HPP

#include <stan/math/prim/scal/meta/broadcast_array.hpp>
#include <stan/math/prim/scal/meta/likely.hpp>
#include <stan/math/rev/core/chainablestack.hpp>
#include <stan/math/rev/core/allocator.hpp>
#include <stan/math/rev/mat/fun/typedefs.hpp>
#include <stan/math/rev/scal/meta/operands_and_partials.hpp>
#include <stan/math/prim/arr/meta/length.hpp>
#include <vector>

namespace stan {
namespace math {
namespace internal {
// Vectorized Univariate
template <>
class ops_partials_edge<double, std::vector<var> > {
 public:
  using partials_t = Eigen::Map<Eigen::VectorXd>;
  partials_t partials_;                       // For univariate use-cases
  broadcast_array<partials_t> partials_vec_;  // For multivariate
  explicit ops_partials_edge(const std::vector<var>& op)
      // it would be rude to use the move iterator to move them from op right?
      : partials_(
            ChainableStack::instance().memalloc_.calloc_array<double>(op.size()),
            op.size()),
        partials_vec_(partials_),
        operands_(op.size()) {
    for (size_t i = 0; i < op.size(); ++i) {
      operands_[i] = op[i].vi_;
    }
  }

 private:
  template <typename, typename, typename, typename, typename, typename>
  friend class stan::math::operands_and_partials;
  std::vector<vari*, autodiff_allocator<vari*>> operands_;
  void chain(double adj) {
    for (size_t i = 0; i < operands_.size(); ++i) {
      operands_[i]->adj_ += adj * partials_[i];
    }
  }
};

template <int R, int C>
class ops_partials_edge<double, Eigen::Matrix<var, R, C>> {
 public:
  typedef Eigen::Map<Eigen::Matrix<double, R, C>> partials_t;
  partials_t partials_;                       // For univariate use-cases
  broadcast_array<partials_t> partials_vec_;  // For multivariate
  explicit ops_partials_edge(const Eigen::Matrix<var, R, C> ops)
      : partials_(
            ChainableStack::instance().memalloc_.calloc_array<double>(ops.size()),
            ops.rows(), ops.cols()),
        partials_vec_(partials_), operands_(ops) {}

 private:
  template <typename, typename, typename, typename, typename, typename>
  friend class stan::math::operands_and_partials;
  Eigen::Matrix<var, R, C> operands_;
  void chain(double adj) {
    for (int i = 0; i < operands_.size(); ++i) {
      operands_(i)->adj_ += adj * partials_(i);
    }
  }
};

// SPECIALIZATIONS FOR MULTIVARIATE VECTORIZATIONS
// (i.e. nested containers)
template <int R, int C>
class ops_partials_edge<double, std::vector<Eigen::Matrix<var, R, C> > > {
 public:
  typedef std::vector<Eigen::Matrix<var, R, C> > Op;
  typedef Eigen::Map<Eigen::Matrix<double, -1, -1>> partial_t;
  std::vector<partial_t, autodiff_allocator<partial_t>> partials_vec_;
  explicit ops_partials_edge(const Op& ops) {
    for (size_t i = 0; i < ops.size(); ++i) {
      partials_vec_[i] = partial_t(ChainableStack::instance().memalloc_.calloc_array<double>(ops.size()), ops[i].rows(), ops[i].cols());
      for (int j = 0; j < ops[i].size(); ++j) {
        operands_[i][j] = ops[i](j).vi_;
      }
    }
  }

 private:
  template <typename, typename, typename, typename, typename, typename>
  friend class stan::math::operands_and_partials;
  using nested_ops_t = std::vector<vari*, autodiff_allocator<vari*>>;
  std::vector<nested_ops_t, autodiff_allocator<nested_ops_t>> operands_;
  void chain(double adj) {
    for (size_t i = 0; i < partials_vec_.size(); ++i) {
      for (int j = 0; j < partials_vec_[i].size(); ++j) {
        operands_[i][j].adj_ += adj * partials_vec_[i](j);
      }
    }
  }
};

// template <>
// class ops_partials_edge<double, std::vector<std::vector<var> > > {
//  public:
//   typedef std::vector<std::vector<var> > Op;
//   typedef std::vector<double> partial_t;
//   std::vector<partial_t> partials_vec_;
//   explicit ops_partials_edge(const Op& ops)
//       : partials_vec_(length(ops)), operands_(ops) {
//     for (size_t i = 0; i < length(ops); ++i) {
//       partials_vec_[i] = partial_t(length(ops[i]), 0.0);
//     }
//   }
//
//  private:
//   template <typename, typename, typename, typename, typename, typename>
//   friend class stan::math::operands_and_partials;
//   const Op& operands_;
//
//   void dump_partials(double* partials) {
//     int p_i = 0;
//     for (size_t i = 0; i < this->partials_vec_.size(); ++i) {
//       for (size_t j = 0; j < this->partials_vec_[i].size(); ++j, ++p_i) {
//         partials[p_i] = this->partials_vec_[i][j];
//       }
//     }
//   }
//   void dump_operands(vari** varis) {
//     int p_i = 0;
//     for (size_t i = 0; i < this->operands_.size(); ++i) {
//       for (size_t j = 0; j < this->operands_[i].size(); ++j, ++p_i) {
//         varis[p_i] = this->operands_[i][j].vi_;
//       }
//     }
//   }
//   int size() {
//     if (unlikely(this->operands_.size() == 0))
//       return 0;
//     return this->operands_.size() * this->operands_[0].size();
//   }
// };
}  // namespace internal
}  // namespace math
}  // namespace stan
#endif
