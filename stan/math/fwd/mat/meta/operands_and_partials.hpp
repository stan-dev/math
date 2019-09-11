#ifndef STAN_MATH_FWD_MAT_META_OPERANDS_AND_PARTIALS_HPP
#define STAN_MATH_FWD_MAT_META_OPERANDS_AND_PARTIALS_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/fwd/scal/meta/is_fvar.hpp>
#include <stan/math/fwd/scal/meta/value_type.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/arr/meta/length.hpp>
#include <stan/math/prim/scal/meta/broadcast_array.hpp>
#include <stan/math/fwd/scal/meta/operands_and_partials.hpp>
#include <type_traits>
#include <vector>

namespace stan {
namespace math {
namespace internal {
// Vectorized Univariate
template <typename ViewElt, typename Op>
class ops_partials_edge<
    ViewElt, Op,
    std::enable_if_t<
        is_fvar<value_type_t<Op>>::value && is_std_vector<Op>::value
        && std::is_same<std::decay_t<ViewElt>,
                        std::decay_t<value_type_t<value_type_t<Op>>>>::value>> {
 public:
  using Dx = value_type_t<value_type_t<Op>>;
  using partials_t = Eigen::Matrix<Dx, -1, 1>;
  partials_t partials_;                       // For univariate use-cases
  broadcast_array<partials_t> partials_vec_;  // For multivariate
  explicit ops_partials_edge(const Op& ops)
      : partials_vec_(partials_),
        operands_(ops) {
          partials_ = partials_t::Zero(ops.size());
        }

 private:
  template <typename, typename, typename, typename, typename, typename>
  friend class stan::math::operands_and_partials;
  Op operands_;

  Dx dx() {
    Dx derivative(0);
    for (size_t i = 0; i < this->operands_.size(); ++i) {
      derivative += this->partials_[i] * this->operands_[i].d_;
    }
    return derivative;
  }
};

template <typename ViewElt, typename Op>
class ops_partials_edge<
    ViewElt, Op,
    std::enable_if_t<
        is_fvar<value_type_t<Op>>::value && is_eigen<Op>::value
        && std::is_same<std::decay_t<ViewElt>,
                        std::decay_t<value_type_t<value_type_t<Op>>>>::value>> {
 public:
  using Dx = value_type_t<value_type_t<Op>>;
  using partials_t = Eigen::Matrix<Dx, -1, 1>;
  partials_t partials_;                       // For univariate use-cases
  broadcast_array<partials_t> partials_vec_;  // For multivariate
  explicit ops_partials_edge(const Op& ops)
      : partials_(partials_t::Zero(ops.rows(), ops.cols())),
        partials_vec_(partials_),
        operands_(ops) {}

 private:
  template <typename, typename, typename, typename, typename, typename>
  friend class stan::math::operands_and_partials;
  Op operands_;

  Dx dx() {
    Dx derivative(0);
    for (int i = 0; i < this->operands_.eval().size(); ++i) {
      derivative += this->partials_(i) * this->operands_(i).d_;
    }
    return derivative;
  }
};

// Multivariate; vectors of eigen types
template <typename ViewElt, typename Op>
class ops_partials_edge<
    ViewElt, Op,
    std::enable_if_t<
        is_std_vector<Op>::value && is_eigen<value_type_t<Op>>::value
        && is_fvar<value_type_t<value_type_t<Op>>>::value
        && std::is_same<std::decay_t<ViewElt>,
                        std::decay_t<value_type_t<
                            value_type_t<value_type_t<Op>>>>>::value>> {
 public:
  using Dx = value_type_t<value_type_t<value_type_t<Op>>>;
  using partial_t = Eigen::Matrix<Dx, -1, 1>;
  std::vector<partial_t> partials_vec_;
  explicit ops_partials_edge(const Op& ops)
      : partials_vec_(ops.size()), operands_(ops) {
    for (size_t i = 0; i < ops.size(); ++i) {
      partials_vec_[i] = partial_t::Zero(ops[i].rows(), ops[i].cols());
    }
  }

 private:
  template <typename, typename, typename, typename, typename, typename,
            typename>
  friend class stan::math::operands_and_partials;
  Op operands_;

  Dx dx() {
    Dx derivative(0);
    for (size_t i = 0; i < this->operands_.size(); ++i) {
      for (int j = 0; j < this->operands_[i].eval().size(); ++j) {
        derivative += this->partials_vec_[i](j) * this->operands_[i](j).d_;
      }
    }
    return derivative;
  }
};

template <typename ViewElt, typename Op>
class ops_partials_edge<
    ViewElt, Op,
    std::enable_if_t<
        is_std_vector<Op>::value && is_std_vector<value_type_t<Op>>::value
        && is_fvar<value_type_t<value_type_t<Op>>>::value
        && std::is_same<std::decay_t<ViewElt>,
                        std::decay_t<value_type_t<
                            value_type_t<value_type_t<Op>>>>>::value>> {
 public:
  using Dx = value_type_t<value_type_t<value_type_t<Op>>>;
  using partial_t = std::vector<Dx>;
  std::vector<partial_t> partials_vec_;
  explicit ops_partials_edge(const Op& ops)
      : partials_vec_(length(ops)), operands_(ops) {
    for (size_t i = 0; i < length(ops); ++i) {
      partials_vec_[i] = partial_t(length(ops[i]), 0.0);
    }
  }

 private:
  template <typename, typename, typename, typename, typename, typename>
  friend class stan::math::operands_and_partials;
  Op operands_;

  Dx dx() {
    Dx derivative(0);
    for (size_t i = 0; i < this->operands_.size(); ++i) {
      for (size_t j = 0; j < this->operands_[i].size(); ++j) {
        derivative += this->partials_vec_[i][j] * this->operands_[i][j].d_;
      }
    }
    return derivative;
  }
};
}  // namespace internal
}  // namespace math
}  // namespace stan
#endif
