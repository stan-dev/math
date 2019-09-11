#ifndef STAN_MATH_REV_MAT_META_OPERANDS_AND_PARTIALS_HPP
#define STAN_MATH_REV_MAT_META_OPERANDS_AND_PARTIALS_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/rev/scal/meta/is_var.hpp>
#include <stan/math/rev/scal/meta/operands_and_partials.hpp>
#include <stan/math/rev/mat/fun/typedefs.hpp>
#include <stan/math/prim/scal/meta/broadcast_array.hpp>
#include <stan/math/prim/scal/meta/likely.hpp>
#include <stan/math/prim/arr/meta/length.hpp>
#include <vector>
#include <type_traits>
namespace stan {
namespace math {
namespace internal {
// Vectorized Univariate
template <typename ViewElt, typename Op>
class ops_partials_edge<
    ViewElt, Op,
    std::enable_if_t<std::is_floating_point<std::decay_t<ViewElt>>::value
                     && is_std_vector<Op>::value
                     && is_var<value_type_t<Op>>::value>> {
 public:
  using partials_t = Eigen::VectorXd;
  partials_t partials_;                       // For univariate use-cases
  broadcast_array<partials_t> partials_vec_;  // For multivariate
  explicit ops_partials_edge(const Op& op)
      : partials_(partials_t::Zero(op.size())),
        partials_vec_(partials_),
        operands_(op) {}

 private:
  template <typename, typename, typename, typename, typename, typename,
            typename>
  friend class stan::math::operands_and_partials;
  Op operands_;

  void dump_partials(ViewElt* partials) {
    for (int i = 0; i < this->partials_.size(); ++i) {
      partials[i] = this->partials_[i];
    }
  }
  void dump_operands(vari** varis) {
    for (size_t i = 0; i < this->operands_.size(); ++i) {
      varis[i] = this->operands_[i].vi_;
    }
  }
  int size() { return this->operands_.size(); }
};

template <typename ViewElt, typename Op>
class ops_partials_edge<
    ViewElt, Op,
    std::enable_if_t<std::is_floating_point<std::decay_t<ViewElt>>::value
                     && is_eigen<Op>::value
                     && is_var<value_type_t<Op>>::value>> {
 public:
  using partials_t = Eigen::Matrix<ViewElt, Op::RowsAtCompileTime, Op::ColsAtCompileTime>;
  partials_t partials_;                       // For univariate use-cases
  broadcast_array<partials_t> partials_vec_;  // For multivariate
  explicit ops_partials_edge(const Op& ops)
      : partials_(partials_t::Zero(ops.rows(), ops.cols())),
        partials_vec_(partials_),
        operands_(ops) {}

 private:
  template <typename, typename, typename, typename, typename, typename,
            typename>
  friend class stan::math::operands_and_partials;
  Op operands_;

  void dump_operands(vari** varis) {
    for (int i = 0; i < this->operands_.size(); ++i) {
      varis[i] = this->operands_(i).vi_;
    }
  }
  void dump_partials(ViewElt* partials) {
    for (int i = 0; i < this->partials_.size(); ++i) {
      partials[i] = this->partials_(i);
    }
  }
  int size() { return this->operands_.size(); }
};

// SPECIALIZATIONS FOR MULTIVARIATE VECTORIZATIONS
// (i.e. nested containers)
template <typename ViewElt, typename Op>
class ops_partials_edge<
    ViewElt, Op,
    std::enable_if_t<std::is_floating_point<std::decay_t<ViewElt>>::value
                     && is_std_vector<Op>::value
                     && is_eigen<value_type_t<Op>>::value
                     && is_var<value_type_t<value_type_t<Op>>>::value>> {
 public:
  using inner_type = value_type_t<Op>;
  using partial_t = Eigen::Matrix<ViewElt, inner_type::RowsAtCompileTime, inner_type::ColsAtCompileTime>;
  std::vector<partial_t> partials_vec_;
  explicit ops_partials_edge(const Op& ops)
      : partials_vec_(ops.size()), operands_(ops) {
    for (size_t i = 0; i < ops.size(); ++i) {
      partials_vec_[i] = partial_t::Zero(ops[i].rows(), ops[i].cols()).eval();
    }
  }

 private:
  template <typename, typename, typename, typename, typename, typename,
            typename>
  friend class stan::math::operands_and_partials;
  Op operands_;

  void dump_partials(ViewElt* partials) {
    int p_i = 0;
    for (size_t i = 0; i < this->partials_vec_.size(); ++i) {
      for (int j = 0; j < this->partials_vec_[i].eval().size(); ++j, ++p_i) {
        partials[p_i] = this->partials_vec_[i].eval()(j);
      }
    }
  }
  void dump_operands(vari** varis) {
    int p_i = 0;
    for (size_t i = 0; i < this->operands_.size(); ++i) {
      for (int j = 0; j < this->operands_[i].eval().size(); ++j, ++p_i) {
        varis[p_i] = this->operands_[i].eval()(j).vi_;
      }
    }
  }
  int size() {
    if (unlikely(this->operands_.size() == 0)) {
      return 0;
    }
    return this->operands_.size() * this->operands_[0].size();
  }
};

template <typename ViewElt, typename Op>
class ops_partials_edge<
    ViewElt, Op,
    std::enable_if_t<
        std::is_floating_point<std::decay_t<ViewElt>>::value
        && is_std_vector<Op>::value
        && is_std_vector<value_type_t<Op>>::value
        && is_var<value_type_t<value_type_t<Op>>>::value>> {
 public:
  using partial_t = std::vector<ViewElt>;
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

  void dump_partials(ViewElt* partials) {
    int p_i = 0;
    for (size_t i = 0; i < this->partials_vec_.size(); ++i) {
      for (size_t j = 0; j < this->partials_vec_[i].size(); ++j, ++p_i) {
        partials[p_i] = this->partials_vec_[i][j];
      }
    }
  }
  void dump_operands(vari** varis) {
    int p_i = 0;
    for (size_t i = 0; i < this->operands_.size(); ++i) {
      for (size_t j = 0; j < this->operands_[i].size(); ++j, ++p_i) {
        varis[p_i] = this->operands_[i][j].vi_;
      }
    }
  }
  int size() {
    if (unlikely(this->operands_.size() == 0)) {
      return 0;
    }
    return this->operands_.size() * this->operands_[0].size();
  }
};
}  // namespace internal
}  // namespace math
}  // namespace stan
#endif
