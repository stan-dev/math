#ifndef STAN_MATH_REV_META_OPERANDS_AND_PARTIALS_HPP
#define STAN_MATH_REV_META_OPERANDS_AND_PARTIALS_HPP

#include <stan/math/rev/core/chainablestack.hpp>
#include <stan/math/rev/core/precomputed_gradients.hpp>
#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/core/vari.hpp>
#include <stan/math/rev/fun/typedefs.hpp>
#include <stan/math/prim/meta/require_generics.hpp>
#include <stan/math/prim/meta/broadcast_array.hpp>
#include <stan/math/prim/meta/operands_and_partials.hpp>
#include <stan/math/prim/meta/is_vector_like.hpp>
#include <stan/math/prim/meta/likely.hpp>
#include <stan/math/prim/meta/promote_scalar_type.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <vector>
#include <iostream>


namespace stan {
namespace math {
namespace internal {

/** \ingroup type_trait
 * \callergraph
 */
template <>
class ops_partials_edge<double, var> {
 public:
  double partial_{0};
  broadcast_array<double> partials_{partial_};
  explicit ops_partials_edge(const var& op):
    partial_(0), partials_(partial_), operand_(op) {}
 private:
   template <typename...>
   friend class stan::math::operands_and_partials;
   template <typename, typename, typename...>
   friend class stan::math::operands_and_partials_impl;
  const var& operand_;

  void dump_partials(double* partials) { partials[0] = this->partial_; }
  void dump_operands(vari** varis) { varis[0] = this->operand_.vi_; }
  int size() const { return 1; }
};
}  // namespace internal

/** \ingroup type_trait
 * \callergraph
 * This class builds partial derivatives with respect to a set of
 * operands. There are two reason for the generality of this
 * class. The first is to handle vector and scalar arguments
 * without needing to write additional code. The second is to use
 * this class for writing probability distributions that handle
 * primitives, reverse mode, and forward mode variables
 * seamlessly.
 *
 * Conceptually, this class is used when we want to manually calculate
 * the derivative of a function and store this manual result on the
 * autodiff stack in a sort of "compressed" form. Think of it like an
 * easy-to-use interface to rev/core/precomputed_gradients.
 *
 * This class now supports multivariate use-cases as well by
 * exposing edge#_.partials_vec
 *
 * This is the specialization for when the return type is var,
 * which should be for all of the reverse mode cases.
 *
 * NB: since ops_partials_edge.partials_ and ops_partials_edge.partials_vec
 * are sometimes represented internally as a broadcast_array, we need to take
 * care with assignments to them. Indeed, we can assign any right hand side
 * which allows for indexing to a broadcast_array. The resulting behaviour is
 * that the entry for the first index is what gets assigned. The most common
 * use-case should be where the rhs is some container of length 1.
 *
 * @tparam Op1 type of the first operand
 * @tparam Op2 type of the second operand
 * @tparam Op3 type of the third operand
 * @tparam Op4 type of the fourth operand
 * @tparam Op5 type of the fifth operand
 */
 template <typename ReturnType, typename... Ops>
 class operands_and_partials_impl<ReturnType, require_var_t<ReturnType>, Ops...> {
 public:
  std::tuple<internal::ops_partials_edge<double, Ops>...> edges_;
  template <int id>
  inline auto& edge() {
    return std::get<id - 1>(edges_);
  }
  explicit operands_and_partials_impl(const Ops&... ops) :
   edges_(ops...) {}

  /**
   * Build the node to be stored on the autodiff graph.
   * This should contain both the value and the tangent.
   *
   * For scalars, we don't calculate any tangents.
   * For reverse mode, we end up returning a type of var that will calculate
   * the appropriate adjoint using the stored operands and partials.
   * Forward mode just calculates the tangent on the spot and returns it in
   * a vanilla fvar.
   *
   * @param value the return value of the function we are compressing
   * @return the node to be stored in the expression graph for autodiff
   */
  template <typename T, typename... Args>
  inline void dump_operands_and_partials(int idx, vari** vari_stuff, double* partials_stuff, T& op1, Args&&... args) {
    op1.dump_operands(&vari_stuff[idx]);
    op1.dump_partials(&partials_stuff[idx]);
    dump_operands_and_partials(idx + op1.size(), vari_stuff, partials_stuff, args...);
  }
  inline void dump_operands_and_partials(int idx, vari** vari_stuff, double* partials_stuff) {}

  inline var build(double value) {
    size_t edges_size = apply(
    [&](auto&&... args) {
      return sum(args.size()...);
    },
    this->edges_);
    vari** varis
        = ChainableStack::instance_->memalloc_.alloc_array<vari*>(edges_size);
    double* partials
        = ChainableStack::instance_->memalloc_.alloc_array<double>(edges_size);
    apply([&](auto&&... args) {
      dump_operands_and_partials(0, varis, partials, args...);
    }, this->edges_);
    return var(
        new precomputed_gradients_vari(value, edges_size, varis, partials));
  }
};

namespace internal {
// Vectorized Univariate
template <>
class ops_partials_edge<double, std::vector<var>> {
 public:
  using Op = std::vector<var>;
  using partials_t = Eigen::VectorXd;
  partials_t partials_;                       // For univariate use-cases
  broadcast_array<partials_t> partials_vec_{partials_};  // For multivariate
  explicit ops_partials_edge(const Op& op)
      : partials_(partials_t::Zero(op.size())),
      partials_vec_(partials_),
        operands_(op) {}

 private:
   template <typename...>
   friend class stan::math::operands_and_partials;
   template <typename, typename, typename...>
   friend class stan::math::operands_and_partials_impl;
  const Op& operands_;

  void dump_partials(double* partials) {
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

template <typename Op>
class ops_partials_edge<double, Op, require_eigen_st<is_var, Op>> {
 public:
  using partials_t = promote_scalar_t<double, Op>;
  partials_t partials_;                       // For univariate use-cases
  broadcast_array<partials_t> partials_vec_{partials_};  // For multivariate
  explicit ops_partials_edge(const Op& ops)
      : partials_(partials_t::Zero(ops.rows(), ops.cols())),
        partials_vec_(partials_),
        operands_(ops) {
        }

 private:
   template <typename...>
   friend class stan::math::operands_and_partials;
   template <typename, typename, typename...>
   friend class stan::math::operands_and_partials_impl;
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
// (i.e. nested containers)
template <int R, int C>
class ops_partials_edge<double, std::vector<Eigen::Matrix<var, R, C>>> {
 public:
  using Op = std::vector<Eigen::Matrix<var, R, C>>;
  using partial_t = Eigen::Matrix<double, -1, -1>;
  std::vector<partial_t> partials_vec_;
  explicit ops_partials_edge(const Op& ops)
      : partials_vec_(ops.size()), operands_(ops) {
    for (size_t i = 0; i < ops.size(); ++i) {
      partials_vec_[i] = partial_t::Zero(ops[i].rows(), ops[i].cols());
    }
  }

 private:
   template <typename...>
   friend class stan::math::operands_and_partials;
   template <typename, typename, typename...>
   friend class stan::math::operands_and_partials_impl;
  const Op& operands_;

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
    if (unlikely(this->operands_.size() == 0)) {
      return 0;
    }
    return this->operands_.size() * this->operands_[0].size();
  }
};

template <>
class ops_partials_edge<double, std::vector<std::vector<var>>> {
 public:
  using Op = std::vector<std::vector<var>>;
  using partial_t = std::vector<double>;
  std::vector<partial_t> partials_vec_;
  explicit ops_partials_edge(const Op& ops)
      : partials_vec_(stan::math::size(ops)), operands_(ops) {
    for (size_t i = 0; i < stan::math::size(ops); ++i) {
      partials_vec_[i] = partial_t(stan::math::size(ops[i]), 0.0);
    }
  }

 private:
   template <typename...>
   friend class stan::math::operands_and_partials;
   template <typename, typename, typename...>
   friend class stan::math::operands_and_partials_impl;
  const Op& operands_;

  void dump_partials(double* partials) {
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
