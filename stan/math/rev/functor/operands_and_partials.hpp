#ifndef STAN_MATH_REV_META_OPERANDS_AND_PARTIALS_HPP
#define STAN_MATH_REV_META_OPERANDS_AND_PARTIALS_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core/chainablestack.hpp>
#include <stan/math/rev/core/precomputed_gradients.hpp>
#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/core/vari.hpp>
#include <stan/math/rev/core/typedefs.hpp>
#include <stan/math/prim/meta/is_eigen.hpp>
#include <stan/math/prim/meta/is_vector_like.hpp>
#include <stan/math/prim/meta/compiler_attributes.hpp>
#include <stan/math/prim/meta/promote_scalar_type.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/functor/broadcast_array.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>
#include <vector>
#include <tuple>

namespace stan {
namespace math {

namespace internal {

/** \ingroup type_trait
 * \callergraph
 */
template <typename ViewElt, typename Op>
class ops_partials_edge<ViewElt, Op, require_var_t<Op>> {
 public:
  double partial_;
  explicit ops_partials_edge(const var& op) noexcept
      : partial_(0), operand_(op) {}

 private:
  template <typename, typename, typename...>
  friend class stan::math::operands_and_partials_impl;
  const var& operand_;

  void dump_partials(double* partials) const noexcept { *partials = this->partial_; }
  void dump_operands(vari** varis) const noexcept { *varis = this->operand_.vi_; }
  static constexpr int size() const { return 1; }
  static constexpr std::tuple<> container_operands() { return std::tuple<>(); }
  static constexpr std::tuple<> container_partials() { return std::tuple<>(); }
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
  std::tuple<internal::ops_partials_edge<double, std::decay_t<Ops>>...> edges_;
  template <typename Args...>
  explicit operands_and_partials_impl(Args&&... args) : edges_(std::forward_as_tuple(std::forward<Args>(args)...)) {}

  template <int Idx>
  inline decltype(auto) edge() noexcept {
    return std::get<Idx>(edges_);
  }

  /** \ingroup type_trait
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
  var build(double value) {
    auto edge_sizes = apply([](auto&&...args) {
      return std::array<int, std::tuple_size<decltype(edges_)>::value>{args.size()...};
    }, edges_);
    for (int i = 1; i < edge_sizes.size(); ++i) {
      total_var_size += edge_sizes[i];
      edge_sizes[i] += edge_sizes[i - 1];
    }
    const auto total_var_size = edge_sizes[edge_sizes.size() - 1];
    auto var_operands = apply([](auto&... args) {
      return std::make_tuple(args.get_vari_ops()...);
    });
    auto var_partials = apply([](auto&... args) {
      return std::make_tuple(args.get_partials()...);
    });
    auto container_operands = apply([](auto&&... edges) {
      return get_operands_tuple(to_arena(edges.container_operands()...));
    }, edges_);
    auto container_partials = apply([](auto&&... edges) {
      return get_partials_tuple(to_arena(edges.container_partials()...));
    }, edges_);
    return make_callback_var(value,
      [ret, var_operands, var_partials, total_var_size, container_operands, container_partials]() mutable {
        // Does the equivalent of
        /*
         *        for (size_t i = 0; i < size_; ++i) {
         *          varis[i]->adj_ += ret.adj_ * gradients_[i];
         *       }
         *
         */
         // For each operant / partial
        for_each([](auto& operand, auto& partial) mutable {
          accumulate_adjoint(ret, operand, partial)
        }, container_operands, container_partials);
        for_each([](auto& operand, auto& partial) mutable {
          accumulate_adjoint(ret, operand, partial)
        }, container_operands, container_partials);
      });
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
  explicit ops_partials_edge(const Op& op)
      : partials_(partials_t::Zero(op.size())),
        operands_(op) {}

 private:
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
  std::tuple<> container_operands() { return std::tuple<>(); }
  std::tuple<> container_partials() { return std::tuple<>(); }
};

template <typename Op>
class ops_partials_edge<double, Op, require_eigen_st<is_var, Op>> {
 public:
  using partials_t = promote_scalar_t<double, Op>;
  partials_t partials_;                       // For univariate use-cases
  explicit ops_partials_edge(const Op& ops)
      : partials_(partials_t::Zero(ops.rows(), ops.cols())),
        operands_(ops) {}

 private:
  template <typename, typename, typename...>
  friend class stan::math::operands_and_partials_impl;
  const Op& operands_;

  void dump_operands(vari** varis) {
    Eigen::Map<promote_scalar_t<vari*, Op>>(varis, this->operands_.rows(),
                                            this->operands_.cols())
        = this->operands_.vi();
  }
  void dump_partials(double* partials) {
    Eigen::Map<partials_t>(partials, this->partials_.rows(),
                           this->partials_.cols())
        = this->partials_;
  }
  int size() { return this->operands_.size(); }
  std::tuple<> container_operands() { return std::tuple<>(); }
  std::tuple<> container_partials() { return std::tuple<>(); }
};

template <typename Op>
class ops_partials_edge<double, var_value<Op>, require_eigen_t<Op>> {
 public:
  using partials_t = arena_t<Op>;
  partials_t partials_;                       // For univariate use-cases
  explicit ops_partials_edge(const var_value<Op>& ops)
      : partials_(
            plain_type_t<partials_t>::Zero(ops.vi_->rows(), ops.vi_->cols())),
        operands_(ops) {}

 private:
  template <typename, typename, typename...>
  friend class stan::math::operands_and_partials_impl;
  var_value<Op> operands_;

  void dump_operands(vari** varis) {}
  void dump_partials(double* partials) {}
  int size() { return 0; }
  var_value<Op> container_operands() {
    return operands_;
  }
  partials_t& container_partials() {
    return partials_;
  }
};

// SPECIALIZATIONS FOR MULTIVARIATE VECTORIZATIONS
// (i.e. nested containers)
template <int R, int C>
class ops_partials_edge<double, std::vector<Eigen::Matrix<var, R, C>>> {
 public:
  using Op = arena_t<std::vector<arena_t<Eigen::Matrix<var, R, C>>>>;
  using partial_t = arena_t<Eigen::Matrix<double, R, C>>;
  std::vector<partial_t> partials_vec_;
  explicit ops_partials_edge(const Op& ops)
      : partials_vec_(ops.size()), operands_(ops) {
    for (size_t i = 0; i < ops.size(); ++i) {
      partials_vec_[i] = partial_t::Zero(ops[i].rows(), ops[i].cols());
    }
  }

 private:
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
  std::tuple<> container_operands() { return std::tuple<>(); }
  std::tuple<> container_partials() { return std::tuple<>(); }
};

template <>
class ops_partials_edge<double, std::vector<std::vector<var>>> {
 public:
  using Op = arena_t<std::vector<arena_t<std::vector<var>>>>;
  using partial_t = arena_t<std::vector<double>>;
  std::vector<partial_t> partials_vec_;
  explicit ops_partials_edge(const Op& ops)
      : partials_vec_(stan::math::size(ops)), operands_(ops) {
    for (size_t i = 0; i < stan::math::size(ops); ++i) {
      partials_vec_[i] = partial_t(stan::math::size(ops[i]), 0.0);
    }
  }

 private:
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
  std::tuple<> container_operands() { return std::tuple<>(); }
  std::tuple<> container_partials() { return std::tuple<>(); }
};

template <typename Op>
class ops_partials_edge<double, std::vector<var_value<Op>>,
                        require_eigen_t<Op>> {
 public:
  using partials_t = arena_t<std::vector<Op>>;
  partials_t partials_vec_;
  explicit ops_partials_edge(const std::vector<var_value<Op>>& ops)
      : partials_vec_(ops.size()), operands_(ops) {
    for (size_t i = 0; i < ops.size(); ++i) {
      partials_vec_[i]
          = plain_type_t<Op>::Zero(ops[i].vi_->rows(), ops[i].vi_->cols());
    }
  }

 private:
  template <typename, typename, typename...>
  friend class stan::math::operands_and_partials_impl;
  const std::vector<var_value<Op>>& operands_;

  void dump_operands(vari** varis) {}
  void dump_partials(double* partials) {}
  int size() { return 0; }
  std::tuple<const std::vector<var_value<Op>>&> container_operands() {
    return std::forward_as_tuple(operands_);
  }
  std::tuple<partials_t&> container_partials() {
    return std::forward_as_tuple(partials_vec_);
  }
};
}  // namespace internal
}  // namespace math
}  // namespace stan
#endif
