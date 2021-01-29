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
template <>
class ops_partials_edge<double, var> {
 public:
  double partial_{0};
  broadcast_array<double> partials_{partial_};
  explicit ops_partials_edge(const var& op) noexcept
      : operands_(op) {}

  inline auto& partial() {
    return partial_;
  }
  inline auto& operand() noexcept {
    return operands_;
  }

  var operands_;
  static constexpr int size() { return 1; }
};

template <typename Scalar1, typename Scalar2, require_var_t<Scalar1>* = nullptr>
inline void accumulate_adjoints(Scalar1&& x, Scalar2&& y, const var& z) {
  x.adj() += z.adj() * y;
}
template <typename EigT1, typename EigT2, require_rev_matrix_t<EigT1>* = nullptr>
inline void accumulate_adjoints(EigT1&& x, EigT2&& y, const var& z) {
  x.adj().array() += z.adj() * y.array();
}
template <typename Arith, typename Alt, require_st_arithmetic<Arith>* = nullptr>
inline void accumulate_adjoints(Arith&& /* x */, Alt&& /* y */, const var& z) {
}

template <typename StdVec1, typename StdVec2, require_std_vector_t<StdVec1>* = nullptr>
inline void accumulate_adjoints(StdVec1&& x, StdVec2&& y, const var& z) {
  for (size_t i = 0; i < x.size(); ++i) {
    accumulate_adjoints(x[i], y[i], z);
  }
}

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
 * @tparam Ops Type of the operands placed into the edges
 * @tparam ReturnType The type returned from the `build()` method.
 */
 template <typename ReturnType, typename... Ops>
 class operands_and_partials_impl<ReturnType, require_var_t<ReturnType>, Ops...> {
 public:
  std::tuple<internal::ops_partials_edge<double, plain_type_t<std::decay_t<Ops>>>...> edges_;
  template <typename... Types>
  explicit operands_and_partials_impl(Types&&... ops) :
    edges_(ops...) {}

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
  inline var build(double value) {
    auto operands_tuple = stan::math::apply([](auto&&... edges) {
      return std::make_tuple(edges.operand()...);
    }, edges_);
    auto partials_tuple = stan::math::apply([](auto&&... edges) {
      return std::make_tuple(edges.partial()...);
    }, edges_);
    var ret(value);
    for_each([ret](auto& operand, auto& partial) mutable {
      reverse_pass_callback([operand, partial, ret]() mutable {
        accumulate_adjoints(operand, partial, ret);
      });
    }, operands_tuple, partials_tuple);
    return ret;
  }

};

// Vectorized Univariate
template <>
class ops_partials_edge<double, std::vector<var>> {
 public:
  using Op = arena_t<std::vector<var>>;
  using partials_t = arena_t<Eigen::VectorXd>;
  partials_t partials_;                       // For univariate use-cases
  broadcast_array<partials_t> partials_vec_;  // For multivariate
  explicit ops_partials_edge(const std::vector<var>& op)
      : partials_(Eigen::VectorXd::Zero(op.size())),
        partials_vec_(partials_),
        operands_(to_arena(op)) {}
  inline auto& partial() noexcept {
    return partials_;
  }
  inline auto& operand() noexcept {
    return operands_;
  }
  inline int size() const noexcept { return this->operands_.size(); }
  Op operands_;
};

template <typename Op>
class ops_partials_edge<double, Op, require_eigen_st<is_var, Op>> {
 public:
  using partials_t = arena_t<promote_scalar_t<double, Op>>;
  partials_t partials_;                       // For univariate use-cases
  broadcast_array<partials_t> partials_vec_;  // For multivariate
  explicit ops_partials_edge(const Op& ops)
      : partials_(partials_t::Zero(ops.rows(), ops.cols())),
        partials_vec_(partials_),
        operands_(to_arena(ops)) {}

  inline auto& partial() noexcept {
    return partials_;
  }
  inline auto& operand() noexcept {
    return operands_;
  }
  arena_t<Op> operands_;
  inline auto size() const noexcept { return this->operands_.size(); }
};

template <typename Op>
class ops_partials_edge<double, var_value<Op>, require_eigen_t<Op>> {
 public:
  using partials_t = arena_t<Op>;
  partials_t partials_;                       // For univariate use-cases
  broadcast_array<partials_t> partials_vec_;  // For multivariate
  explicit ops_partials_edge(const var_value<Op>& ops)
      : partials_(
            plain_type_t<partials_t>::Zero(ops.vi_->rows(), ops.vi_->cols())),
        partials_vec_(partials_),
        operands_(ops) {}
  inline auto& partial() noexcept {
    return partials_;
  }
  inline auto& operand() noexcept {
    return operands_;
  }

  const var_value<Op>& operands_;
  static constexpr int size() { return 0; }
};

// SPECIALIZATIONS FOR MULTIVARIATE VECTORIZATIONS
// (i.e. nested containers)
template <int R, int C>
class ops_partials_edge<double, std::vector<Eigen::Matrix<var, R, C>>> {
 public:
  using Op = arena_t<std::vector<arena_t<Eigen::Matrix<var, R, C>>>>;
  using partial_t = arena_t<Eigen::Matrix<double, R, C>>;
  arena_t<std::vector<partial_t>> partials_vec_;
  explicit ops_partials_edge(const std::vector<Eigen::Matrix<var, R, C>>& ops)
      : partials_vec_(ops.size()), operands_(to_arena(ops)) {
    for (size_t i = 0; i < ops.size(); ++i) {
      partials_vec_[i] = partial_t::Zero(ops[i].rows(), ops[i].cols());
    }
  }

  inline auto& partial() noexcept {
    return partials_vec_;
  }
  inline auto& operand() noexcept {
    return operands_;
  }
  inline int size() const noexcept {
    if (unlikely(this->operands_.size() == 0)) {
      return 0;
    }
    return this->operands_.size() * this->operands_[0].size();
  }
  Op operands_;
};

template <>
class ops_partials_edge<double, std::vector<std::vector<var>>> {
 public:
  using Op = arena_t<std::vector<arena_t<std::vector<var>>>>;
  using partial_t = arena_t<Eigen::Matrix<double, -1, 1>>;
  arena_t<std::vector<partial_t>> partials_vec_;
  explicit ops_partials_edge(const std::vector<std::vector<var>>& ops)
      : partials_vec_(stan::math::size(ops)), operands_(to_arena(ops)) {
    for (size_t i = 0; i < stan::math::size(ops); ++i) {
      partials_vec_[i] = to_arena(Eigen::Matrix<double, -1, 1>(stan::math::size(ops[i])));
    }
  }
  inline auto& partial() noexcept {
    return partials_vec_;
  }

  inline auto& operand() noexcept {
    return operands_;
  }
  inline int size() const noexcept {
    if (unlikely(this->operands_.size() == 0)) {
      return 0;
    }
    return this->operands_.size() * this->operands_[0].size();
  }

  arena_t<Op> operands_;
};

template <typename Op>
class ops_partials_edge<double, std::vector<var_value<Op>>,
                        require_eigen_t<Op>> {
 public:
  using partials_t = arena_t<std::vector<Op>>;
  partials_t partials_vec_;
  explicit ops_partials_edge(const std::vector<var_value<Op>>& ops)
      : partials_vec_(ops.size()), operands_(to_arena(ops)) {
    for (size_t i = 0; i < ops.size(); ++i) {
      partials_vec_[i]
          = plain_type_t<Op>::Zero(ops[i].vi_->rows(), ops[i].vi_->cols());
    }
  }
  inline auto& partial() noexcept {
    return partials_vec_;
  }
  inline auto& operand() noexcept {
    return operands_;
  }
  static constexpr int size() noexcept { return 0; }
  arena_t<std::vector<var_value<Op>>> operands_;
};
}  // namespace internal
}  // namespace math
}  // namespace stan
#endif
