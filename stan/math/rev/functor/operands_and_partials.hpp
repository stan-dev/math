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
#include <stan/math/prim/functor/apply.hpp>
#include <stan/math/prim/functor/for_each.hpp>
#include <vector>
#include <tuple>

namespace stan {
namespace math {

namespace internal {

// OpenCL
template <typename T1, typename T2,
          require_all_kernel_expressions_and_none_scalar_t<T1, T2>* = nullptr>
inline void update_adjoints(var_value<T1>& x, const T2& y, const vari& z) {
  x.adj() += z.adj() * y;
}

template <typename T1, typename T2,
          require_all_kernel_expressions_and_none_scalar_t<T1, T2>* = nullptr>
inline void update_adjoints(var_value<T1>& x, const T2& y, const var& z) {
  x.adj() += z.adj() * y;
}

// Scalars
template <typename Scalar1, typename Scalar2, require_var_t<Scalar1>* = nullptr,
          require_not_var_matrix_t<Scalar1>* = nullptr,
          require_arithmetic_t<Scalar2>* = nullptr>
inline void update_adjoints(Scalar1 x, Scalar2 y, const vari& z) noexcept {
  x.adj() += z.adj() * y;
}

template <typename Scalar1, typename Scalar2, require_var_t<Scalar1>* = nullptr,
          require_not_var_matrix_t<Scalar1>* = nullptr,
          require_arithmetic_t<Scalar2>* = nullptr>
inline void update_adjoints(Scalar1 x, Scalar2 y, const var& z) noexcept {
  x.adj() += z.adj() * y;
}

// Matrix
template <typename Matrix1, typename Matrix2,
          require_rev_matrix_t<Matrix1>* = nullptr,
          require_st_arithmetic<Matrix2>* = nullptr>
inline void update_adjoints(Matrix1& x, const Matrix2& y, const vari& z) {
  x.adj().array() += z.adj() * y.array();
}

template <typename Matrix1, typename Matrix2,
          require_rev_matrix_t<Matrix1>* = nullptr,
          require_st_arithmetic<Matrix2>* = nullptr>
inline void update_adjoints(Matrix1& x, const Matrix2& y, const var& z) {
  x.adj().array() += z.adj() * y.array();
}

template <typename Arith, typename Alt, require_st_arithmetic<Arith>* = nullptr>
inline constexpr void update_adjoints(Arith&& /* x */, Alt&& /* y */,
                                      const vari& /* z */) noexcept {}

template <typename Arith, typename Alt, require_st_arithmetic<Arith>* = nullptr>
inline constexpr void update_adjoints(Arith&& /* x */, Alt&& /* y */,
                                      const var& /* z */) noexcept {}

// Vectors
template <typename StdVec1, typename Vec2,
          require_std_vector_t<StdVec1>* = nullptr,
          require_st_arithmetic<Vec2>* = nullptr>
inline void update_adjoints(StdVec1& x, const Vec2& y, const vari& z) {
  for (size_t i = 0; i < x.size(); ++i) {
    update_adjoints(x[i], y[i], z);
  }
}

template <typename StdVec1, typename Vec2,
          require_std_vector_t<StdVec1>* = nullptr,
          require_st_arithmetic<Vec2>* = nullptr>
inline void update_adjoints(StdVec1& x, const Vec2& y, const var& z) {
  for (size_t i = 0; i < x.size(); ++i) {
    update_adjoints(x[i], y[i], z);
  }
}

/** \ingroup type_trait
 * \callergraph
 */
template <>
class ops_partials_edge<double, var> {
 public:
  double partial_{0};
  broadcast_array<double> partials_{partial_};
  explicit ops_partials_edge(const var& op) noexcept : operands_(op) {}
  explicit ops_partials_edge(const ops_partials_edge<double, var>& other)
      : partial_(other.partial_),
        partials_(partial_),
        operands_(other.operands_) {}
  inline auto& partial() { return partial_; }
  inline auto& operand() noexcept { return operands_; }

  var operands_;
  static constexpr int size() { return 1; }
};
// Vectorized Univariate
// Vectorized Univariate
template <>
class ops_partials_edge<double, std::vector<var>> {
 public:
  using Op = std::vector<var, arena_allocator<var>>;
  using partials_t = arena_t<Eigen::VectorXd>;
  partials_t partials_;                       // For univariate use-cases
  broadcast_array<partials_t> partials_vec_;  // For multivariate
  explicit ops_partials_edge(const std::vector<var>& op)
      : partials_(partials_t::Zero(op.size())),
        partials_vec_(partials_),
        operands_(op.begin(), op.end()) {}

  Op operands_;

  inline int size() const noexcept { return this->operands_.size(); }
  inline auto&& operand() noexcept { return std::move(this->operands_); }
  inline auto& partial() noexcept { return this->partials_; }
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

  explicit ops_partials_edge(
      const ops_partials_edge<double, Op, require_eigen_st<is_var, Op>>& other)
      : partials_(other.partials_),
        partials_vec_(partials_),
        operands_(other.operands_) {}

  inline auto& partial() noexcept { return partials_; }
  inline auto& operand() noexcept { return operands_; }
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

  explicit ops_partials_edge(
      const ops_partials_edge<double, var_value<Op>, require_eigen_t<Op>>&
          other)
      : partials_(other.partials_),
        partials_vec_(partials_),
        operands_(other.operands_) {}

  inline auto& partial() noexcept { return partials_; }
  inline auto& operand() noexcept { return operands_; }

  var_value<Op> operands_;
  static constexpr int size() { return 0; }
};

// SPECIALIZATIONS FOR MULTIVARIATE VECTORIZATIONS
// (i.e. nested containers)
template <int R, int C>
class ops_partials_edge<double, std::vector<Eigen::Matrix<var, R, C>>> {
 public:
  using inner_op = arena_t<Eigen::Matrix<var, R, C>>;
  using Op = std::vector<inner_op, arena_allocator<inner_op>>;
  using partial_t = arena_t<Eigen::Matrix<double, R, C>>;
  std::vector<partial_t, arena_allocator<partial_t>> partials_vec_;
  explicit ops_partials_edge(const std::vector<Eigen::Matrix<var, R, C>>& ops)
      : partials_vec_(ops.size()), operands_(ops.begin(), ops.end()) {
    for (size_t i = 0; i < ops.size(); ++i) {
      partials_vec_[i] = partial_t::Zero(ops[i].rows(), ops[i].cols());
    }
  }

  Op operands_;

  inline int size() const noexcept {
    if (unlikely(this->operands_.size() == 0)) {
      return 0;
    }
    return this->operands_.size() * this->operands_[0].size();
  }
  inline auto&& operand() noexcept { return std::move(this->operands_); }
  inline auto& partial() noexcept { return this->partials_vec_; }
};

template <>
class ops_partials_edge<double, std::vector<std::vector<var>>> {
 public:
  using inner_vec = std::vector<var, arena_allocator<var>>;
  using Op = std::vector<inner_vec, arena_allocator<inner_vec>>;
  using partial_t = std::vector<double, arena_allocator<double>>;
  std::vector<partial_t, arena_allocator<partial_t>> partials_vec_;
  explicit ops_partials_edge(const std::vector<std::vector<var>>& ops)
      : partials_vec_(stan::math::size(ops)), operands_(ops.size()) {
    for (size_t i = 0; i < stan::math::size(ops); ++i) {
      operands_[i] = inner_vec(ops[i].begin(), ops[i].end());
      partials_vec_[i] = partial_t(stan::math::size(ops[i]), 0.0);
    }
  }

  Op operands_;
  inline int size() const noexcept {
    return this->operands_.size() * this->operands_[0].size();
  }
  inline auto&& operand() noexcept { return std::move(this->operands_); }
  inline auto&& partial() noexcept { return std::move(this->partials_vec_); }
};

template <typename Op>
class ops_partials_edge<double, std::vector<var_value<Op>>,
                        require_eigen_t<Op>> {
 public:
  using partials_t = std::vector<arena_t<Op>, arena_allocator<arena_t<Op>>>;
  partials_t partials_vec_;
  explicit ops_partials_edge(const std::vector<var_value<Op>>& ops)
      : partials_vec_(ops.size()), operands_(ops.begin(), ops.end()) {
    for (size_t i = 0; i < ops.size(); ++i) {
      partials_vec_[i]
          = plain_type_t<Op>::Zero(ops[i].vi_->rows(), ops[i].vi_->cols());
    }
  }

  std::vector<var_value<Op>, arena_allocator<var_value<Op>>> operands_;

  static constexpr int size() noexcept { return 0; }
  inline auto&& operand() noexcept { return std::move(this->operands_); }
  inline auto&& partial() noexcept { return std::move(this->partials_vec_); }
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
 * @tparam Op6 type of the sixth operand
 * @tparam Op7 type of the seventh operand
 * @tparam Op8 type of the eighth operand
 */
template <typename Op1, typename Op2, typename Op3, typename Op4, typename Op5,
          typename Op6, typename Op7, typename Op8>
class operands_and_partials<Op1, Op2, Op3, Op4, Op5, Op6, Op7, Op8, var> {
 public:
  internal::ops_partials_edge<double, std::decay_t<Op1>> edge1_;
  internal::ops_partials_edge<double, std::decay_t<Op2>> edge2_;
  internal::ops_partials_edge<double, std::decay_t<Op3>> edge3_;
  internal::ops_partials_edge<double, std::decay_t<Op4>> edge4_;
  internal::ops_partials_edge<double, std::decay_t<Op5>> edge5_;
  internal::ops_partials_edge<double, std::decay_t<Op6>> edge6_;
  internal::ops_partials_edge<double, std::decay_t<Op7>> edge7_;
  internal::ops_partials_edge<double, std::decay_t<Op8>> edge8_;

  explicit operands_and_partials(const Op1& o1) : edge1_(o1) {}
  operands_and_partials(const Op1& o1, const Op2& o2)
      : edge1_(o1), edge2_(o2) {}
  operands_and_partials(const Op1& o1, const Op2& o2, const Op3& o3)
      : edge1_(o1), edge2_(o2), edge3_(o3) {}
  operands_and_partials(const Op1& o1, const Op2& o2, const Op3& o3,
                        const Op4& o4)
      : edge1_(o1), edge2_(o2), edge3_(o3), edge4_(o4) {}
  operands_and_partials(const Op1& o1, const Op2& o2, const Op3& o3,
                        const Op4& o4, const Op5& o5)
      : edge1_(o1), edge2_(o2), edge3_(o3), edge4_(o4), edge5_(o5) {}
  operands_and_partials(const Op1& o1, const Op2& o2, const Op3& o3,
                        const Op4& o4, const Op5& o5, const Op6& o6)
      : edge1_(o1),
        edge2_(o2),
        edge3_(o3),
        edge4_(o4),
        edge5_(o5),
        edge6_(o6) {}
  operands_and_partials(const Op1& o1, const Op2& o2, const Op3& o3,
                        const Op4& o4, const Op5& o5, const Op6& o6,
                        const Op7& o7)
      : edge1_(o1),
        edge2_(o2),
        edge3_(o3),
        edge4_(o4),
        edge5_(o5),
        edge6_(o6),
        edge7_(o7) {}
  operands_and_partials(const Op1& o1, const Op2& o2, const Op3& o3,
                        const Op4& o4, const Op5& o5, const Op6& o6,
                        const Op7& o7, const Op8& o8)
      : edge1_(o1),
        edge2_(o2),
        edge3_(o3),
        edge4_(o4),
        edge5_(o5),
        edge6_(o6),
        edge7_(o7),
        edge8_(o8) {}

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
    return make_callback_var(
        value, [operand1 = edge1_.operand(), partial1 = edge1_.partial(),
                operand2 = edge2_.operand(), partial2 = edge2_.partial(),
                operand3 = edge3_.operand(), partial3 = edge3_.partial(),
                operand4 = edge4_.operand(), partial4 = edge4_.partial(),
                operand5 = edge5_.operand(), partial5 = edge5_.partial(),
                operand6 = edge6_.operand(), partial6 = edge6_.partial(),
                operand7 = edge7_.operand(), partial7 = edge7_.partial(),
                operand8 = edge8_.operand(),
                partial8 = edge8_.partial()](const auto& vi) mutable {
          if (!is_constant<Op1>::value) {
            internal::update_adjoints(operand1, partial1, vi);
          }
          if (!is_constant<Op2>::value) {
            internal::update_adjoints(operand2, partial2, vi);
          }
          if (!is_constant<Op3>::value) {
            internal::update_adjoints(operand3, partial3, vi);
          }
          if (!is_constant<Op4>::value) {
            internal::update_adjoints(operand4, partial4, vi);
          }
          if (!is_constant<Op5>::value) {
            internal::update_adjoints(operand5, partial5, vi);
          }
          if (!is_constant<Op6>::value) {
            internal::update_adjoints(operand6, partial6, vi);
          }
          if (!is_constant<Op7>::value) {
            internal::update_adjoints(operand7, partial7, vi);
          }
          if (!is_constant<Op8>::value) {
            internal::update_adjoints(operand8, partial8, vi);
          }
        });
  }
};

}  // namespace math
}  // namespace stan
#endif
