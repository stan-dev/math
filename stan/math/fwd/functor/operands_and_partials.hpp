#ifndef STAN_MATH_FWD_META_OPERANDS_AND_PARTIALS_HPP
#define STAN_MATH_FWD_META_OPERANDS_AND_PARTIALS_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/fwd/core/fvar.hpp>
#include <stan/math/fwd/meta.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/elt_multiply.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/functor/broadcast_array.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>
#include <stan/math/prim/functor/apply.hpp>
#include <stan/math/prim/functor/for_each.hpp>
#include <vector>

namespace stan {
namespace math {
namespace internal {

/**
 * End of recursion for summing `.dx()` for `fvar<T>` ops and partials.
 */
static constexpr auto sum_dx() { return static_cast<double>(0.0); }

/**
 * End of recursion for summing `.dx()` for `fvar<T>` ops and partials.
 * @tparam T1 a type with a `.dx()` method.
 * @param a an edge from `operands_and_partials_impl`
 */
template <typename T1>
inline auto sum_dx(T1& a) {
  return a.dx();
}

/**
 * Accumulate the `.dx()` from each of the edges in `ops_and_partials_impl`.
 * @tparam T1 a type with a `.dx()` method.
 * @tparam T2 a type with a `.dx()` method.
 * @tparam Types a type with a `.dx()` method.
 * @param a an edge from `operands_and_partials_impl`
 * @param b an edge from `operands_and_partials_impl`
 * @param args edges from `operands_and_partials_impl`
 */
template <typename T1, typename T2, typename... Types>
inline auto sum_dx(T1& a, T2& b, Types&... args) {
  return a.dx() + b.dx() + sum_dx(args...);
}

template <typename InnerType, typename T>
class ops_partials_edge<InnerType, T, require_fvar_t<T>> {
 public:
  using Op = std::decay_t<T>;
  using Dx = std::decay_t<InnerType>;
  Dx partial_{0};
  broadcast_array<Dx> partials_{partial_};

  explicit ops_partials_edge(const T& op)
      : partial_(0), partials_(partial_), operands_(op) {}

  explicit ops_partials_edge(
      const ops_partials_edge<InnerType, T, require_fvar_t<T>>& other)
      : partial_(other.partial_),
        partials_(partial_),
        operands_(other.operands_) {}

  explicit ops_partials_edge(
      ops_partials_edge<InnerType, T, require_fvar_t<T>>&& other)
      : partial_(other.partial_),
        partials_(partial_),
        operands_(other.operands_) {}

  const Op& operands_;

  inline Dx dx() { return this->partial_ * this->operands_.d_; }
};

// Vectorized Univariate
template <typename InnerType, typename T>
class ops_partials_edge<InnerType, T, require_std_vector_vt<is_fvar, T>> {
 public:
  using Op = std::decay_t<T>;
  using Dx = std::decay_t<InnerType>;
  using partials_t = Eigen::Matrix<Dx, -1, 1>;
  partials_t partials_;  // For univariate use-cases
  broadcast_array<partials_t> partials_vec_{partials_};  // For multivariate
  explicit ops_partials_edge(const Op& ops)
      : partials_(partials_t::Zero(ops.size()).eval()), operands_(ops) {}

  explicit ops_partials_edge(
      const ops_partials_edge<InnerType, T, require_std_vector_vt<is_fvar, T>>&
          other)
      : partials_(other.partials_),
        partials_vec_(partials_),
        operands_(other.operands_) {}

  explicit ops_partials_edge(
      ops_partials_edge<InnerType, T, require_std_vector_vt<is_fvar, T>>&&
          other)
      : partials_(std::move(other.partials_)),
        partials_vec_(partials_),
        operands_(other.operands_) {}

  const Op& operands_;
  inline Dx dx() {
    return dot_product(as_column_vector_or_scalar(this->partials_),
                       as_column_vector_or_scalar(this->operands_).d());
  }
};

template <typename Dx, typename ViewElt>
class ops_partials_edge<Dx, ViewElt, require_eigen_vt<is_fvar, ViewElt>> {
 public:
  using partials_t = promote_scalar_t<Dx, ViewElt>;
  using Op = plain_type_t<ViewElt>;
  partials_t partials_;  // For univariate use-cases
  broadcast_array<partials_t> partials_vec_{partials_};  // For multivariate
  template <typename OpT, require_eigen_vt<is_fvar, OpT>* = nullptr>
  explicit ops_partials_edge(const OpT& ops)
      : partials_(partials_t::Zero(ops.rows(), ops.cols())), operands_(ops) {}

  explicit ops_partials_edge(
      const ops_partials_edge<Dx, ViewElt, require_eigen_vt<is_fvar, ViewElt>>&
          other)
      : partials_(other.partials_),
        partials_vec_(partials_),
        operands_(other.operands_) {}

  explicit ops_partials_edge(
      ops_partials_edge<Dx, ViewElt, require_eigen_vt<is_fvar, ViewElt>>&&
          other)
      : partials_(std::move(other.partials_)),
        partials_vec_(partials_),
        operands_(other.operands_) {}

  Op operands_;

  inline Dx dx() {
    return sum(elt_multiply(this->partials_, this->operands_.d()));
  }
};

// Multivariate; vectors of eigen types
template <typename Dx, int R, int C>
class ops_partials_edge<Dx, std::vector<Eigen::Matrix<fvar<Dx>, R, C>>> {
 public:
  using Op = std::vector<Eigen::Matrix<fvar<Dx>, R, C>>;
  using partial_t = Eigen::Matrix<Dx, R, C>;
  std::vector<partial_t> partials_vec_;
  explicit ops_partials_edge(const Op& ops)
      : partials_vec_(ops.size()), operands_(ops) {
    for (size_t i = 0; i < ops.size(); ++i) {
      partials_vec_[i] = partial_t::Zero(ops[i].rows(), ops[i].cols());
    }
  }

  Op operands_;

  inline Dx dx() {
    Dx derivative(0);
    for (size_t i = 0; i < this->operands_.size(); ++i) {
      derivative
          += sum(elt_multiply(this->partials_vec_[i], this->operands_[i].d()));
    }
    return derivative;
  }
};

template <typename Dx>
class ops_partials_edge<Dx, std::vector<std::vector<fvar<Dx>>>> {
 public:
  using Op = std::vector<std::vector<fvar<Dx>>>;
  using partial_t = std::vector<Dx>;
  std::vector<partial_t> partials_vec_;
  explicit ops_partials_edge(const Op& ops)
      : partials_vec_(stan::math::size(ops)), operands_(ops) {
    for (size_t i = 0; i < stan::math::size(ops); ++i) {
      partials_vec_[i] = partial_t(stan::math::size(ops[i]), 0.0);
    }
  }

  Op operands_;
  inline Dx dx() {
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

/** \ingroup type_trait
 * This class builds partial derivatives with respect to a set of
 * operands. There are two reason for the generality of this
 * class. The first is to handle vector and scalar arguments
 * without needing to write additional code. The second is to use
 * this class for writing probability distributions that handle
 * primitives, reverse mode, and forward mode variables
 * seamlessly.
 *
 * Conceptually, this class is used when we want to calculate manually
 * the derivative of a function and store this manual result on the
 * autodiff stack in a sort of "compressed" form. Think of it like an
 * easy-to-use interface to rev/core/precomputed_gradients.
 *
 * This class now supports multivariate use-cases as well by
 * exposing edge#_.partials_vec
 *
 * This is the specialization for when the return type is fvar,
 * which should be for forward mode and all higher-order cases.
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
 * @tparam T_return_type return type of the expression. This defaults
 *   to a template metaprogram that calculates the scalar promotion of
 *   Op1 -- Op8
 */
template <typename Op1, typename Op2, typename Op3, typename Op4, typename Op5,
          typename Op6, typename Op7, typename Op8, typename Dx>
class operands_and_partials<Op1, Op2, Op3, Op4, Op5, Op6, Op7, Op8, fvar<Dx>> {
 public:
  internal::ops_partials_edge<Dx, std::decay_t<Op1>> edge1_;
  internal::ops_partials_edge<Dx, std::decay_t<Op2>> edge2_;
  internal::ops_partials_edge<Dx, std::decay_t<Op3>> edge3_;
  internal::ops_partials_edge<Dx, std::decay_t<Op4>> edge4_;
  internal::ops_partials_edge<Dx, std::decay_t<Op5>> edge5_;
  internal::ops_partials_edge<Dx, std::decay_t<Op6>> edge6_;
  internal::ops_partials_edge<Dx, std::decay_t<Op7>> edge7_;
  internal::ops_partials_edge<Dx, std::decay_t<Op8>> edge8_;
  using T_return_type = fvar<Dx>;
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
   * @return the value with its derivative
   */
  T_return_type build(Dx value) {
    Dx deriv = edge1_.dx() + edge2_.dx() + edge3_.dx() + edge4_.dx()
               + edge5_.dx() + edge6_.dx() + edge7_.dx() + edge8_.dx();
    return T_return_type(value, deriv);
  }
};

}  // namespace math
}  // namespace stan
#endif
