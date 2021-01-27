#ifndef STAN_MATH_FWD_META_OPERANDS_AND_PARTIALS_HPP
#define STAN_MATH_FWD_META_OPERANDS_AND_PARTIALS_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/functor/broadcast_array.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>
#include <stan/math/fwd/core/fvar.hpp>
#include <stan/math/fwd/meta.hpp>
#include <vector>

namespace stan {
namespace math {
namespace internal {

template <typename InnerType, typename T>
class ops_partials_edge<InnerType, T, require_fvar_t<T>> {
 public:
  using Op = std::decay_t<T>;
  std::decay_t<InnerType> partial_{0};
  broadcast_array<std::decay_t<InnerType>> partials_{partial_};
  template <typename S, require_stan_scalar_t<S>* = nullptr>
  void update_partial(const S& x) {
    partial_ += x;
  }
  template <typename S, require_not_stan_scalar_t<S>* = nullptr>
  void update_partial(const S& x) {
    partial_ += sum(x);
  }
  template <typename OpT, require_same_t<OpT, Op>* = nullptr>
  explicit ops_partials_edge(const OpT& op)
      : operand_(op) {}

 private:
  template <typename, typename, typename...>
  friend class stan::math::internal::operands_and_partials_impl;
  const Op& operand_;

  auto dx() { return this->partial_ * this->operand_.d_; }
};



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
 * @tparam T_return_type return type of the expression. This defaults
 *   to a template metaprogram that calculates the scalar promotion of
 *   Op1 -- Op5
 */
 template <typename ReturnType, typename... Ops>
 class operands_and_partials_impl<ReturnType, require_fvar_t<ReturnType>, Ops...> {
   auto sum_dx() {
     return static_cast<double>(0.0);
   }

   template <typename T1>
   auto sum_dx(T1& a) {
     std::cout << "\n A: " << a.dx();
     return a.dx();
   }

   template <typename T1, typename T2>
   auto sum_dx(T1& a, T2& b) {
     std::cout << "\n A: " << a.dx();
     std::cout << "\n B: " << b.dx();
     return a.dx() + b.dx();
   }

   template <typename T1, typename T2, typename T3>
   auto sum_dx(T1& a, T2& b, T3& c) {
     std::cout << "\n A: " << a.dx();
     std::cout << "\n B: " << b.dx();
     std::cout << "\n C: " << c.dx();
     return a.dx() + b.dx() + c.dx();
   }

   template <typename T1, typename T2, typename T3, typename T4>
   auto sum_dx(T1& a, T2& b, T3& c, T4& d) {
     return a.dx() + b.dx() + c.dx() + d.dx();
   }

   template <typename T1, typename T2, typename T3, typename T4, typename T5>
   auto sum_dx(T1& a, T2& b, T3& c, T4& d, T5& ff) {
     return a.dx() + b.dx() + c.dx() + d.dx() + ff.dx();
   }

 public:
  using Dx = partials_type_t<return_type_t<Ops...>>;
  std::tuple<internal::ops_partials_edge<partials_type_t<return_type_t<Ops>>, std::decay_t<Ops>>...> edges_;
  using T_return_type = fvar<Dx>;
  template <typename... Types>
  explicit operands_and_partials_impl(Types&&... ops) :
   edges_(std::make_tuple(internal::ops_partials_edge<partials_type_t<return_type_t<Ops>>, std::decay_t<Ops>>(ops)...)) {}

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
    auto deriv = apply([&](auto&... args) {
      return this->sum_dx(args...);
    }, edges_);
    return T_return_type(value, deriv);
  }


};

// Vectorized Univariate
template <typename InnerType, typename T>
class ops_partials_edge<InnerType, T, require_std_vector_vt<is_fvar, T>> {
 public:
  using Op = std::decay_t<T>;
  using Dx = std::decay_t<InnerType>;
  using partials_t = Eigen::Matrix<Dx, -1, 1>;
  partials_t partials_;                       // For univariate use-cases
  broadcast_array<partials_t> partials_vec_;  // For multivariate
  explicit ops_partials_edge(const Op& ops)
      : partials_(partials_t::Zero(ops.size()).eval()),
        partials_vec_(partials_),
        operands_(ops) {}

 private:
  template <typename, typename, typename...>
  friend class stan::math::internal::operands_and_partials_impl;
  const Op& operands_;

  Dx dx() {
    Dx derivative(0);
    for (size_t i = 0; i < this->operands_.size(); ++i) {
      derivative += this->partials_[i] * this->operands_[i].d_;
    }
    return derivative;
  }
};

template <typename Dx, int R, int C>
class ops_partials_edge<Dx, Eigen::Matrix<fvar<Dx>, R, C>> {
 public:
  using partials_t = Eigen::Matrix<Dx, R, C>;
  using Op = Eigen::Matrix<fvar<Dx>, R, C>;
  partials_t partials_;                       // For univariate use-cases
  broadcast_array<partials_t> partials_vec_;  // For multivariate
  explicit ops_partials_edge(const Op& ops)
      : partials_(partials_t::Zero(ops.rows(), ops.cols())),
        partials_vec_(partials_),
        operands_(ops) {}

 private:
  template <typename, typename, typename...>
  friend class stan::math::internal::operands_and_partials_impl;
  const Op& operands_;

  Dx dx() {
    Dx derivative(0);
    for (int i = 0; i < this->operands_.size(); ++i) {
      derivative += this->partials_(i) * this->operands_(i).d_;
    }
    return derivative;
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

 private:
  template <typename, typename, typename...>
  friend class stan::math::internal::operands_and_partials_impl;
  const Op& operands_;

  Dx dx() {
    Dx derivative(0);
    for (size_t i = 0; i < this->operands_.size(); ++i) {
      for (int j = 0; j < this->operands_[i].size(); ++j) {
        derivative += this->partials_vec_[i](j) * this->operands_[i](j).d_;
      }
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

 private:
  template <typename, typename, typename...>
  friend class stan::math::internal::operands_and_partials_impl;
  const Op& operands_;

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
