#ifndef STAN_MATH_REV_CORE_EIGEN_PRODUCT_KERNEL_HPP
#define STAN_MATH_REV_CORE_EIGEN_PRODUCT_KERNEL_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/rev/core/arena_matrix.hpp>
#include <stan/math/rev/core/reverse_pass_callback.hpp>
#include <stan/math/rev/core/std_complex.hpp>
#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/core/vari.hpp>
#include <complex>
#include <type_traits>

/**
 * @file
 * Reverse-mode kernels for Eigen products whose operands have `var` or
 * `std::complex<var>` scalars. Each product puts one callback on the
 * reverse-mode stack plus one value-only (nochain) vari per result
 * coefficient, instead of one vari per scalar operation.
 *
 * All arithmetic inside the kernels is on `double` or `std::complex<double>`
 * values. The kernels must never form an Eigen product of `var` scalars
 * themselves, since that would dispatch back into these kernels.
 */

namespace stan {
namespace math {
namespace internal {

/**
 * True for `var` and `std::complex<var>`.
 */
template <typename T>
struct is_rev_product_scalar : std::false_type {};
template <>
struct is_rev_product_scalar<var> : std::true_type {};
template <>
struct is_rev_product_scalar<std::complex<var>> : std::true_type {};

/**
 * True for scalars a reverse-mode product may mix with `var`.
 */
template <typename T>
struct is_rev_product_operand_scalar
    : bool_constant<std::is_arithmetic<T>::value
                    || is_rev_product_scalar<T>::value> {};
template <typename T>
struct is_rev_product_operand_scalar<std::complex<T>>
    : bool_constant<std::is_arithmetic<T>::value
                    || std::is_same<T, var>::value> {};

/**
 * The value type used for the arithmetic of a product whose result has
 * scalar type `T`: `double`, or `std::complex<double>` for complex results.
 */
template <typename T>
struct rev_product_value {
  using type = double;
};
template <typename T>
struct rev_product_value<std::complex<T>> {
  using type = std::complex<double>;
};
template <typename T>
using rev_product_value_t = typename rev_product_value<T>::type;

template <typename T, require_arithmetic_t<T>* = nullptr>
inline double rev_product_value_of(const T& x) {
  return x;
}
template <typename T, require_arithmetic_t<T>* = nullptr>
inline std::complex<double> rev_product_value_of(const std::complex<T>& z) {
  return {static_cast<double>(z.real()), static_cast<double>(z.imag())};
}
inline double rev_product_value_of(const var& x) { return x.val(); }
inline std::complex<double> rev_product_value_of(const std::complex<var>& z) {
  return {z.real().val(), z.imag().val()};
}

inline double rev_product_conj(double x) { return x; }
inline std::complex<double> rev_product_conj(const std::complex<double>& z) {
  return std::conj(z);
}

/**
 * Add the adjoint `g` of a product result to operand `x`. For a real `x` and
 * a complex `g`, only the real part contributes.
 */
inline void rev_product_add_adjoint(const var& x, double g) {
  x.vi_->adj_ += g;
}
inline void rev_product_add_adjoint(const var& x,
                                    const std::complex<double>& g) {
  x.vi_->adj_ += g.real();
}
inline void rev_product_add_adjoint(const std::complex<var>& z,
                                    const std::complex<double>& g) {
  z.real().vi_->adj_ += g.real();
  z.imag().vi_->adj_ += g.imag();
}

/**
 * A new value-only vari holding `x`. Results of a product kernel get their
 * adjoints from the kernel's callback, so they need no `chain()`.
 */
inline var rev_product_result(double x) { return var(new vari(x, false)); }
inline std::complex<var> rev_product_result(const std::complex<double>& z) {
  return {var(new vari(z.real(), false)), var(new vari(z.imag(), false))};
}

/**
 * Operand of a reverse-mode product kernel: its values, and for autodiff
 * scalars its vars, copied to the arena once.
 *
 * @tparam Scalar scalar type of the operand expression
 */
template <typename Scalar,
          bool IsAutodiff = is_rev_product_scalar<Scalar>::value>
struct rev_product_operand {
  static constexpr bool is_autodiff = false;
  using value_t = rev_product_value_t<Scalar>;
  arena_matrix<Eigen::Matrix<value_t, Eigen::Dynamic, Eigen::Dynamic>> val_;
  template <typename X>
  explicit rev_product_operand(const X& x) : val_(x.template cast<value_t>()) {}
  const auto& val() const { return val_; }
  template <typename G>
  void add_adjoint(const G&) const {}
};

template <typename Scalar>
struct rev_product_operand<Scalar, true> {
  static constexpr bool is_autodiff = true;
  using value_t = rev_product_value_t<Scalar>;
  arena_matrix<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>> vars_;
  arena_matrix<Eigen::Matrix<value_t, Eigen::Dynamic, Eigen::Dynamic>> val_;
  template <typename X>
  explicit rev_product_operand(const X& x)
      : vars_(x), val_(vars_.unaryExpr([](const Scalar& s) {
          return rev_product_value_of(s);
        })) {}
  const auto& val() const { return val_; }
  template <typename G>
  void add_adjoint(const G& g) const {
    using g_scalar = typename Eigen::internal::traits<G>::Scalar;
    if constexpr (std::is_same<Scalar, var>::value
                  && std::is_same<g_scalar, double>::value) {
      // a shallow copy of the arena map, whose adjoint view is writable
      arena_matrix<Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>> vars
          = vars_;
      vars.adj() += g;
    } else {
      const auto& g_ref = g.eval();
      for (Eigen::Index j = 0; j < vars_.cols(); ++j) {
        for (Eigen::Index i = 0; i < vars_.rows(); ++i) {
          rev_product_add_adjoint(vars_.coeff(i, j), g_ref.coeff(i, j));
        }
      }
    }
  }
};

/**
 * Scalar factors folded into a product kernel: the product of the
 * arithmetic factors, and the autodiff factors.
 *
 * @tparam ValueT `double` or `std::complex<double>`
 */
template <typename ValueT>
struct rev_product_scales {
  // three peels (expression, left and right operand) of up to
  // max_peel_depth factors each
  static constexpr int max_autodiff = 3 * max_peel_depth;
  struct autodiff_scale {
    vari* re_;
    vari* im_;  // nullptr for a real scale
    ValueT value() const {
      if constexpr (std::is_same<ValueT, double>::value) {
        return re_->val_;
      } else {
        return {re_->val_, im_ == nullptr ? 0.0 : im_->val_};
      }
    }
  };
  ValueT constant_{1.0};
  // allocated on the arena at the first autodiff factor, so a kernel's
  // callback captures a pointer instead of storage for every factor
  autodiff_scale* autodiff_{nullptr};
  int n_autodiff_{0};

  template <typename T, require_arithmetic_t<T>* = nullptr>
  void push(const T& x) {
    constant_ *= static_cast<double>(x);
  }
  template <typename T, require_arithmetic_t<T>* = nullptr>
  void push(const std::complex<T>& x) {
    constant_ *= std::complex<double>(x.real(), x.imag());
  }
  void push(const var& x) { push_autodiff({x.vi_, nullptr}); }
  void push(const std::complex<var>& x) {
    push_autodiff({x.real().vi_, x.imag().vi_});
  }

  /** Product of every factor. */
  ValueT value() const {
    ValueT v = constant_;
    for (int i = 0; i < n_autodiff_; ++i) {
      v *= autodiff_[i].value();
    }
    return v;
  }

  /**
   * Given `t = sum(G .* conj(M))` for the unscaled product `M` and result
   * adjoints `G`, add each autodiff factor's adjoint.
   */
  void add_adjoints(const ValueT& t) const {
    for (int i = 0; i < n_autodiff_; ++i) {
      ValueT others = constant_;
      for (int k = 0; k < n_autodiff_; ++k) {
        if (k != i) {
          others *= autodiff_[k].value();
        }
      }
      const ValueT g = rev_product_conj(others) * t;
      if constexpr (std::is_same<ValueT, double>::value) {
        autodiff_[i].re_->adj_ += g;
      } else {
        autodiff_[i].re_->adj_ += g.real();
        if (autodiff_[i].im_ != nullptr) {
          autodiff_[i].im_->adj_ += g.imag();
        }
      }
    }
  }

 private:
  void push_autodiff(const autodiff_scale& s) {
    if (autodiff_ == nullptr) {
      autodiff_
          = ChainableStack::instance_->memalloc_.alloc_array<autodiff_scale>(
              max_autodiff);
    }
    autodiff_[n_autodiff_++] = s;
  }
};

template <typename Result>
using rev_result_value_t = rev_product_value_t<typename Result::Scalar>;

/**
 * Build the plain result of a product kernel from its values, with one
 * nochain vari per coefficient, and keep an arena copy for the callback.
 */
template <typename Result, typename Values>
inline Result rev_product_make_result(
    const Values& values,
    arena_matrix<Eigen::Matrix<typename Result::Scalar, Eigen::Dynamic,
                               Eigen::Dynamic>>& res_arena) {
  res_arena = arena_matrix<
      Eigen::Matrix<typename Result::Scalar, Eigen::Dynamic, Eigen::Dynamic>>(
      values.rows(), values.cols());
  Result res(values.rows(), values.cols());
  for (Eigen::Index j = 0; j < values.cols(); ++j) {
    for (Eigen::Index i = 0; i < values.rows(); ++i) {
      res_arena.coeffRef(i, j) = rev_product_result(values.coeff(i, j));
      res.coeffRef(i, j) = res_arena.coeff(i, j);
    }
  }
  return res;
}

/**
 * Adjoints of a kernel result as a lazy `double` or `std::complex<double>`
 * expression.
 */
template <typename ResScalar>
inline auto rev_product_result_adjoints(
    const arena_matrix<
        Eigen::Matrix<ResScalar, Eigen::Dynamic, Eigen::Dynamic>>& res) {
  if constexpr (std::is_same<ResScalar, var>::value) {
    return res.adj();
  } else {
    return res.unaryExpr([](const ResScalar& z) {
      return std::complex<double>(z.real().adj(), z.imag().adj());
    });
  }
}

template <typename X>
using rev_operand_t
    = rev_product_operand<typename Eigen::internal::traits<X>::Scalar>;

/**
 * `sigma * X * Y` for matrix operands `X` and `Y`, with one callback.
 *
 * @tparam Result plain Eigen type of the product
 */
template <typename Result, typename Scales, typename XOp, typename YOp>
inline Result rev_matmul(const Scales& scales, const XOp& x, const YOp& y) {
  using value_t = rev_result_value_t<Result>;
  using res_scalar = typename Result::Scalar;
  const value_t sigma = scales.value();
  arena_matrix<Eigen::Matrix<res_scalar, Eigen::Dynamic, Eigen::Dynamic>>
      res_arena;
  arena_matrix<Eigen::Matrix<value_t, Eigen::Dynamic, Eigen::Dynamic>> unscaled;
  Result res;
  if (scales.n_autodiff_ > 0) {
    unscaled = x.val() * y.val();
    res = rev_product_make_result<Result>((sigma * unscaled).eval(), res_arena);
  } else {
    res = rev_product_make_result<Result>((sigma * (x.val() * y.val())).eval(),
                                          res_arena);
  }
  reverse_pass_callback([res_arena, x, y, scales, sigma, unscaled]() {
    // evaluated once: both operand gradients are products with it
    const auto g = rev_product_result_adjoints(res_arena).eval();
    const value_t sigma_conj = rev_product_conj(sigma);
    if constexpr (XOp::is_autodiff) {
      x.add_adjoint(sigma_conj * (g * y.val().adjoint()));
    }
    if constexpr (YOp::is_autodiff) {
      y.add_adjoint(sigma_conj * (x.val().adjoint() * g));
    }
    if (scales.n_autodiff_ > 0) {
      scales.add_adjoints((g.array() * unscaled.array().conjugate()).sum());
    }
  });
  return res;
}

/**
 * `sigma * (X .* Y)` for operands of equal size, with one callback.
 */
template <typename Result, typename Scales, typename XOp, typename YOp>
inline Result rev_cwise(const Scales& scales, const XOp& x, const YOp& y) {
  using value_t = rev_result_value_t<Result>;
  using res_scalar = typename Result::Scalar;
  const value_t sigma = scales.value();
  arena_matrix<Eigen::Matrix<res_scalar, Eigen::Dynamic, Eigen::Dynamic>>
      res_arena;
  Result res = rev_product_make_result<Result>(
      (sigma * (x.val().array() * y.val().array())).matrix().eval(), res_arena);
  reverse_pass_callback([res_arena, x, y, scales, sigma]() {
    const auto g = rev_product_result_adjoints(res_arena);
    const value_t sigma_conj = rev_product_conj(sigma);
    if constexpr (XOp::is_autodiff) {
      x.add_adjoint(
          (sigma_conj * (g.array() * y.val().array().conjugate())).matrix());
    }
    if constexpr (YOp::is_autodiff) {
      y.add_adjoint(
          (sigma_conj * (g.array() * x.val().array().conjugate())).matrix());
    }
    if (scales.n_autodiff_ > 0) {
      scales.add_adjoints(
          (g.array() * (x.val().array() * y.val().array()).conjugate()).sum());
    }
  });
  return res;
}

/**
 * `sigma * X` for a matrix operand `X`, with one callback.
 */
template <typename Result, typename Scales, typename XOp>
inline Result rev_scale(const Scales& scales, const XOp& x) {
  using value_t = rev_result_value_t<Result>;
  using res_scalar = typename Result::Scalar;
  const value_t sigma = scales.value();
  arena_matrix<Eigen::Matrix<res_scalar, Eigen::Dynamic, Eigen::Dynamic>>
      res_arena;
  Result res
      = rev_product_make_result<Result>((sigma * x.val()).eval(), res_arena);
  reverse_pass_callback([res_arena, x, scales, sigma]() {
    const auto g = rev_product_result_adjoints(res_arena);
    if constexpr (XOp::is_autodiff) {
      x.add_adjoint(rev_product_conj(sigma) * g);
    }
    if (scales.n_autodiff_ > 0) {
      scales.add_adjoints((g.array() * x.val().array().conjugate()).sum());
    }
  });
  return res;
}

/**
 * Matrix product `lhs * rhs` of Eigen expressions. Scalar factors of either
 * operand are folded into the product's callback.
 *
 * @tparam Result plain Eigen type of the product
 */
template <typename Result, typename Lhs, typename Rhs>
inline Result rev_matrix_product(const Lhs& lhs, const Rhs& rhs) {
  rev_product_scales<rev_result_value_t<Result>> scales;
  const auto& lhs_core = peel_scales(scales, lhs);
  const auto& rhs_core = peel_scales(scales, rhs);
  using lhs_core_t = std::decay_t<decltype(lhs_core)>;
  using rhs_core_t = std::decay_t<decltype(rhs_core)>;
  return rev_matmul<Result>(scales, rev_operand_t<lhs_core_t>(lhs_core),
                            rev_operand_t<rhs_core_t>(rhs_core));
}

/**
 * Coefficient-wise product expression `xpr` (including a matrix times a
 * scalar). A scaled matrix product `(A * B) * s` becomes one matrix-product
 * callback.
 *
 * @tparam Result plain Eigen type of the expression
 */
template <typename Result, typename Xpr>
inline Result rev_cwise_product(const Xpr& xpr) {
  using value_t = rev_result_value_t<Result>;
  rev_product_scales<value_t> scales;
  const auto& core = peel_scales(scales, xpr);
  using core_t = std::decay_t<decltype(core)>;
  if constexpr (std::is_same<core_t, Xpr>::value) {
    // a coefficient-wise product of two non-constant operands
    const auto& lhs_core = peel_scales(scales, xpr.lhs());
    const auto& rhs_core = peel_scales(scales, xpr.rhs());
    using lhs_core_t = std::decay_t<decltype(lhs_core)>;
    using rhs_core_t = std::decay_t<decltype(rhs_core)>;
    return rev_cwise<Result>(scales, rev_operand_t<lhs_core_t>(lhs_core),
                             rev_operand_t<rhs_core_t>(rhs_core));
  } else if constexpr (is_dense_product<core_t>::value) {
    const auto& lhs_core = peel_scales(scales, core.lhs());
    const auto& rhs_core = peel_scales(scales, core.rhs());
    using lhs_core_t = std::decay_t<decltype(lhs_core)>;
    using rhs_core_t = std::decay_t<decltype(rhs_core)>;
    return rev_matmul<Result>(scales, rev_operand_t<lhs_core_t>(lhs_core),
                              rev_operand_t<rhs_core_t>(rhs_core));
  } else {
    return rev_scale<Result>(scales, rev_operand_t<core_t>(core));
  }
}

}  // namespace internal
}  // namespace math
}  // namespace stan

#endif
