#ifndef STAN_MATH_FWD_CORE_EIGEN_PRODUCT_KERNEL_HPP
#define STAN_MATH_FWD_CORE_EIGEN_PRODUCT_KERNEL_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/fwd/core/fvar.hpp>
#include <stan/math/fwd/core/operator_addition.hpp>
#include <stan/math/fwd/core/operator_multiplication.hpp>
#include <stan/math/fwd/core/std_complex.hpp>
#include <complex>
#include <optional>
#include <type_traits>

/**
 * @file
 * Forward-mode kernels for Eigen products whose operands have `fvar<T>` or
 * `std::complex<fvar<T>>` scalars. A product of `fvar` matrices is split
 * into products of their values and tangents, which are Eigen products of
 * the inner scalar `T`. For `T = var` those go through the reverse-mode
 * kernels, so an `fvar<var>` matrix product puts two callbacks on the
 * reverse-mode stack instead of several varis per scalar operation; for
 * `T = double` they are ordinary double products.
 *
 * The kernels must never form a coefficient-wise or matrix product with
 * `fvar<T>` scalars themselves, since that would dispatch back into these
 * kernels. Lazy results use `fwd_lazy_product_op` instead of Eigen's
 * `scalar_product_op`.
 */

namespace stan {
namespace math {
namespace internal {

/**
 * The inner scalar type `T` of `fvar<T>` and `std::complex<fvar<T>>`, and
 * `void` for other scalars.
 */
template <typename S>
struct fwd_product_inner {
  using type = void;
};
template <typename T>
struct fwd_product_inner<fvar<T>> {
  using type = T;
};
template <typename T>
struct fwd_product_inner<std::complex<fvar<T>>> {
  using type = T;
};
template <typename S>
using fwd_product_inner_t = typename fwd_product_inner<S>::type;

/**
 * True if the innermost scalar of `fvar<...<T>...>` is a reverse-mode
 * autodiff type, so products at that level record reverse-mode nodes.
 */
template <typename T>
struct fwd_product_uses_tape : is_var<T> {};
template <typename T>
struct fwd_product_uses_tape<fvar<T>> : fwd_product_uses_tape<T> {};
template <typename T>
struct fwd_product_uses_tape<std::complex<T>> : fwd_product_uses_tape<T> {};

/**
 * The value type used for an operand of a product with inner scalar `T`:
 * `T` itself for `fvar<T>`, `std::complex<T>` for `std::complex<fvar<T>>`,
 * and `double` or `std::complex<double>` for arithmetic operands.
 */
template <typename S, typename = void>
struct fwd_product_value {
  using type = double;
};
template <typename S>
struct fwd_product_value<std::complex<S>,
                         std::enable_if_t<std::is_arithmetic<S>::value>> {
  using type = std::complex<double>;
};
template <typename T>
struct fwd_product_value<fvar<T>> {
  using type = T;
};
template <typename T>
struct fwd_product_value<std::complex<fvar<T>>> {
  using type = std::complex<T>;
};

template <typename T>
inline const T& fwd_product_val(const fvar<T>& x) {
  return x.val_;
}
template <typename T>
inline const T& fwd_product_tangent(const fvar<T>& x) {
  return x.d_;
}
template <typename T>
inline std::complex<T> fwd_product_val(const std::complex<fvar<T>>& z) {
  return {z.real().val_, z.imag().val_};
}
template <typename T>
inline std::complex<T> fwd_product_tangent(const std::complex<fvar<T>>& z) {
  return {z.real().d_, z.imag().d_};
}

template <typename T>
inline const T& fwd_product_real(const T& x) {
  return x;
}
template <typename T>
inline T fwd_product_real(const std::complex<T>& z) {
  return z.real();
}
template <typename T>
inline double fwd_product_imag(const T&) {
  return 0.0;
}
template <typename T>
inline T fwd_product_imag(const std::complex<T>& z) {
  return z.imag();
}

/**
 * Build a result coefficient from its value and tangent, which may have a
 * lower scalar type (e.g. `double` for products of arithmetic operands
 * scaled by an `fvar`).
 */
template <typename T, typename V, typename D>
inline void fwd_product_assign(fvar<T>& res, const V& val, const D& d) {
  res = fvar<T>(T(val), T(d));
}
template <typename T, typename V, typename D>
inline void fwd_product_assign(std::complex<fvar<T>>& res, const V& val,
                               const D& d) {
  res = std::complex<fvar<T>>(
      fvar<T>(T(fwd_product_real(val)), T(fwd_product_real(d))),
      fvar<T>(T(fwd_product_imag(val)), T(fwd_product_imag(d))));
}

template <typename X, typename Scalar>
using fwd_operand_plain_t
    = Eigen::Matrix<Scalar, X::RowsAtCompileTime, X::ColsAtCompileTime,
                    (X::RowsAtCompileTime == 1 && X::ColsAtCompileTime != 1)
                        ? Eigen::RowMajor
                        : Eigen::ColMajor,
                    X::MaxRowsAtCompileTime, X::MaxColsAtCompileTime>;

/**
 * Operand of a forward-mode product kernel: its values and, for autodiff
 * scalars, its tangents, each a plain matrix of the inner scalar type.
 *
 * @tparam X operand expression type
 */
template <typename X,
          typename Scalar = typename Eigen::internal::traits<X>::Scalar,
          bool HasTangent = !std::is_void<fwd_product_inner_t<Scalar>>::value>
struct fwd_product_operand {
  static constexpr bool has_tangent = false;
  using value_t = typename fwd_product_value<Scalar>::type;
  fwd_operand_plain_t<X, value_t> val_;
  explicit fwd_product_operand(const X& x) : val_(x.template cast<value_t>()) {}
};

template <typename X, typename Scalar>
struct fwd_product_operand<X, Scalar, true> {
  static constexpr bool has_tangent = true;
  using value_t = typename fwd_product_value<Scalar>::type;
  fwd_operand_plain_t<X, value_t> val_;
  fwd_operand_plain_t<X, value_t> d_;
  explicit fwd_product_operand(const X& x) {
    const fwd_operand_plain_t<X, Scalar> x_eval = x;
    val_ = x_eval.unaryExpr(
        [](const Scalar& s) -> value_t { return fwd_product_val(s); });
    d_ = x_eval.unaryExpr(
        [](const Scalar& s) -> value_t { return fwd_product_tangent(s); });
  }
};

/**
 * Scalar factors folded into a forward-mode product kernel, multiplied
 * together as a scalar of the result type.
 *
 * @tparam ResScalar scalar type of the product result
 */
template <typename ResScalar>
struct fwd_product_scales {
  // empty until the first factor: constructing an `fvar<var>`, even a zero
  // or one, puts varis on the reverse-mode stack
  std::optional<ResScalar> value_;
  bool any() const { return value_.has_value(); }
  template <typename S>
  void push(const S& s) {
    if (value_) {
      *value_ *= s;
    } else {
      value_.emplace(s);
    }
  }
};

/**
 * Write the value and tangent matrices of a kernel result into its plain
 * result type.
 */
template <typename Result, typename Val, typename Tangent>
inline Result fwd_product_make_result(const Val& val, const Tangent& d) {
  Result res(val.rows(), val.cols());
  for (Eigen::Index j = 0; j < val.cols(); ++j) {
    for (Eigen::Index i = 0; i < val.rows(); ++i) {
      fwd_product_assign(res.coeffRef(i, j), val.coeff(i, j), d.coeff(i, j));
    }
  }
  return res;
}

/**
 * The value and tangent of `sigma * M` from those of `M`, where `sigma` is
 * the product of the folded scalar factors.
 */
template <typename Result, typename Scales, typename MVal, typename MTangent>
inline Result fwd_product_scaled_result(const Scales& scales, const MVal& m_val,
                                        const MTangent& m_d) {
  if (!scales.any()) {
    return fwd_product_make_result<Result>(m_val, m_d);
  }
  const auto s_val = fwd_product_val(*scales.value_);
  const auto s_d = fwd_product_tangent(*scales.value_);
  const auto val = (s_val * m_val).eval();
  const auto d = (s_d * m_val + s_val * m_d).eval();
  return fwd_product_make_result<Result>(val, d);
}

/**
 * `sigma * X * Y` for matrix operands. The tangent is `X.val * Y.d + X.d *
 * Y.val`, computed as one product of the concatenations `[X.val, X.d]` and
 * `[Y.d; Y.val]` when both operands have tangents.
 */
template <typename Result, typename Scales, typename XOp, typename YOp>
inline Result fwd_matmul(const Scales& scales, const XOp& x, const YOp& y) {
  const auto m_val = (x.val_ * y.val_).eval();
  if constexpr (XOp::has_tangent && YOp::has_tangent) {
    using x_value_t = typename XOp::value_t;
    using y_value_t = typename YOp::value_t;
    Eigen::Matrix<x_value_t, Eigen::Dynamic, Eigen::Dynamic> x_cat(
        x.val_.rows(), 2 * x.val_.cols());
    x_cat << x.val_, x.d_;
    Eigen::Matrix<y_value_t, Eigen::Dynamic, Eigen::Dynamic> y_cat(
        2 * y.val_.rows(), y.val_.cols());
    y_cat << y.d_, y.val_;
    const auto m_d = (x_cat * y_cat).eval();
    return fwd_product_scaled_result<Result>(scales, m_val, m_d);
  } else if constexpr (XOp::has_tangent) {
    const auto m_d = (x.d_ * y.val_).eval();
    return fwd_product_scaled_result<Result>(scales, m_val, m_d);
  } else if constexpr (YOp::has_tangent) {
    const auto m_d = (x.val_ * y.d_).eval();
    return fwd_product_scaled_result<Result>(scales, m_val, m_d);
  } else {
    // arithmetic operands: the tangent comes from the scalar factors only
    return fwd_product_scaled_result<Result>(
        scales, m_val,
        std::decay_t<decltype(m_val)>::Zero(m_val.rows(), m_val.cols()));
  }
}

/**
 * `sigma * (X .* Y)` for operands of equal size.
 */
template <typename Result, typename Scales, typename XOp, typename YOp>
inline Result fwd_cwise(const Scales& scales, const XOp& x, const YOp& y) {
  const auto m_val = (x.val_.array() * y.val_.array()).matrix().eval();
  if constexpr (XOp::has_tangent && YOp::has_tangent) {
    const auto m_d
        = (x.val_.array() * y.d_.array() + x.d_.array() * y.val_.array())
              .matrix()
              .eval();
    return fwd_product_scaled_result<Result>(scales, m_val, m_d);
  } else if constexpr (XOp::has_tangent) {
    const auto m_d = (x.d_.array() * y.val_.array()).matrix().eval();
    return fwd_product_scaled_result<Result>(scales, m_val, m_d);
  } else if constexpr (YOp::has_tangent) {
    const auto m_d = (x.val_.array() * y.d_.array()).matrix().eval();
    return fwd_product_scaled_result<Result>(scales, m_val, m_d);
  } else {
    // arithmetic operands: the tangent comes from the scalar factors only
    return fwd_product_scaled_result<Result>(
        scales, m_val,
        std::decay_t<decltype(m_val)>::Zero(m_val.rows(), m_val.cols()));
  }
}

/**
 * `sigma * X` for a matrix operand.
 */
template <typename Result, typename Scales, typename XOp>
inline Result fwd_scale(const Scales& scales, const XOp& x) {
  if constexpr (XOp::has_tangent) {
    return fwd_product_scaled_result<Result>(scales, x.val_, x.d_);
  } else {
    return fwd_product_scaled_result<Result>(
        scales, x.val_,
        std::decay_t<decltype(x.val_)>::Zero(x.val_.rows(), x.val_.cols())
            .eval());
  }
}

/**
 * Matrix product `lhs * rhs` of Eigen expressions with `fvar` scalars.
 *
 * @tparam Result plain Eigen type of the product
 */
template <typename Result, typename Lhs, typename Rhs>
inline Result fwd_matrix_product(const Lhs& lhs, const Rhs& rhs) {
  fwd_product_scales<typename Result::Scalar> scales;
  const auto& lhs_core = peel_scales(scales, lhs);
  const auto& rhs_core = peel_scales(scales, rhs);
  using lhs_core_t = std::decay_t<decltype(lhs_core)>;
  using rhs_core_t = std::decay_t<decltype(rhs_core)>;
  return fwd_matmul<Result>(scales, fwd_product_operand<lhs_core_t>(lhs_core),
                            fwd_product_operand<rhs_core_t>(rhs_core));
}

/**
 * Coefficient-wise product functor for lazy forward-mode products. It is
 * not Eigen's `scalar_product_op`, so expressions using it are evaluated by
 * Eigen's generic evaluator instead of Stan's product hooks.
 */
template <typename L, typename R, typename Result>
struct fwd_lazy_product_op {
  using result_type = Result;
  Result operator()(const L& a, const R& b) const { return a * b; }
};

/**
 * Coefficient-wise product expression `xpr` (including a matrix times a
 * scalar) with `fvar` scalars.
 *
 * When no reverse-mode tape is involved (e.g. `fvar<double>`), splitting
 * into value and tangent products would only add temporaries, so the
 * product stays a lazy per-coefficient expression. Otherwise the value and
 * tangent products are evaluated as products of the inner scalar, and a
 * scaled matrix product `(A * B) * s` becomes one matrix-product kernel.
 *
 * @tparam Result plain Eigen type of the expression
 */
template <typename Result, typename Xpr>
inline auto fwd_cwise_product(const Xpr& xpr) {
  using res_scalar = typename Result::Scalar;
  if constexpr (!fwd_product_uses_tape<res_scalar>::value) {
    using lhs_scalar = typename Eigen::internal::traits<
        std::decay_t<decltype(xpr.lhs())>>::Scalar;
    using rhs_scalar = typename Eigen::internal::traits<
        std::decay_t<decltype(xpr.rhs())>>::Scalar;
    return xpr.lhs().binaryExpr(
        xpr.rhs(), fwd_lazy_product_op<lhs_scalar, rhs_scalar, res_scalar>());
  } else {
    fwd_product_scales<res_scalar> scales;
    const auto& core = peel_scales(scales, xpr);
    using core_t = std::decay_t<decltype(core)>;
    if constexpr (std::is_same<core_t, Xpr>::value) {
      const auto& lhs_core = peel_scales(scales, xpr.lhs());
      const auto& rhs_core = peel_scales(scales, xpr.rhs());
      using lhs_core_t = std::decay_t<decltype(lhs_core)>;
      using rhs_core_t = std::decay_t<decltype(rhs_core)>;
      return fwd_cwise<Result>(scales,
                               fwd_product_operand<lhs_core_t>(lhs_core),
                               fwd_product_operand<rhs_core_t>(rhs_core));
    } else if constexpr (is_dense_product<core_t>::value) {
      const auto& lhs_core = peel_scales(scales, core.lhs());
      const auto& rhs_core = peel_scales(scales, core.rhs());
      using lhs_core_t = std::decay_t<decltype(lhs_core)>;
      using rhs_core_t = std::decay_t<decltype(rhs_core)>;
      return fwd_matmul<Result>(scales,
                                fwd_product_operand<lhs_core_t>(lhs_core),
                                fwd_product_operand<rhs_core_t>(rhs_core));
    } else {
      return fwd_scale<Result>(scales, fwd_product_operand<core_t>(core));
    }
  }
}

}  // namespace internal
}  // namespace math
}  // namespace stan

#endif
