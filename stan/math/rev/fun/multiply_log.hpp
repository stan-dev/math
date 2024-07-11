#ifndef STAN_MATH_REV_FUN_MULTIPLY_LOG_HPP
#define STAN_MATH_REV_FUN_MULTIPLY_LOG_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/log.hpp>
#include <stan/math/rev/fun/elt_multiply.hpp>
#include <stan/math/rev/fun/multiply.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/multiply_log.hpp>
#include <stan/math/prim/fun/is_any_nan.hpp>
#include <cmath>

namespace stan {
namespace math {

namespace internal {
class multiply_log_vv_vari : public op_vv_vari {
 public:
  multiply_log_vv_vari(vari* avi, vari* bvi)
      : op_vv_vari(multiply_log(avi->val_, bvi->val_), avi, bvi) {}
  void chain() {
    using std::log;
    avi_->adj_ += adj_ * log(bvi_->val_);
    bvi_->adj_ += adj_ * avi_->val_ / bvi_->val_;
  }
};
class multiply_log_vd_vari : public op_vd_vari {
 public:
  multiply_log_vd_vari(vari* avi, double b)
      : op_vd_vari(multiply_log(avi->val_, b), avi, b) {}
  void chain() {
    using std::log;
    avi_->adj_ += adj_ * log(bd_);
  }
};
class multiply_log_dv_vari : public op_dv_vari {
 public:
  multiply_log_dv_vari(double a, vari* bvi)
      : op_dv_vari(multiply_log(a, bvi->val_), a, bvi) {}
  void chain() { bvi_->adj_ += adj_ * ad_ / bvi_->val_; }
};
}  // namespace internal

/**
 * Return the value of a*log(b).
 *
 * When both a and b are 0, the value returned is 0.
 * The partial derivative with respect to a is log(b).
 * The partial derivative with respect to b is a/b.
 *
 * @param a First variable.
 * @param b Second variable.
 * @return Value of a*log(b)
 */
inline var multiply_log(const var& a, const var& b) {
  return var(new internal::multiply_log_vv_vari(a.vi_, b.vi_));
}
/**
 * Return the value of a*log(b).
 *
 * When both a and b are 0, the value returned is 0.
 * The partial derivative with respect to a is log(b).
 *
 * @param a First variable.
 * @param b Second scalar.
 * @return Value of a*log(b)
 */
inline var multiply_log(const var& a, double b) {
  return var(new internal::multiply_log_vd_vari(a.vi_, b));
}
/**
 * Return the value of a*log(b).
 *
 * When both a and b are 0, the value returned is 0.
 * The partial derivative with respect to b is a/b.
 *
 * @param a First scalar.
 * @param b Second variable.
 * @return Value of a*log(b)
 */
inline var multiply_log(double a, const var& b) {
  if (a == 1.0) {
    return log(b);
  }
  return var(new internal::multiply_log_dv_vari(a, b.vi_));
}

/**
 * Return the elementwise product `a * log(b)`.
 *
 * Both `T1` and `T2` are matrices, and one of `T1` or `T2` must be a
 * `var_value`
 *
 * @tparam T1 Type of first argument
 * @tparam T2 Type of second argument
 * @param a First argument
 * @param b Second argument
 * @return elementwise product of `a` and `log(b)`
 */
template <typename T1, typename T2, require_all_matrix_t<T1, T2>* = nullptr,
          require_any_var_matrix_t<T1, T2>* = nullptr>
inline auto multiply_log(const T1& a, const T2& b) {
  check_matching_dims("multiply_log", "a", a, "b", b);
  arena_t<T1> arena_a = a;
  arena_t<T2> arena_b = b;
  if constexpr (!is_constant_v<T1> && !is_constant_v<T2>) {

    return make_callback_var(
        multiply_log(arena_a.val(), arena_b.val()),
        [arena_a, arena_b](const auto& res) mutable {
          arena_a.adj().array()
              += res.adj().array() * arena_b.val().array().log();
          arena_b.adj().array() += res.adj().array() * arena_a.val().array()
                                   / arena_b.val().array();
        });
  } else if constexpr (!is_constant_v<T1>) {
    return make_callback_var(multiply_log(arena_a.val(), arena_b),
                             [arena_a, arena_b](const auto& res) mutable {
                               arena_a.adj().array()
                                   += res.adj().array()
                                      * arena_b.val().array().log();
                             });
  } else {
    return make_callback_var(multiply_log(arena_a, arena_b.val()),
                             [arena_a, arena_b](const auto& res) mutable {
                               arena_b.adj().array() += res.adj().array()
                                                        * arena_a.val().array()
                                                        / arena_b.val().array();
                             });
  }
}

/**
 * Return the product `a * log(b)`.
 *
 * @tparam T1 Type of matrix argument
 * @tparam T2 Type of scalar argument
 * @param a Matrix argument
 * @param b Scalar argument
 * @return Product of `a` and `log(b)`
 */
template <typename T1, typename T2, require_var_matrix_t<T1>* = nullptr,
          require_stan_scalar_t<T2>* = nullptr>
inline auto multiply_log(const T1& a, const T2& b) {
  using std::log;

  arena_t<T1> arena_a = a;
  auto arena_b = b;
  if constexpr (!is_constant_v<T1> && !is_constant_v<T2>) {
    return make_callback_var(
        multiply_log(arena_a.val(), arena_b.val()),
        [arena_a, arena_b](const auto& res) mutable {
          arena_a.adj().array() += res.adj().array() * log(arena_b.val());
          arena_b.adj() += (res.adj().array() * arena_a.val().array()).sum()
                           / arena_b.val();
        });
  } else if constexpr (!is_constant_v<T1>) {
    return make_callback_var(multiply_log(arena_a.val(), b),
                             [arena_a, b](const auto& res) mutable {
                               arena_a.adj().array()
                                   += res.adj().array() * log(b);
                             });
  } else {
    return make_callback_var(
        multiply_log(arena_a, arena_b.val()),
        [arena_a, arena_b](const auto& res) mutable {
          arena_b.adj()
              += (res.adj().array() * arena_a.array()).sum() / arena_b.val();
        });
  }
}

/**
 * Return the product `a * log(b)`.
 *
 * @tparam T1 Type of scalar argument
 * @tparam T2 Type of matrix argument
 * @param a Scalar argument
 * @param b Matrix argument
 * @return Product of `a` and `log(b)`
 */
template <typename T1, typename T2, require_stan_scalar_t<T1>* = nullptr,
          require_var_matrix_t<T2>* = nullptr>
inline auto multiply_log(const T1& a, const T2& b) {
  auto arena_a = a;
  arena_t<T2> arena_b = b;
  if constexpr (!is_constant_v<T1> && !is_constant_v<T2>) {
    return make_callback_var(
        multiply_log(arena_a.val(), arena_b.val()),
        [arena_a, arena_b](const auto& res) mutable {
          arena_a.adj()
              += (res.adj().array() * arena_b.val().array().log()).sum();
          arena_b.adj().array()
              += arena_a.val() * res.adj().array() / arena_b.val().array();
        });
  } else if constexpr (!is_constant_v<T1>) {
    return make_callback_var(
        multiply_log(arena_a.val(), arena_b),
        [arena_a, arena_b](const auto& res) mutable {
          arena_a.adj()
              += (res.adj().array() * arena_b.val().array().log()).sum();
        });
  } else {
    return make_callback_var(multiply_log(value_of(a), arena_b.val()),
                             [a, arena_b](const auto& res) mutable {
                               arena_b.adj().array() += value_of(a)
                                                        * res.adj().array()
                                                        / arena_b.val().array();
                             });
  }
}

}  // namespace math
}  // namespace stan
#endif
