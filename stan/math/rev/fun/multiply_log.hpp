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
    if (unlikely(is_any_nan(avi_->val_, bvi_->val_))) {
      avi_->adj_ = NOT_A_NUMBER;
      bvi_->adj_ = NOT_A_NUMBER;
    } else {
      avi_->adj_ += adj_ * log(bvi_->val_);
      if (bvi_->val_ == 0.0 && avi_->val_ == 0) {
        bvi_->adj_ += adj_ * INFTY;
      } else {
        bvi_->adj_ += adj_ * avi_->val_ / bvi_->val_;
      }
    }
  }
};
class multiply_log_vd_vari : public op_vd_vari {
 public:
  multiply_log_vd_vari(vari* avi, double b)
      : op_vd_vari(multiply_log(avi->val_, b), avi, b) {}
  void chain() {
    using std::log;
    if (unlikely(is_any_nan(avi_->val_, bd_))) {
      avi_->adj_ = NOT_A_NUMBER;
    } else {
      avi_->adj_ += adj_ * log(bd_);
    }
  }
};
class multiply_log_dv_vari : public op_dv_vari {
 public:
  multiply_log_dv_vari(double a, vari* bvi)
      : op_dv_vari(multiply_log(a, bvi->val_), a, bvi) {}
  void chain() {
    if (bvi_->val_ == 0.0 && ad_ == 0.0) {
      bvi_->adj_ += adj_ * INFTY;
    } else {
      bvi_->adj_ += adj_ * ad_ / bvi_->val_;
    }
  }
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
  if (!is_constant<T1>::value && !is_constant<T2>::value) {
    arena_t<promote_scalar_t<var, T1>> arena_a = a;
    arena_t<promote_scalar_t<var, T2>> arena_b = b;

    return make_callback_var(multiply_log(arena_a.val(), arena_b.val()),
			     [arena_a, arena_b](const auto& res) mutable {
			       for(Eigen::Index j = 0; j < res.adj().cols(); ++j) {
				 for(Eigen::Index i = 0; i < res.adj().rows(); ++i) {
				   if (unlikely(is_any_nan(arena_a.val().coeff(i, j), arena_b.val().coeff(i, j)))) {
				     arena_a.adj().coeffRef(i, j) = NOT_A_NUMBER;
				     arena_b.adj().coeffRef(i, j) = NOT_A_NUMBER;
				   } else {
				     arena_a.adj().coeffRef(i, j) += res.adj().coeff(i, j) * log(arena_b.val().coeff(i, j));
				     if (arena_b.val().coeff(i, j) == 0.0 && arena_a.val().coeff(i, j) == 0) {
				       arena_b.adj().coeffRef(i, j) += res.adj().coeff(i, j) * INFTY;
				     } else {
				       arena_b.adj().coeffRef(i, j) += res.adj().coeff(i, j) * arena_a.val().coeff(i, j) / arena_b.val().coeff(i, j);
				     }
				   }
				 }
			       }
			     });
  } else if (!is_constant<T1>::value) {
    arena_t<promote_scalar_t<var, T1>> arena_a = a;
    auto arena_b = to_arena(value_of(b));

    return make_callback_var(multiply_log(arena_a.val(), arena_b),
                             [arena_a, arena_b](const auto& res) mutable {
			       for(Eigen::Index j = 0; j < res.adj().cols(); ++j) {
				 for(Eigen::Index i = 0; i < res.adj().rows(); ++i) {
				   if (unlikely(is_any_nan(arena_a.val().coeff(i, j), arena_b.coeff(i, j)))) {
				     arena_a.adj().coeffRef(i, j) = NOT_A_NUMBER;
				   } else {
				     arena_a.adj().coeffRef(i, j) += res.adj().coeff(i, j) * log(arena_b.coeff(i, j));
				   }
				 }
			       }
                             });
  } else {
    auto arena_a = to_arena(value_of(a));
    arena_t<promote_scalar_t<var, T2>> arena_b = b;

    return make_callback_var(multiply_log(arena_a, arena_b.val()),
			     [arena_a, arena_b](const auto& res) mutable {
			       for(Eigen::Index j = 0; j < res.adj().cols(); ++j) {
				 for(Eigen::Index i = 0; i < res.adj().rows(); ++i) {
				   if (unlikely(is_any_nan(arena_a.val().coeff(i, j), arena_b.val().coeff(i, j)))) {
				     arena_b.adj().coeffRef(i, j) = NOT_A_NUMBER;
				   } else {
				     if (arena_b.val().coeff(i, j) == 0.0 && arena_a.val().coeff(i, j) == 0) {
				       arena_b.adj().coeffRef(i, j) += res.adj().coeff(i, j) * INFTY;
				     } else {
				       arena_b.adj().coeffRef(i, j) += res.adj().coeff(i, j) * arena_a.coeff(i, j) / arena_b.val().coeff(i, j);
				     }
				   }
				 }
			       }
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

  if (!is_constant<T1>::value && !is_constant<T2>::value) {
    arena_t<promote_scalar_t<var, T1>> arena_a = a;
    var arena_b = b;

    return make_callback_var(multiply_log(arena_a.val(), arena_b.val()),
        [arena_a, arena_b](const auto& res) mutable {
			       for(Eigen::Index j = 0; j < res.adj().cols(); ++j) {
				 for(Eigen::Index i = 0; i < res.adj().rows(); ++i) {
				   if (unlikely(is_any_nan(arena_a.val().coeff(i, j), arena_b.val()))) {
				     arena_a.adj().coeffRef(i, j) = NOT_A_NUMBER;
				     arena_b.adj() = NOT_A_NUMBER;
				   } else {
				     arena_a.adj().coeffRef(i, j) += res.adj().coeff(i, j) * log(arena_b.val());
				     if (arena_b.val() == 0.0 && arena_a.val().coeff(i, j) == 0) {
				       arena_b.adj() += res.adj().coeff(i, j) * INFTY;
				     } else {
				       arena_b.adj() += res.adj().coeff(i, j) * arena_a.val().coeff(i, j) / arena_b.val();
				     }
				   }
				 }
			       }
        });
  } else if (!is_constant<T1>::value) {
    arena_t<promote_scalar_t<var, T1>> arena_a = a;

    return make_callback_var(multiply_log(arena_a.val(), value_of(b)),
                             [arena_a, b](const auto& res) mutable {
			       for(Eigen::Index j = 0; j < res.adj().cols(); ++j) {
				 for(Eigen::Index i = 0; i < res.adj().rows(); ++i) {
				   if (unlikely(is_any_nan(arena_a.val().coeff(i, j), value_of(b)))) {
				     arena_a.adj().coeffRef(i, j) = NOT_A_NUMBER;
				   } else {
				     arena_a.adj().coeffRef(i, j) += res.adj().coeff(i, j) * log(value_of(b));
				   }
				 }
			       }
                             });
  } else {
    arena_t<promote_scalar_t<double, T1>> arena_a = value_of(a);
    var arena_b = b;

    return make_callback_var(multiply_log(arena_a, arena_b.val()),
        [arena_a, arena_b](const auto& res) mutable {
			       for(Eigen::Index j = 0; j < res.adj().cols(); ++j) {
				 for(Eigen::Index i = 0; i < res.adj().rows(); ++i) {
				   if (unlikely(is_any_nan(arena_a.val().coeff(i, j), arena_b.val()))) {
				     arena_b.adj() = NOT_A_NUMBER;
				   } else {
				     if (arena_b.val() == 0.0 && arena_a.val().coeff(i, j) == 0) {
				       arena_b.adj() += res.adj().coeff(i, j) * INFTY;
				     } else {
				       arena_b.adj() += res.adj().coeff(i, j) * arena_a.val().coeff(i, j) / arena_b.val();
				     }
				   }
				 }
			       }
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
  if (!is_constant<T1>::value && !is_constant<T2>::value) {
    var arena_a = a;
    arena_t<promote_scalar_t<var, T2>> arena_b = b;

    return make_callback_var(multiply_log(arena_a.val(), arena_b.val()),
        [arena_a, arena_b](const auto& res) mutable {
			       for(Eigen::Index j = 0; j < res.adj().cols(); ++j) {
				 for(Eigen::Index i = 0; i < res.adj().rows(); ++i) {
				   if (unlikely(is_any_nan(arena_a.val(), arena_b.val().coeff(i, j)))) {
				     arena_a.adj() = NOT_A_NUMBER;
				     arena_b.adj().coeffRef(i, j) = NOT_A_NUMBER;
				   } else {
				     arena_a.adj() += res.adj().coeff(i, j) * log(arena_b.val().coeff(i, j));
				     if (arena_b.val().coeff(i, j) == 0.0 && arena_a.val() == 0) {
				       arena_b.adj().coeffRef(i, j) += res.adj().coeff(i, j) * INFTY;
				     } else {
				       arena_b.adj().coeffRef(i, j) += res.adj().coeff(i, j) * arena_a.val() / arena_b.val().coeff(i, j);
				     }
				   }
				 }
			       }
        });
  } else if (!is_constant<T1>::value) {
    var arena_a = a;
    auto arena_b = to_arena(value_of(b));

    return make_callback_var(multiply_log(arena_a.val(), arena_b),
                             [arena_a, arena_b](const auto& res) mutable {
			       for(Eigen::Index j = 0; j < res.adj().cols(); ++j) {
				 for(Eigen::Index i = 0; i < res.adj().rows(); ++i) {
				   if (unlikely(is_any_nan(arena_a.val(), arena_b.val().coeff(i, j)))) {
				     arena_a.adj() = NOT_A_NUMBER;
				   } else {
				     arena_a.adj() += res.adj().coeff(i, j) * log(arena_b.val().coeff(i, j));
				   }
				 }
			       }
			     });
  } else {
    arena_t<promote_scalar_t<var, T2>> arena_b = b;

    return make_callback_var(multiply_log(value_of(a), arena_b.val()),
        [a, arena_b](const auto& res) mutable {
			       for(Eigen::Index j = 0; j < res.adj().cols(); ++j) {
				 for(Eigen::Index i = 0; i < res.adj().rows(); ++i) {
				   if (unlikely(is_any_nan(value_of(a), arena_b.val().coeff(i, j)))) {
				     arena_b.adj().coeffRef(i, j) = NOT_A_NUMBER;
				   } else {
				     if (arena_b.val().coeff(i, j) == 0.0 && value_of(a) == 0) {
				       arena_b.adj().coeffRef(i, j) += res.adj().coeff(i, j) * INFTY;
				     } else {
				       arena_b.adj().coeffRef(i, j) += res.adj().coeff(i, j) * value_of(a) / arena_b.val().coeff(i, j);
				     }
				   }
				 }
			       }
        });
  }
}

}  // namespace math
}  // namespace stan
#endif
