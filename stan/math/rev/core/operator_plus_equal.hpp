#ifndef STAN_MATH_REV_CORE_OPERATOR_PLUS_EQUAL_HPP
#define STAN_MATH_REV_CORE_OPERATOR_PLUS_EQUAL_HPP

#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/core/operator_addition.hpp>
#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

template <typename T>
inline var_value<T>& var_value<T>::operator+=(const var_value<T>& b) {
  auto* old_vi = this->vi_;
  vi_ = new vari(this->val() + b.val(), false);
  if (unlikely(is_any_nan(old_vi->val_, b.val()))) {
    reverse_pass_callback([old_vi, b]() mutable {
      old_vi->adj_ = NOT_A_NUMBER;
      b.adj() = NOT_A_NUMBER;
    });
  } else {
    reverse_pass_callback([new_vi = this->vi_, old_vi, b]() mutable {
      old_vi->adj_ += new_vi->adj_;
      b.adj() += new_vi->adj_;
    });
  }
  return *this;
}

template <typename T>
inline var_value<T>& var_value<T>::operator+=(T b) {
  if (b == 0.0) {
    return *this;
  }
  auto* old_vi = this->vi_;
  vi_ = new vari(this->val() + b, false);
  if (unlikely(is_any_nan(old_vi->val_, b))) {
    reverse_pass_callback(
        [old_vi, b]() mutable { old_vi->adj_ = NOT_A_NUMBER; });
  } else {
    reverse_pass_callback([new_vi = this->vi_, old_vi, b]() mutable {
      old_vi->adj_ += new_vi->adj_;
    });
  }
  return *this;
}

}  // namespace math
}  // namespace stan
#endif
