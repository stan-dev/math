#ifndef STAN_MATH_REV_FUN_DEEP_COPY_HPP
#define STAN_MATH_REV_FUN_DEEP_COPY_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <cmath>

namespace stan {
namespace math {

namespace internal {

template <typename T, require_arithmetic_t<scalar_type_t<T>>...>
const T& deep_copy(const T& arg) {
  return arg;
}

var deep_copy(const var& arg) {
  return var(new vari(arg.val(), false));
}

std::vector<var> deep_copy(const std::vector<var>& arg) {
  std::vector<var> copy(arg.size());
  for (size_t i = 0; i < arg.size(); ++i) {
    copy[i] = new vari(arg[i].val(), false);
  }
  return copy;
}

template <typename T, require_t<is_var<scalar_type_t<T>>>...>
std::vector<T> deep_copy(const std::vector<T>& arg) {
  std::vector<T> copy(arg.size());
  for (size_t i = 0; i < arg.size(); ++i) {
    copy[i] = deep_copy(arg[i]);
  }
  return copy;
}

template <int RowType, int ColType>
Eigen::Matrix<var, RowType, ColType> deep_copy(const Eigen::Matrix<var, RowType, ColType>& arg) {
  Eigen::Matrix<var, RowType, ColType> copy(arg.rows(), arg.cols());
  for (size_t i = 0; i < arg.size(); ++i) {
    copy(i) = new vari(arg(i).val(), false);
  }
  return copy;
}
}  // namespace internal

}  // namespace math
}  // namespace stan
#endif
