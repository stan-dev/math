#ifndef STAN_MATH_PRIM_MAT_ERR_CHECK_FINITE_HPP
#define STAN_MATH_PRIM_MAT_ERR_CHECK_FINITE_HPP

#include <Eigen/Dense>
#include <stan/math/prim/scal/err/check_finite.hpp>

namespace stan {
namespace math {
namespace {
template <typename T, int R, int C>
struct finite<Eigen::Matrix<T, R, C>, true> {
  static void check(const char* function, const char* name,
                    const Eigen::Matrix<T, R, C>& y) {
    if (!y.allFinite())
      domain_error(function, name, y, "is ", ", but must be finite!");
  }
};
}
}
}

#endif
