#ifndef STAN_MATH_FWD_FUN_LOG_DETERMINANT_HPP
#define STAN_MATH_FWD_FUN_LOG_DETERMINANT_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/fun/determinant.hpp>
#include <stan/math/fwd/fun/fabs.hpp>
#include <stan/math/fwd/fun/log.hpp>

namespace stan {
namespace math {

template <typename EigMat, require_eigen_vt<is_fvar, EigMat>* = nullptr>
inline value_type_t<EigMat> log_determinant(const EigMat& m) {
  check_square("log_determinant", "m", m);

  return log(fabs(determinant(m)));
}

}  // namespace math
}  // namespace stan
#endif
