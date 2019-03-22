#ifndef STAN_MATH_PRIM_SCAL_FUN_AS_SCALAR_HPP
#define STAN_MATH_PRIM_SCAL_FUN_AS_SCALAR_HPP

#include <Eigen\Dense>

namespace stan {
namespace math {

/**
 * Converts input to a scalar. For scalar arguments this is an identity function.
 * @param a Input value
 * @return Same value
 */
double as_scalar(double a){
  return a;
}

}
}

#endif
