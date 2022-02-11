#ifndef STAN_MATH_REV_FUNCTOR_CVODES_UTILS_HPP
#define STAN_MATH_REV_FUNCTOR_CVODES_UTILS_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <cvodes/cvodes.h>
#include <sstream>
#include <stdexcept>

namespace stan {
namespace math {

inline void cvodes_set_options(void* cvodes_mem,
                               // NOLINTNEXTLINE(runtime/int)
                               long int max_num_steps) {
  // Initialize solver parameters
  CHECK_CVODES_CALL(CVodeSetMaxNumSteps(cvodes_mem, max_num_steps));

  double init_step = 0;
  CHECK_CVODES_CALL(CVodeSetInitStep(cvodes_mem, init_step));

  long int max_err_test_fails = 20;  // NOLINT(runtime/int)
  CHECK_CVODES_CALL(CVodeSetMaxErrTestFails(cvodes_mem, max_err_test_fails));

  long int max_conv_fails = 50;  // NOLINT(runtime/int)
  CHECK_CVODES_CALL(CVodeSetMaxConvFails(cvodes_mem, max_conv_fails));
}

}  // namespace math
}  // namespace stan
#endif
