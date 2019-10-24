#ifndef STAN_MATH_REV_MAT_FUNCTOR_CVODES_UTILS_HPP
#define STAN_MATH_REV_MAT_FUNCTOR_CVODES_UTILS_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/mat/err/check_flag_sundials.hpp>
#include <cvodes/cvodes.h>
#include <sstream>
#include <stdexcept>

namespace stan {
namespace math {

// no-op error handler to silence CVodes error output;  errors handled
// directly by Stan
extern "C" inline void cvodes_silent_err_handler(int error_code,
                                                 const char* module,
                                                 const char* function,
                                                 char* msg, void* eh_data) {}

inline void cvodes_set_options(void* cvodes_mem, double rel_tol, double abs_tol,
                               // NOLINTNEXTLINE(runtime/int)
                               long int max_num_steps) {
  // forward CVode errors to noop error handler
  CVodeSetErrHandlerFn(cvodes_mem, cvodes_silent_err_handler, nullptr);

  // Initialize solver parameters
  check_flag_sundials(CVodeSStolerances(cvodes_mem, rel_tol, abs_tol),
                      "CVodeSStolerances");

  check_flag_sundials(CVodeSetMaxNumSteps(cvodes_mem, max_num_steps),
                      "CVodeSetMaxNumSteps");

  double init_step = 0;
  check_flag_sundials(CVodeSetInitStep(cvodes_mem, init_step),
                      "CVodeSetInitStep");

  long int max_err_test_fails = 20;  // NOLINT(runtime/int)
  check_flag_sundials(CVodeSetMaxErrTestFails(cvodes_mem, max_err_test_fails),
                      "CVodeSetMaxErrTestFails");

  long int max_conv_fails = 50;  // NOLINT(runtime/int)
  check_flag_sundials(CVodeSetMaxConvFails(cvodes_mem, max_conv_fails),
                      "CVodeSetMaxConvFails");
}

}  // namespace math
}  // namespace stan
#endif
