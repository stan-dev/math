#ifndef STAN_MATH_REV_MAT_FUNCTOR_INTEGRATE_ODE_ADAMS_HPP
#define STAN_MATH_REV_MAT_FUNCTOR_INTEGRATE_ODE_ADAMS_HPP

#include <stan/math/prim/arr/fun/value_of.hpp>
#include <stan/math/prim/scal/err/check_less.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/arr/err/check_nonzero_size.hpp>
#include <stan/math/prim/arr/err/check_ordered.hpp>
#include <stan/math/rev/scal/meta/is_var.hpp>
#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/rev/mat/functor/cvodes_utils.hpp>
#include <stan/math/rev/mat/functor/cvodes_ode_data.hpp>
#include <stan/math/rev/mat/functor/cvodes_integrator.hpp>
#include <stan/math/rev/arr/fun/decouple_ode_states.hpp>
#include <cvodes/cvodes.h>
#include <cvodes/cvodes_band.h>
#include <cvodes/cvodes_dense.h>
#include <nvector/nvector_serial.h>
#include <algorithm>
#include <ostream>
#include <vector>

namespace stan {
namespace math {

template <typename F, typename T_initial, typename T_param>
std::vector<std::vector<typename stan::return_type<T_initial, T_param>::type> >
integrate_ode_adams(const F& f, const std::vector<T_initial>& y0, double t0,
                  const std::vector<double>& ts,
                  const std::vector<T_param>& theta,
                  const std::vector<double>& x, const std::vector<int>& x_int,
                  std::ostream* msgs = 0, double relative_tolerance = 1e-10,
                  double absolute_tolerance = 1e-10,
                  long int max_num_steps = 1e8) {  // NOLINT(runtime/int)
  stan::math::cvodes_integrator<CV_ADAMS> integrator;
  return integrator.integrate(f, y0, t0,
                              ts,
                              theta,
                              x, x_int,
                              msgs, relative_tolerance,
                              absolute_tolerance,
                              max_num_steps);
}

}  // namespace math
}  // namespace stan
#endif
