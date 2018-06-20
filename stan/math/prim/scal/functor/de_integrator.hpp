/*
 * The original code on which this file is based is from John Cook,
 * which is licensed under the 2-clause BSD license, reproduced below:
 *
 * Copyright (c) 2015, John D. Cook
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 * ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 */
#ifndef STAN_MATH_PRIM_SCAL_FUNCTOR_DEINTEGRATOR_HPP
#define STAN_MATH_PRIM_SCAL_FUNCTOR_DEINTEGRATOR_HPP

#include <stan/math/rev/mat/functor/de_integrator_constants.hpp>
#include <cmath>
#include <cfloat>

namespace stan {

namespace math {

/**
 * Double Exponential Integrator
 *
 * @tparam T Type of f
 * @param f a functor with signature double (double)
 * @param a lower limit of integration, must be double type
 * @param b upper limit of integration, must be double type
 * @param tolerance target absolute error
 * @return numeric integral of function f
 */
template <typename F>
inline double de_integrator(const F& f, double a, double b, double tolerance) {
  using std::fabs;
  using std::log;

  // Apply the linear change of variables x = ct + d
  // $$\int_a^b f(x) dx = c \int_{-1}^1 f( ct + d ) dt$$
  // c = (b-a)/2, d = (a+b)/2
  double c = 0.5 * (b - a);

  // If size of integral is zero, handle this specially to avoid divide by zero
  if (c == 0.0) {
    return 0.0;
  }

  double d = 0.5 * (a + b);
  int num_function_evaluations;

  tolerance /= c;

  // Offsets to where each level's integration constants start.
  // The last element is not a beginning but an end.
  static int offsets[] = {1, 4, 7, 13, 25, 49, 97, 193};
  int num_levels = sizeof(offsets) / sizeof(int) - 1;

  double new_contribution = 0.0;
  double integral = 0.0;
  double error_estimate = DBL_MAX;
  double h = 1.0;
  double previous_delta, current_delta = DBL_MAX;

  integral = f(c * de_abcissas[0] + d) * de_weights[0];

  int i;
  for (i = offsets[0]; i != offsets[1]; ++i)
    integral += de_weights[i]
                * (f(c * de_abcissas[i] + d) + f(-c * de_abcissas[i] + d));

  int level;
  for (level = 1; level != num_levels; ++level) {
    h *= 0.5;
    new_contribution = 0.0;
    for (i = offsets[level]; i != offsets[level + 1]; ++i)
      new_contribution
          += de_weights[i]
             * (f(c * de_abcissas[i] + d) + f(-c * de_abcissas[i] + d));
    new_contribution *= h;

    // difference in consecutive integral estimates
    previous_delta = current_delta;
    current_delta = fabs(0.5 * integral - new_contribution);
    integral = 0.5 * integral + new_contribution;

    // Once convergence kicks in, error is approximately squared
    // at each step.
    // Determine whether we're in the convergent region by looking
    // at the trend in the error.
    if (level == 1)
      // previous_delta meaningless, so cannot check
      // convergence.
      continue;

    // Exact comparison with zero is harmless here. Could possibly
    // be replaced with a small positive upper limit on the size
    // of current_delta, but determining that upper limit would be
    // difficult. At worse, the loop is executed more times than
    // necessary. But no infinite loop can result since there is
    // an upper bound on the loop variable.
    if (current_delta == 0.0)
      break;
    // previous_delta != 0 or would have been kicked out
    // previously
    double r = log(current_delta) / log(previous_delta);

    if (r > 1.9 && r < 2.1) {
      // If convergence theory applied perfectly, r would be 2 in
      // the convergence region. r close to 2 is good enough. We
      // expect the difference between this integral estimate and
      // the next one to be roughly delta^2.
      error_estimate = current_delta * current_delta;
    } else {
      // Not in the convergence region. Assume only that error is
      // decreasing.
      error_estimate = current_delta;
    }

    if (error_estimate < tolerance)
      break;
  }

  if (level == num_levels) {
    check_less_or_equal("de_integrator",
                        "max integration level reached, but error estimate",
                        error_estimate * c, tolerance * c);
  } else {
    check_less_or_equal("de_integrator", "error estimate", error_estimate * c,
                        tolerance * c);
  }

  num_function_evaluations = 2 * i - 1;
  error_estimate *= c;
  return c * integral;
}
}  // namespace math
}  // namespace stan
#endif
