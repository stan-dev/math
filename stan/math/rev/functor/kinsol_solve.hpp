#ifndef STAN_MATH_REV_FUNCTOR_KINSOL_SOLVE_HPP
#define STAN_MATH_REV_FUNCTOR_KINSOL_SOLVE_HPP

#include <stan/math/rev/functor/kinsol_data.hpp>
#include <stan/math/rev/functor/algebra_system.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/to_array_1d.hpp>
#include <stan/math/prim/fun/to_vector.hpp>
#include <sundials/sundials_context.h>
#include <kinsol/kinsol.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <nvector/nvector_serial.h>
#include <vector>

namespace stan {
namespace math {

/**
 * Return the solution to the specified algebraic system,
 * given an initial guess. Invokes the Kinsol solver from Sundials.
 *
 * @tparam F type of equation system function.
 * @tparam T type of initial guess vector.
 *
 * @param[in] f Functor that evaluated the system of equations.
 * @param[in] x Vector of starting values.
 * @param[in] y Parameter vector for the equation system. The function
 *            is overloaded to treat y as a vector of doubles or of a
 *            a template type T.
 * @param[in] dat Continuous data vector for the equation system.
 * @param[in] dat_int Integer data vector for the equation system.
 * @param[in, out] msgs The print stream for warning messages.
 * @param[in] scaling_step_tol Scaled-step stopping tolerance. If
 *            a Newton step is smaller than the scaling step
 *            tolerance, the code breaks, assuming the solver is no
 *            longer making significant progress (i.e. is stuck)
 * @param[in] function_tolerance determines whether roots are acceptable.
 * @param[in] max_num_steps Maximum number of function evaluations.
 * @param[in] custom_jacobian If 0, use Kinsol's default to compute the
 *            jacobian for the Newton step, namely Quotient Difference
 *            (finite difference). If 1, use reverse-mode AD, unless
 *            the user specifies their own method.
 * @param[in] J_f A functor that computes the Jacobian for the Newton step.
 *            Defaults to reverse-mode AD.
 * @param[in] steps_eval_jacobian Maximum number of steps before the
 *            Jacobian gets recomputed. Note that Kinsol's default is 10.
 *            If equal to 1, the algorithm computes exact Newton steps.
 * @param[in] global_line_search does the solver use a global line search?
 *            If equal to KIN_NONE, no, if KIN_LINESEARCH, yes.
 * @return x_solution Vector of solutions to the system of equations.
 * @throw <code>std::invalid_argument</code> if Kinsol returns a negative
 *        flag when setting up the solver.
 * @throw <code>std::domain_error</code> if Kinsol fails to solve
 *        equation in max_num_steps iterations.
 * @throw <code>std::runtime_error</code> if Kinsol returns a
 *        negative flag that is not due to hitting max_num_steps.
 */
template <typename F1, typename... Args>
Eigen::VectorXd kinsol_solve(const F1& f, const Eigen::VectorXd& x,
                             const double scaling_step_tol,    // = 1e-3
                             const double function_tolerance,  // = 1e-6
                             const int64_t max_num_steps,      // = 200
                             const bool custom_jacobian,       // = 1
                             const int steps_eval_jacobian,    // = 10
                             const int global_line_search,  // = KIN_LINESEARCH
                             std::ostream* const msgs, const Args&... args) {
  int N = x.size();
  typedef kinsol_system_data<F1, Args...> system_data;
  system_data kinsol_data(f, x, msgs, args...);

  CHECK_KINSOL_CALL(KINInit(kinsol_data.kinsol_memory_,
                            &system_data::kinsol_f_system, kinsol_data.nv_x_));

  N_Vector scaling = N_VNew_Serial(N, kinsol_data.sundials_context_);
  N_Vector nv_x = N_VNew_Serial(N, kinsol_data.sundials_context_);
  Eigen::VectorXd x_solution(N);

  try {
    N_VConst_Serial(1.0, scaling);  // no scaling

    CHECK_KINSOL_CALL(
        KINSetNumMaxIters(kinsol_data.kinsol_memory_, max_num_steps));
    CHECK_KINSOL_CALL(
        KINSetFuncNormTol(kinsol_data.kinsol_memory_, function_tolerance));
    CHECK_KINSOL_CALL(
        KINSetScaledStepTol(kinsol_data.kinsol_memory_, scaling_step_tol));
    CHECK_KINSOL_CALL(
        KINSetMaxSetupCalls(kinsol_data.kinsol_memory_, steps_eval_jacobian));

    // CHECK
    // The default value is 1000 * ||u_0||_D where ||u_0|| is the initial guess.
    // So we run into issues if ||u_0|| = 0.
    // If the norm is non-zero, use kinsol's default (accessed with 0),
    // else use the dimension of x -- CHECK - find optimal length.
    double max_newton_step = (x.norm() == 0) ? x.size() : 0;
    CHECK_KINSOL_CALL(
        KINSetMaxNewtonStep(kinsol_data.kinsol_memory_, max_newton_step));
    CHECK_KINSOL_CALL(KINSetUserData(kinsol_data.kinsol_memory_,
                                     static_cast<void*>(&kinsol_data)));

    // construct Linear solver
    CHECK_KINSOL_CALL(KINSetLinearSolver(kinsol_data.kinsol_memory_,
                                         kinsol_data.LS_, kinsol_data.J_));

    if (custom_jacobian)
      CHECK_KINSOL_CALL(KINSetJacFn(kinsol_data.kinsol_memory_,
                                    &system_data::kinsol_jacobian));

    for (int i = 0; i < N; i++)
      NV_Ith_S(nv_x, i) = x(i);

    kinsol_check(KINSol(kinsol_data.kinsol_memory_, nv_x, global_line_search,
                        scaling, scaling),
                 "KINSol", max_num_steps);

    for (int i = 0; i < N; i++)
      x_solution(i) = NV_Ith_S(nv_x, i);
  } catch (const std::exception& e) {
    N_VDestroy(nv_x);
    N_VDestroy(scaling);
    throw;
  }

  N_VDestroy(nv_x);
  N_VDestroy(scaling);

  return x_solution;
}

}  // namespace math
}  // namespace stan
#endif
