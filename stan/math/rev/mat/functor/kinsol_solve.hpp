#ifndef STAN_MATH_REV_MAT_FUNCTOR_KINSOL_SOLVE_HPP
#define STAN_MATH_REV_MAT_FUNCTOR_KINSOL_SOLVE_HPP

#include <stan/math/prim/mat/fun/to_array_1d.hpp>
#include <stan/math/prim/mat/fun/to_vector.hpp>
#include <stan/math/rev/mat/functor/kinsol_data.hpp>
#include <stan/math/rev/mat/functor/algebra_system.hpp>
#include <stan/math/prim/mat/err/check_flag.hpp>

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
   * @param[in] f Functor that evaluated the system of equations.
   * @param[in] x Vector of starting values.
   * @param[in] y Parameter vector for the equation system. The function
   *            is overloaded to treat y as a vector of doubles or of a
   *            a template type T.
   * @param[in] dat Continuous data vector for the equation system.
   * @param[in] dat_int Integer data vector for the equation system.
   * @param[in, out] msgs The print stream for warning messages.
   * @param[in] scaling_step_size Scaled-step stopping tolerance. If
   *            a Newton step is smaller than the scaling step
   *            tolerance, the code breaks, assuming the solver is no
   *            longer making significant progress (i.e. is stuck)
   * @param[in] function_tolerance determines whether roots are acceptable.
   * @param[in] max_num_steps Maximum number of function evaluations.
   * @param[in] custom_jacobian If 0, use Kinsol's default to compute the
   *            jacobian for the Newton step, namely Quotient Difference
   *            (finite difference). If 1, use reverse-mode AD, unless
   *            the user specifies their own method.
   * @param[in] J_F A functor that computes the Jacobian for the Newton step.
   *            Defaults to reverse-mode AD.
   * @param[in] steps_eval_jacobian Maximum number of steps before the
   *            Jacobian gets recomputed. Note that Kinsol's default is 10.
   *            If equal to 1, the algorithm computes exact Newton steps.
   * @param[in] global_line_search does the solver use a global line search?
   *            If equal to KIN_NONE, no, if KIN_LINESEARCH, yes.
   * @return x_solution Vector of solutions to the system of equations.
   * @throw <code>std::invalid_argument</code> if Kinsol returns a negative
   *        flag when setting up the solver.
   * @throw <code>boost::math::evaluation_error</code> if Kinsol returns a
   *        negative flag after attempting to solve the equation.
   */
  template <typename F1, typename F2 = kinsol_J_f>
  Eigen::VectorXd
  kinsol_solve(const F1& f,
               const Eigen::VectorXd& x,
               const Eigen::VectorXd& y,
               const std::vector<double>& dat,
               const std::vector<int>& dat_int,
               std::ostream* msgs = nullptr,
               double scaling_step_tol = 1e-3,
               double function_tolerance = 1e-6,
               long int max_num_steps = 1e+3,  // NOLINT(runtime/int)
               bool custom_jacobian = 1,
               const F2& J_f = kinsol_J_f(),
               int steps_eval_jacobian = 10,
               int global_line_search = KIN_LINESEARCH) {
    int N = x.size();
    typedef kinsol_system_data<F1, F2> system_data;
    system_data kinsol_data(f, J_f, x, y, dat, dat_int, msgs);

    void* kinsol_memory = KINCreate();

    check_flag(KINInit(kinsol_memory, &system_data::kinsol_f_system,
                       kinsol_data.nv_x_), "KINInit");

    N_Vector scaling = N_VNew_Serial(N);
    N_VConst_Serial(1.0, scaling);  // no scaling

    check_flag(KINSetFuncNormTol(kinsol_memory, function_tolerance),
               "KINSetFuncNormTol");
    check_flag(KINSetScaledStepTol(kinsol_memory, scaling_step_tol),
               "KINSetScaledStepTol");
    check_flag(KINSetMaxSetupCalls(kinsol_memory, steps_eval_jacobian),
               "KINSetMaxSetupCalls");

    // CHECK
    // The default value is 1000 * ||u_0||_D where ||u_0|| is the initial guess.
    // So we run into issues if ||u_0|| = 0.
    // If the norm is non-zero, use kinsol's default (accessed with 0),
    // else use the dimension of x -- CHECK - find optimal length.
    double max_newton_step = (x.norm() == 0) ? x.size() : 0;
    check_flag(KINSetMaxNewtonStep(kinsol_memory, max_newton_step),
               "KINSetMaxNewtonStep");
    check_flag(KINSetUserData(kinsol_memory,
                              reinterpret_cast<void*>(&kinsol_data)),
                 "KINSetUserData");

    // construct Linear solver
    check_flag(KINSetLinearSolver(kinsol_memory, kinsol_data.LS_,
                                  kinsol_data.J_),
               "KINSetLinearSolver");

    if (custom_jacobian) check_flag(KINSetJacFn(kinsol_memory,
                                    &system_data::kinsol_jacobian),
                                    "KINSetJacFn");

    N_Vector nv_x = N_VNew_Serial(N);
    realtype* nv_x_data = N_VGetArrayPointer_Serial(nv_x);
    for (int i = 0; i < N; i++) nv_x_data[i] = x(i);

    check_flag_kinsol(KINSol(kinsol_memory, nv_x,
                             global_line_search, scaling, scaling),
                      max_num_steps);

    KINFree(&kinsol_memory);

    Eigen::VectorXd x_solution(N);
    for (int i = 0; i < N; i++) x_solution(i) = nv_x_data[i];

    return x_solution;
  }

}  // namespace math
}  // namespace stan
#endif
