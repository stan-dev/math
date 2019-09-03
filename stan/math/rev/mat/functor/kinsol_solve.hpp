#ifndef STAN_MATH_REV_MAT_FUNCTOR_KINSOL_SOLVE_HPP
#define STAN_MATH_REV_MAT_FUNCTOR_KINSOL_SOLVE_HPP

#include <stan/math/prim/mat/fun/to_array_1d.hpp>
#include <stan/math/prim/mat/fun/to_vector.hpp>
#include <stan/math/rev/mat/functor/kinsol_data.hpp>
#include <stan/math/rev/mat/functor/algebra_system.hpp>

#include <kinsol/kinsol.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <nvector/nvector_serial.h>

#include <vector>

namespace stan {
namespace math {
  /**
   * Return the solution to the specified algebraic system,
   * given an initial guess.
   * 
   * . . .
   * @param[in] function_tolerance determines how small ||f(x)|| needs to be.
   * @param[in] max_num_steps maximum number of iterations.
   * @param[in] scaling_step_tol if a Newton step is smaller than the scaling
   *            step tolerance, the code breaks, assuming the solver is no
   *            longer making significant progress (i.e. is stuck).
   * @param[in] custom_jacobian. If 0, use Kinsol's quotient differentiation to
   *            do the linear solve. If 1, either use reverse-mode autodiff, or
   *            a method specified by the user.
   * @param[in] J_f user supplied method for computing the Jacobian of f
   *            w.r.t x. Defaults to reverse mode autodiff.
   * @param[in] steps_eval_jacobian maximum number of steps before the
   *            Jacobian gets recomputed. Note that Kinsol's default is 10.
   *            If equal to 1, the algorithm computes exact Newton steps.
   * @param[in] global_line_search does the solver use a global line search?
   *            If equal to KIN_NONE, no, if KIN_LINESEARCH, yes.  
   */
  template <typename F1, typename F2 = kinsol_J_f>
  Eigen::VectorXd 
  kinsol_solve(const F1& f,
               const Eigen::VectorXd& x,
               const Eigen::VectorXd& y,
               const std::vector<double>& dat,
               const std::vector<int>& dat_int,
               std::ostream* msgs = nullptr,
               double function_tolerance = 1e-6,
               long int max_num_steps = 1e+3,
               double scaling_step_tol = 1e-3,
               bool custom_jacobian = 1,  // TEST - should be 0.
               const F2& J_f = kinsol_J_f(),
               int steps_eval_jacobian = 10,
               int global_line_search = KIN_LINESEARCH) {
    int N = x.size();
    typedef kinsol_system_data<F1, F2> system_data;
    system_data kinsol_data(f, J_f, x, y, dat, dat_int, msgs);

    void* kinsol_memory = KINCreate();

    int flag;  // TO DO -- replace with a checkflag procedure
    flag = KINInit(kinsol_memory, &system_data::kinsol_f_system,
                   kinsol_data.nv_x_);

    // TO DO - construct a set option procedure
    N_Vector scaling = N_VNew_Serial(N);
    N_VConst_Serial(1.0, scaling);  // no scaling
    flag = KINSetFuncNormTol(kinsol_memory, function_tolerance);
    flag = KINSetScaledStepTol(kinsol_memory, scaling_step_tol);
    flag = KINSetMaxSetupCalls(kinsol_memory, steps_eval_jacobian);

    // FIX ME
    // The default value is 1000 * ||u_0||_D where ||u_0|| is the initial guess.
    // So we run into issues if ||u_0|| = 0.
    // If the norm is non-zero, use kinsol's default (accessed with 0),
    // else use the dimension of x -- CHECK - find optimal length.
    double max_newton_step = (x.norm() == 0) ? x.size() : 0;
    flag = KINSetMaxNewtonStep(kinsol_memory, max_newton_step);  // scaled length of Newton step

    flag = KINSetUserData(kinsol_memory,
                          reinterpret_cast<void*>(&kinsol_data));

    // construct Linear solver
    flag = KINSetLinearSolver(kinsol_memory, kinsol_data.LS_, kinsol_data.J_);

    // FOR TESTS: comment this out to use Kinsol's default methods for Jacobian,
    // i.e. finite differentiation.
    if (!custom_jacobian) flag = KINSetJacFn(kinsol_memory, 0);
    else flag = KINSetJacFn(kinsol_memory, &system_data::kinsol_jacobian);

    // TO DO - a better way to do this conversion.
    N_Vector nv_x = N_VNew_Serial(N);
    realtype* nv_x_data = N_VGetArrayPointer_Serial(nv_x);
    for (int i = 0; i < N; i++) nv_x_data[i] = x(i);

    flag = KINSol(kinsol_memory, nv_x,
                  global_line_search, scaling, scaling);

    // TO DO - return an exception when the flag is negative.
    // std::cout << "Kinsol flag: " << flag << std::endl;
    if (flag < 0) {
      std::ostringstream message;
      message << "lgp_solver: the kinsol solver encounter an error"
              << " with flag = " << flag;
      throw boost::math::evaluation_error(message.str());
    }

    // keep track of how many iterations are used.
    // Useful when running tests.
    // long int nniters;
    // KINGetNumNonlinSolvIters(kinsol_memory, &nniters);
    // std::cout << "number of iterations: " << nniters << std::endl;

    KINFree(&kinsol_memory);

    // CHECK - avoid / simplifies this conversion step?
    Eigen::VectorXd x_solution(N);
    for (int i = 0; i < N; i++) x_solution(i) = nv_x_data[i];

    return x_solution;
  }

}  // namespace math
}  // namespace stan
#endif
