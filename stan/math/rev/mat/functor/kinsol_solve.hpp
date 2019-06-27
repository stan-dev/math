#ifndef STAN_MATH_REV_MAT_FUNCTOR_KINSOL_SOLVE_HPP
#define STAN_MATH_REV_MAT_FUNCTOR_KINSOL_SOLVE_HPP

#include <stan/math/prim/mat/fun/to_array_1d.hpp>
#include <stan/math/prim/mat/fun/to_vector.hpp>
// #include <stan/math/prim/mat/functor/Eigen.hpp>
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
   * [...]
   */
  template <typename F>
  Eigen::VectorXd 
  kinsol_solve(const F& f, const Eigen::VectorXd& x,
               const Eigen::VectorXd& y, const std::vector<double>& dat,
               const std::vector<int>& dat_int, std::ostream* msgs = nullptr,
               double function_tolerance = 1e-6,
               long int max_num_steps = 1e+3,
               int global_line_search = KIN_NONE, 
               int steps_eval_jacobian = 1,
               double scaling_step_tol = 1e-5) {
    // CHECK -- what tuning parameters do we want to include?
    // E.g scaling_step_tol, scaling, etc.
    int N = x.size();
    typedef kinsol_system_data<F> system_data;
    system_data kinsol_data(f, x, y, dat, dat_int, 0);

    void* kinsol_memory = KINCreate();

    int flag;  // FIX ME -- replace with a checkflag procedure
    flag = KINInit(kinsol_memory, &system_data::kinsol_f_system,
                   kinsol_data.nv_x_);

    // FIX ME - construct a set option procedure
    N_Vector scaling = N_VNew_Serial(N);
    N_VConst_Serial(1.0, scaling);  // no scaling
    flag = KINSetFuncNormTol(kinsol_memory, function_tolerance);
    flag = KINSetScaledStepTol(kinsol_memory, scaling_step_tol);
    flag = KINSetMaxSetupCalls(kinsol_memory, steps_eval_jacobian);

    flag = KINSetUserData(kinsol_memory,
                          reinterpret_cast<void*>(&kinsol_data));

    // construct Linear solver
    flag = KINSetLinearSolver(kinsol_memory, kinsol_data.LS_, kinsol_data.J_);
    flag = KINSetJacFn(kinsol_memory, &system_data::kinsol_jacobian);

    // CHECK - a better way to do this conversion.
    N_Vector nv_x = N_VNew_Serial(N);
    realtype* nv_x_data = N_VGetArrayPointer_Serial(nv_x);
    for (int i = 0; i < N; i++) nv_x_data[i] = x(i);

    flag = KINSol(kinsol_memory, nv_x,
                  global_line_search, scaling, scaling);

    KINFree(&kinsol_memory);

    // CHECK - avoid / simplifies this conversion step?
    Eigen::VectorXd x_solution(N);
    for (int i = 0; i < N; i++) x_solution(i) = nv_x_data[i];

    return x_solution;
  }

}  // namespace math
}  // namespace stan
#endif
