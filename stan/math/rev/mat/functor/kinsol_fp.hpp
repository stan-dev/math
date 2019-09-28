#ifndef STAN_MATH_REV_MAT_FUNCTOR_FP_SOLVER_HPP
#define STAN_MATH_REV_MAT_FUNCTOR_FP_SOLVER_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/mat/fun/to_array_1d.hpp>
#include <stan/math/prim/mat/fun/to_vector.hpp>
#include <stan/math/prim/mat/fun/value_of.hpp>
#include <stan/math/prim/mat/fun/mdivide_left.hpp>
#include <stan/math/prim/scal/err/check_nonnegative.hpp>
#include <stan/math/rev/mat/functor/algebra_system.hpp>
#include <stan/math/rev/mat/functor/kinsol_data.hpp>
#include <stan/math/rev/mat/functor/algebra_system.hpp>
#include <stan/math/prim/mat/err/check_flag_sundials.hpp>

#include <kinsol/kinsol.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <nvector/nvector_serial.h>

#include <iostream>
#include <string>
#include <vector>

namespace stan {
namespace math {

  /**
   * KINSOL algebraic system data holder.
   *
   * @tparam F1 functor type for system function.
   * @tparam F2 functor type for jacobian function. Default is 0.
   *         If 0, use rev mode autodiff to compute the Jacobian.
   */
  template <typename F>
  struct KinsolFixedPointEnv {
    const F& f_;
    const Eigen::VectorXd& y_;
    const size_t N_;
    const size_t M_;
    const std::vector<double>& x_r_;
    const std::vector<int>& x_i_;
    std::ostream* msgs_;

    void* mem_;
    N_Vector nv_x_;
    N_Vector nv_u_scal_;
    N_Vector nv_f_scal_;

    /* Constructor */
    KinsolFixedPointEnv(const F& f, const Eigen::VectorXd& x,
                       const Eigen::VectorXd& y, const std::vector<double>& x_r,
                       const std::vector<int>& x_i, std::ostream* msgs,
                       const std::vector<double>& u_scale,
                       const std::vector<double>& f_scale)
      : f_(f),
        y_(y),
        x_r_(x_r),
        x_i_(x_i),
        msgs_(msgs),
        N_(x.size()),
        M_(y.size()),
        mem_(KINCreate()),
        nv_x_(N_VNew_Serial(N_)),
        nv_u_scal_(N_VNew_Serial(N_)),
        nv_f_scal_(N_VNew_Serial(N_))
    {
      for (int i = 0; i < N_; ++i) {
        NV_Ith_S(nv_x_, i) = x(i);
        NV_Ith_S(nv_u_scal_, i) = u_scale[i];
        NV_Ith_S(nv_f_scal_, i) = f_scale[i];
      }
    }

    ~KinsolFixedPointEnv() {
      N_VDestroy_Serial(nv_x_);
      N_VDestroy_Serial(nv_u_scal_);
      N_VDestroy_Serial(nv_f_scal_);
      KINFree(&mem_);
    }

    /* Implements the user-defined function passed to KINSOL. */
    static int kinsol_f_system(N_Vector x, N_Vector f, void* user_data) {
      auto g = static_cast<const KinsolFixedPointEnv<F>*>(user_data);
      Eigen::VectorXd x_eigen(Eigen::Map<Eigen::VectorXd>(NV_DATA_S(x), g->N_));
      Eigen::Map<Eigen::VectorXd>(N_VGetArrayPointer(f), g->N_)
        = g->f_(x_eigen, g->y_, g->x_r_, g->x_i_, g->msgs_);
      return 0;
    }
  };

  struct FixedPointSolver {

    template <typename F>
    void kinsol_solve_fp(Eigen::VectorXd& x,
                         KinsolFixedPointEnv<F>& env, 
                         double f_tol,
                         long int max_num_steps) {
      int N = env.N_;
      void* mem = env.mem_;

      const int default_anderson_depth = 4;
      int anderson_depth = std::min(N, default_anderson_depth);

      check_flag_sundials(KINSetNumMaxIters(mem, max_num_steps),
                          "KINSetNumMaxIters");
      check_flag_sundials(KINSetMAA(mem, anderson_depth),
                          "KINSetMAA");
      check_flag_sundials(KINInit(mem, &env.kinsol_f_system, env.nv_x_), "KINInit");

      check_flag_sundials(KINSetFuncNormTol(mem, f_tol),
                          "KINSetFuncNormTol");

      check_flag_sundials(KINSetUserData(mem, static_cast<void*>(&env)),
                          "KINSetUserData");

      KINSol(mem, env.nv_x_, KIN_FP, env.nv_u_scal_, env.nv_f_scal_);

      long int nniter;
      double norm;
      KINGetNumNonlinSolvIters(mem, &nniter);
      KINGetFuncNorm(mem, &norm);

      for (int i = 0; i < N; ++i) {
        x(i) = NV_Ith_S(env.nv_x_, i);
      }
    }

    template <typename F, typename T1>
    Eigen::Matrix<double, -1, 1>
    solve(const Eigen::Matrix<T1, -1, 1>& x,
          const Eigen::Matrix<double, -1, 1>& y,
          KinsolFixedPointEnv<F>& env,
          double f_tol,
          long int max_num_steps) {
      Eigen::VectorXd xd(stan::math::value_of(x));
      kinsol_solve_fp(xd, env, f_tol, max_num_steps);
      return xd;
    }

    template <typename F, typename T1>
    Eigen::Matrix<stan::math::var, -1, 1>
    solve(const Eigen::Matrix<T1, -1, 1>& x,
          const Eigen::Matrix<stan::math::var, -1, 1>& y,
          KinsolFixedPointEnv<F>& env,
          double f_tol,
          long int max_num_steps) {
      using stan::math::value_of;
      using stan::math::var;

      // FP solution
      Eigen::VectorXd xd(value_of(x));
      kinsol_solve_fp(xd, env, f_tol, max_num_steps);

      auto g = [&env](const Eigen::Matrix<var, -1, 1>& par_) {
        Eigen::Matrix<var, -1, 1> x_(par_.head(env.N_));
        Eigen::Matrix<var, -1, 1> y_(par_.tail(env.M_));
        return env.f_(x_, y_, env.x_r_, env.x_i_, env.msgs_);
      };

      Eigen::VectorXd theta(x.size() + y.size());
      for (int i = 0; i < env.N_; ++i) {
        theta(i) = xd(i);
      }
      for (int i = 0; i < env.M_; ++i) {
        theta(i + env.N_) = env.y_(i);
      }
      Eigen::Matrix<double, -1, 1> fx;
      Eigen::Matrix<double, -1, -1> J;
      stan::math::jacobian(g, theta, fx, J);
      Eigen::MatrixXd A(J.block(0, 0, env.N_, env.N_));
      Eigen::MatrixXd b(J.block(0, env.N_, env.N_, env.M_));
      A = Eigen::MatrixXd::Identity(env.N_, env.N_) - A;
      Eigen::MatrixXd Jxy = A.colPivHouseholderQr().solve(b);
      std::vector<double> gradients(env.M_);
      Eigen::Matrix<var, -1, 1> x_sol_var(env.N_);
      for (int i = 0; i < env.N_; ++i) {
        gradients = stan::math::to_array_1d(Eigen::VectorXd(Jxy.row(i)));
        x_sol_var[i] = stan::math::precomputed_gradients(xd(i), stan::math::to_array_1d(y), gradients);
      }
      return x_sol_var;
    }
  };
}
}

#endif
