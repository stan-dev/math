#ifndef STAN_MATH_REV_MAT_FUN_ODE_VARI_HPP
#define STAN_MATH_REV_MAT_FUN_ODE_VARI_HPP

#include <stan/math/rev/mat/functor/cvodes_ode_data.hpp>
#include <stan/math/prim/scal/meta/return_type.hpp>
#include <cvodes/cvodes.h>
#include <cvodes/cvodes_spils.h>
#include <nvector/nvector_serial.h>
#include <ostream>
#include <vector>

namespace stan {
namespace math {
template <typename T_ode_data>
class cvodes_adj_ode_alloc : public chainable_alloc {
 public:
  void* cvodes_mem_;
  T_ode_data* cvodes_data_;

  explicit cvodes_adj_ode_alloc(void* cvodes_mem, T_ode_data* cvodes_data)
      : cvodes_mem_(cvodes_mem), cvodes_data_(cvodes_data) {}

  ~cvodes_adj_ode_alloc() {
    CVodeFree(&cvodes_mem_);
    delete cvodes_data_;
  }
};

/*
 * ode_vari is the special vari that will handle running the adjoint ODE
 * for all the other varis
 *
 * @tparam T_ode_data Type of cvodes_ode_data structure
 */
template <typename T_ode_data>
class cvodes_adj_ode_vari : public vari {
 protected:
  double t0_;
  int N_;
  int M_;
  double relative_tolerance_;
  double absolute_tolerance_;
  vari** initial_v_;
  vari** theta_v_;
  std::vector<double> ts_;
  vari** non_chaining_varis_;
  void* cvodes_mem_;
  T_ode_data* cvodes_data_;
  cvodes_adj_ode_alloc<T_ode_data>* alloc_;

 public:
  /*
   * ode_vari is the special vari that handles running the adjoint ODE.
   * For each ODE solve, there is only one ode_vari
   *
   * @param t0 Initial time for forward ODE solve
   * @param relative_tolerance Relative tolerance of adjoint ODE solve
   * @param absolute_tolerance Absolute tolerance of adjoint ODE solve
   * @param y Output of forward ODE solve
   * @param y0_ Initial conditions for forward ODE solve
   * @param theta_ Parameters for forward ODE solve
   * @param ts Time steps to produce output in forward ODE solve
   * @param cvodes_mem Pointer to CVODES internal memory space
   * @param cvodes_data Pointer to cvodes_ode_data struct that computes
   *   jacobians and such
   */
  explicit cvodes_adj_ode_vari(double t0, double relative_tolerance,
                               double absolute_tolerance, double value,
                               const std::vector<var>& initial,
                               const std::vector<var>& theta,
                               const std::vector<double>& ts,
                               vari** non_chaining_varis, void* cvodes_mem,
                               T_ode_data* cvodes_data)
      : vari(value),
        t0_(t0),
        N_(initial.size()),
        M_(theta.size()),
        relative_tolerance_(relative_tolerance),
        absolute_tolerance_(absolute_tolerance),
        initial_v_(ChainableStack::instance().memalloc_.alloc_array<vari*>(
            initial.size())),
        theta_v_(ChainableStack::instance().memalloc_.alloc_array<vari*>(
            theta.size())),
        ts_(ts),
        non_chaining_varis_(non_chaining_varis),
        cvodes_mem_(cvodes_mem),
        cvodes_data_(cvodes_data),
        // I'm not 100% sure I'm allocing alloc_ right
        alloc_(new cvodes_adj_ode_alloc<T_ode_data>(cvodes_mem, cvodes_data)) {
    // These are the input parameters, so we'll need to increment these
    // adjoints
    for (size_t i = 0; i < N_; i++)
      initial_v_[i] = initial[i].vi_;
    for (size_t i = 0; i < M_; i++)
      theta_v_[i] = theta[i].vi_;
  }

  virtual void chain() {
    // std::cout << "chain" << std::endl; <-- Good way to verify it's only
    //  being called once
    N_VConst(0.0, cvodes_data_->nv_state_sens_adj_);

    try {
      int indexB;
      // This is all boilerplate CVODES setting up the adjoint ODE to solve
      // CV_ADAMS seemed to work better than CV_BDF on the toy problem
      // I was playing with.
      cvodes_check_flag(CVodeCreateB(cvodes_mem_, CV_BDF, CV_NEWTON, &indexB),
                        "CVodeCreateB");

      cvodes_check_flag(
          CVodeSetUserDataB(cvodes_mem_, indexB,
                            reinterpret_cast<void*>(cvodes_data_)),
          "CVodeSetUserDataB");

      // The ode_rhs_adj_sense functions passed in here cause problems with
      // the autodiff stack (they can cause reallocations of the internal
      // vectors and cause segfaults)
      cvodes_check_flag(
          CVodeInitB(cvodes_mem_, indexB, &T_ode_data::ode_rhs_adj_sens,
                     ts_.back(), cvodes_data_->nv_state_sens_adj_),
          "CVodeInitB");

      cvodes_check_flag(
          CVodeSStolerancesB(cvodes_mem_, indexB, relative_tolerance_,
                             absolute_tolerance_),
          "CVodeSStolerancesB");

      // cvodes_check_flag(CVSpilsSetLinearSolverB(cvodes_mem_, indexB,
      // cvodes_data_->LS_adj_), "CVSpilsSetLinearSolverB");
      cvodes_check_flag(
          CVDlsSetLinearSolverB(cvodes_mem_, indexB, cvodes_data_->LS_adj_,
                                cvodes_data_->A_adj_),
          "CVDlsSetLinearSolverB");

      // At every time step, collect the adjoints from the output
      // variables and re-initialize the solver
      for (int i = ts_.size() - 1; i >= 0; --i) {
        // Take in the adjoints from all the output variables at this point
        // in time
        for (int j = 0; j < N_; j++) {
          // std::cout << "j " << j << std::endl;
          if (i == ts_.size() - 1 && j == N_ - 1) {
            NV_Ith_S(cvodes_data_->nv_state_sens_adj_, j) += adj_;
          } else {
            NV_Ith_S(cvodes_data_->nv_state_sens_adj_, j)
                += non_chaining_varis_[i * N_ + j]->adj_;
          }
        }

        cvodes_check_flag(CVodeReInitB(cvodes_mem_, indexB, ts_[i],
                                       cvodes_data_->nv_state_sens_adj_),
                          "CVodeB");

        // std::cout << i << std::endl;
        cvodes_check_flag(
            CVodeB(cvodes_mem_, (i > 0) ? ts_[i - 1] : t0_, CV_NORMAL),
            "CVodeB");

        // Currently unused, but I should probably check it. CVodeGetB sets
        // it to the time that the ode successfully integrated back to.
        // This should be equal to (i > 0) ? ts_[i - 1] : t0_
        double tret;

        cvodes_check_flag(CVodeGetB(cvodes_mem_, indexB, &tret,
                                    cvodes_data_->nv_state_sens_adj_),
                          "CVodeGetB");
      }

      // After integrating all the way back to t0, we finally have the
      // the adjoints we wanted
      // These are the dlog_density / d(initial_conditions[s]) adjoints
      for (size_t s = 0; s < N_; s++) {
        initial_v_[s]->adj_ += NV_Ith_S(cvodes_data_->nv_state_sens_adj_, s);
      }

      // These are the dlog_density / d(parameters[s]) adjoints
      for (size_t s = 0; s < M_; s++) {
        theta_v_[s]->adj_ += NV_Ith_S(cvodes_data_->nv_state_sens_adj_, N_ + s);
      }
    } catch (const std::exception& e) {
      throw;
    }
  }
};

/*
 * build_ode_varis takes the CVODES output and wraps it in a 2D Stan array
 * of vars. To avoiding chain for every individual output, we create
 * one giant ode_vari which does the actual adjoint calculation and then
 * num_timesteps * num_states - 1 varis which are just placeholders.
 *
 * Because all the varis are pushed onto the stack at once, whenever any of
 * them are called, all the adjoints at these outputs are ready to be
 * chained back up.
 *
 * By putting all but one vari on the autodiff stack, only one vari gets
 * called and one one backwards ODE solve needs to happen
 *
 * @tparam T_ode_data Convenience template to hide true type of cvodes_data
 * @tparam T_initial Type of initial conditions, can be double or var
 * @tparam T_param Type of parameters, can be double or var
 * @param t0 Initial time for forward ODE solve
 * @param relative_tolerance Relative tolerance of adjoint ODE solve
 * @param absolute_tolerance Absolute tolerance of adjoint ODE solve
 * @param y Output of forward ODE solve
 * @param y0_ Initial conditions for forward ODE solve
 * @param theta_ Parameters for forward ODE solve
 * @param ts Time steps to produce output in forward ODE solve
 * @param cvodes_mem Pointer to CVODES internal memory space
 * @param cvodes_data Pointer to cvodes_ode_data struct that computes
 *   jacobians and such
 */
template <typename T_ode_data, typename T_initial, typename T_param>
inline std::vector<
    std::vector<typename stan::return_type<T_initial, T_param>::type> >
build_cvodes_adj_ode_vari(double t0, double relative_tolerance,
                          double absolute_tolerance,
                          const std::vector<std::vector<double> >& y,
                          const std::vector<T_initial>& y0_,
                          const std::vector<T_param>& theta_,
                          const std::vector<double>& ts, void* cvodes_mem,
                          T_ode_data* cvodes_data) {
  using std::vector;
  const size_t N = y0_.size();
  vector<var> y0(y0_.begin(), y0_.end());
  vector<var> theta(theta_.begin(), theta_.end());
  vector<vector<var> > y_return(y.size(), vector<var>(N, 0));
  vari** non_chaining_varis
      = ChainableStack::instance().memalloc_.alloc_array<vari*>(y.size() * N
                                                                - 1);
  for (size_t i = 0; i < y.size(); i++) {
    for (size_t j = 0; j < N; j++) {
      // The special vari that will handle the adjoint solve corresponds to
      // the output y[y.size() - 1][N - 1]
      if (i == y.size() - 1 && j == N - 1)
        continue;
      // non_chaining_varis[i * N + j] corresponds to the vari attached to
      // the ode output at time t[i] and state j
      non_chaining_varis[i * N + j] = new vari(y[i][j], false);
    }
  }
  cvodes_adj_ode_vari<T_ode_data>* ode_vari_
      = new cvodes_adj_ode_vari<T_ode_data>(
          t0, relative_tolerance, absolute_tolerance, y[y.size() - 1][N - 1],
          y0, theta, ts, non_chaining_varis, cvodes_mem, cvodes_data);
  for (size_t i = 0; i < y.size(); i++)
    for (size_t j = 0; j < N; j++)
      // Inject our special vari at y[y.size() - 1][N - 1]
      if (i == y.size() - 1 && j == N - 1)
        y_return[i][j] = var(ode_vari_);
      else
        y_return[i][j] = var(non_chaining_varis[i * N + j]);
  return y_return;
}

/*
 * If theta and y are both doubles, just pass the values through (there's
 * no autodiff to handle here).
 */
template <typename T_ode_data>
inline std::vector<std::vector<double> > build_cvodes_adj_ode_vari(
    double t0, double relative_tolerance, double absolute_tolerance,
    const std::vector<std::vector<double> >& y, const std::vector<double>& y0,
    const std::vector<double>& theta, const std::vector<double>& ts,
    void* cvodes_mem, T_ode_data* cvodes_data) {
  return y;
}
}  // namespace math
}  // namespace stan
#endif
