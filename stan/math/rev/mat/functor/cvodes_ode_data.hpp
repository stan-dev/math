#ifndef STAN_MATH_REV_MAT_FUNCTOR_CVODES_ODE_DATA_HPP
#define STAN_MATH_REV_MAT_FUNCTOR_CVODES_ODE_DATA_HPP

#include <stan/math/rev/mat/functor/ode_system.hpp>
#include <stan/math/rev/scal/meta/is_var.hpp>
#include <cvodes/cvodes.h>
#include <cvodes/cvodes_band.h>
#include <cvodes/cvodes_dense.h>
#include <nvector/nvector_serial.h>
#include <algorithm>
#include <vector>

namespace stan {
namespace math {

/**
 * CVODES ode data holder object which is used during CVODES
 * integration for CVODES callbacks.
 *
 * @tparam F type of functor for the base ode system.
 * @tparam T_initial type of initial values
 * @tparam T_param type of parameters
 */

template <typename F, typename T_initial, typename T_param>
class cvodes_ode_data {
  const F& f_;
  const std::vector<T_initial>& y0_;
  const std::vector<T_param>& theta_;
  const std::vector<double> theta_dbl_;
  const size_t N_;
  const size_t M_;
  const size_t param_var_ind_;
  const std::vector<double>& x_;
  const std::vector<int>& x_int_;
  std::ostream* msgs_;
  const coupled_ode_system<F, T_initial, T_param> coupled_ode_;
  const ode_system<F> ode_system_;

  typedef cvodes_ode_data<F, T_initial, T_param> ode_data;
  typedef stan::is_var<T_initial> initial_var;
  typedef stan::is_var<T_param> param_var;

 public:
  std::vector<double> coupled_state_;
  const size_t S_;
  N_Vector nv_state_;
  N_Vector* nv_state_sens_;
  
  /**
   * Construct CVODES ode data object to enable callbacks from
   * CVODES during ODE integration. Static callbacks are defined
   * for the ODE RHS (<code>ode_rhs</code>), the ODE sensitivity
   * RHS (<code>ode_rhs_sens</code>) and for the ODE Jacobian wrt
   * to the states (<code>dense_jacobian</code>).
   *
   * @param[in] f ode functor.
   * @param[in] y0 initial state of the base ode.
   * @param[in] theta parameters of the base ode.
   * @param[in] x continuous data vector for the ODE.
   * @param[in] x_int integer data vector for the ODE.
   * @param[in] msgs stream to which messages are printed.
   */
  cvodes_ode_data(const F& f, const std::vector<T_initial>& y0,
                  const std::vector<T_param>& theta,
                  const std::vector<double>& x, const std::vector<int>& x_int,
                  std::ostream* msgs)
      : f_(f),
        y0_(y0),
        theta_(theta),
        theta_dbl_(value_of(theta)),
        N_(y0.size()),
        M_(theta.size()),
        param_var_ind_(initial_var::value ? N_ : 0),
        x_(x),
        x_int_(x_int),
        msgs_(msgs),
        coupled_ode_(f, y0, theta, x, x_int, msgs),
        ode_system_(f, value_of(theta), x, x_int, msgs),
        coupled_state_(coupled_ode_.initial_state()),
        S_((initial_var::value ? N_ : 0) + (param_var::value ? M_ : 0)),
        nv_state_(N_VMake_Serial(N_, &coupled_state_[0])),
        nv_state_sens_(NULL)
  {
    if (S_ != 0) {
      nv_state_sens_ = N_VCloneVectorArrayEmpty_Serial(S_, nv_state_);
      for(std::size_t i=0; i < S_; i++) {
        //nv_state_sens_[i] = N_VMake_Serial(N_, &coupled_state_[N_] + i * N_);
        NV_DATA_S(nv_state_sens_[i]) = &coupled_state_[N_] + i * N_;
      }
    }
  }

  /**/
  ~cvodes_ode_data() {
   N_VDestroy_Serial(nv_state_);
   if (nv_state_sens_ != NULL)
     N_VDestroyVectorArray_Serial(nv_state_sens_, S_);
  }
  /**/

  static int ode_rhs(double t, N_Vector y, N_Vector ydot, void* user_data) {
    const ode_data* explicit_ode = static_cast<const ode_data*>(user_data);
    explicit_ode->rhs(NV_DATA_S(y), NV_DATA_S(ydot), t);
    return 0;
  }

  static int ode_rhs_sens(int Ns, realtype t, N_Vector y, N_Vector ydot,
                          N_Vector* yS, N_Vector* ySdot, void* user_data,
                          N_Vector tmp1, N_Vector tmp2) {
    const ode_data* explicit_ode = static_cast<const ode_data*>(user_data);
    //const std::vector<double> y_vec(NV_DATA_S(y),
    //                                NV_DATA_S(y) + explicit_ode->N_);
    //explicit_ode->rhs_sens(explicit_ode->y0_, explicit_ode->theta_, t, y_vec,
    //                       yS, ySdot);
    explicit_ode->rhs_sens2(t, NV_DATA_S(y), yS, ySdot);
    return 0;
  }

  static int dense_jacobian(long int N,  // NOLINT(runtime/int)
                            realtype t, N_Vector y, N_Vector fy, DlsMat J,
                            void* user_data, N_Vector tmp1, N_Vector tmp2,
                            N_Vector tmp3) {
    const ode_data* explicit_ode = static_cast<const ode_data*>(user_data);
    return explicit_ode->dense_jacobian(NV_DATA_S(y), J, t);
  }

 private:
  inline void rhs(const double y[], double dy_dt[], double t) const {
    const std::vector<double> y_vec(y, y + N_);
    const std::vector<double> dy_dt_vec = f_(t, y_vec, theta_dbl_, x_, x_int_, msgs_);
    check_size_match("cvodes_ode_data", "dz_dt", dy_dt_vec.size(),
                     "states", N_);
    std::copy(dy_dt_vec.begin(), dy_dt_vec.end(), dy_dt);
  }

  // the Jacobian of the ODE system is the coupled ode system for
  // varying states evaluated at the state y.
  inline int dense_jacobian(const double* y, DlsMat J, double t) const {
    const std::vector<double> y_vec(y, y + N_);
    std::vector<var> y_vec_var(y_vec.begin(), y_vec.end());
    coupled_ode_system<F,var,double> ode_jacobian(f_, y_vec_var, theta_dbl_, x_, x_int_, msgs_);
    std::vector<double> ode_jacobian_y(ode_jacobian.size(), 0);
    ode_jacobian(ode_jacobian.initial_state(), ode_jacobian_y, t);
    std::copy(ode_jacobian_y.begin() + N_, ode_jacobian_y.end(), J->data);
    //Eigen::Map<Eigen::MatrixXd> Jy_map(J->data, N_, N_);
    //Jy_map = Eigen::Map<Eigen::MatrixXd>(&ode_jacobian_y[0] + N_, N_, N_);
    return 0;
  }

  inline void rhs_sens_initial(const Eigen::MatrixXd& Jy, N_Vector* yS,
                               N_Vector* ySdot) const {
    for (size_t m = 0; m < N_; ++m) {
      Eigen::Map<Eigen::VectorXd> yS_eig(NV_DATA_S(yS[m]), N_);
      Eigen::Map<Eigen::VectorXd> ySdot_eig(NV_DATA_S(ySdot[m]), N_);
      ySdot_eig = Jy * yS_eig;
    }
  }

  inline void rhs_sens_param(const Eigen::MatrixXd& Jy,
                             const Eigen::MatrixXd& Jtheta, N_Vector* yS,
                             N_Vector* ySdot) const {
    using Eigen::Map;
    using Eigen::VectorXd;
    for (size_t m = 0; m < M_; ++m) {
      Map<VectorXd> yS_eig(NV_DATA_S(yS[param_var_ind_ + m]), N_);
      Map<VectorXd> ySdot_eig(NV_DATA_S(ySdot[param_var_ind_ + m]), N_);
      ySdot_eig = Jy * yS_eig + Jtheta.col(m);
    }
  }

  void rhs_sens2(double t, double y[], N_Vector* yS,
                 N_Vector* ySdot) const {
    std::vector<double> z(coupled_state_.size());
    std::vector<double> dz_dt(coupled_state_.size());
    std::copy(y, y + N_, z.begin());
    for(std::size_t s=0; s < S_; s++)
      std::copy(NV_DATA_S(yS[s]), NV_DATA_S(yS[s]) + N_, z.begin() + (s+1) * N_);
    coupled_ode_(z, dz_dt, t);
    for(std::size_t s=0; s < S_; s++)
      std::copy(dz_dt.begin() + (s+1) * N_,  dz_dt.begin() + (s+2) * N_, NV_DATA_S(ySdot[s]));
  }

  /**
   * Calculate the sensitivity RHS for varying initials and parameters.
   *
   * @param[in] initial var vector
   * @param[in] param var vector
   * @param[in] t time
   * @param[in] y state of the base ODE system
   * @param[in] yS array of M N_Vectors of size N, i.e. state of sensitivity
   * RHS
   * @param[out] ySdot array of M N_Vectors of size N of the sensitivity RHS
   */
  void rhs_sens(const std::vector<var>& initial, const std::vector<var>& param,
                double t, const std::vector<double>& y, N_Vector* yS,
                N_Vector* ySdot) const {
    Eigen::VectorXd dy_dt(N_);
    Eigen::MatrixXd Jy(N_, N_);
    Eigen::MatrixXd Jtheta(N_, M_);
    ode_system_.jacobian(t, y, dy_dt, Jy, Jtheta);
    rhs_sens_initial(Jy, yS, ySdot);
    rhs_sens_param(Jy, Jtheta, yS, ySdot);
  }

  /**
   * Calculate the sensitivity RHS for fixed initials and varying parameters.
   *
   * @param[in] initial double vector
   * @param[in] param var vector
   * @param[in] t time
   * @param[in] y state of the base ODE system
   * @param[in] yS array of M N_Vectors of size N, i.e. state of sensitivity
   * RHS
   * @param[out] ySdot array of M N_Vectors of size N of the sensitivity RHS
   */
  void rhs_sens(const std::vector<double>& initial,
                const std::vector<var>& param, double t,
                const std::vector<double>& y, N_Vector* yS,
                N_Vector* ySdot) const {
    Eigen::VectorXd dy_dt(N_);
    Eigen::MatrixXd Jy(N_, N_);
    Eigen::MatrixXd Jtheta(N_, M_);
    ode_system_.jacobian(t, y, dy_dt, Jy, Jtheta);
    rhs_sens_param(Jy, Jtheta, yS, ySdot);
  }

  /**
   * Calculate the sensitivity RHS for varying initials and fixed parameters.
   *
   * @param[in] initial var vector
   * @param[in] param double vector
   * @param[in] t time
   * @param[in] y state of the base ODE system
   * @param[in] yS array of M N_Vectors of size N, i.e. state of sensitivity
   * RHS
   * @param[out] ySdot array of M N_Vectors of size N of the sensitivity RHS
   */
  void rhs_sens(const std::vector<var>& initial,
                const std::vector<double>& param, double t,
                const std::vector<double>& y, N_Vector* yS,
                N_Vector* ySdot) const {
    Eigen::VectorXd dy_dt(N_);
    Eigen::MatrixXd Jy(N_, N_);
    ode_system_.jacobian(t, y, dy_dt, Jy);
    rhs_sens_initial(Jy, yS, ySdot);
  }

  /**
   * Calculate the empty sensitivity RHS.
   *
   * @param[in] initial double vector
   * @param[in] param double vector
   * @param[in] t time
   * @param[in] y state of the base ODE system
   * @param[in] yS array of M N_Vectors of size N, i.e. state of sensitivity
   * RHS
   * @param[out] ySdot array of M N_Vectors of size N of the sensitivity RHS
   */
  void rhs_sens(const std::vector<double>& initial,
                const std::vector<double>& param, double t,
                const std::vector<double>& y, N_Vector* yS,
                N_Vector* ySdot) const {}
};

}  // namespace math
}  // namespace stan
#endif
