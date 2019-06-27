#ifndef STAN_MATH_REV_MAT_FUNCTOR_KINSOL_DATA_HPP
#define STAN_MATH_REV_MAT_FUNCTOR_KINSOL_DATA_HPP

#include <stan/math/prim/mat/fun/to_array_1d.hpp>
#include <stan/math/prim/mat/fun/to_vector.hpp>

#include <kinsol/kinsol.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <nvector/nvector_serial.h>

#include <vector>

namespace stan {
namespace math {

/**
 * KINSOL algebraic system data holder.
 * Based on cvodes_ode_data.
 * (EXPERIMENTAL)
 * 
 * F: structure type of the functor which returns the system function.
 * CHECK -- do we need a flexible mode for the parameters?
 */
template <typename F>
class kinsol_system_data {
  const F& f_;
  const Eigen::VectorXd& x_;
  const Eigen::VectorXd& y_;
  const size_t N_;
  const std::vector<double>& dat_;
  const std::vector<int>& dat_int_;
  std::ostream* msgs_;

  typedef kinsol_system_data<F> system_data;

public:
  N_Vector nv_x_;
  SUNMatrix J_;
  SUNLinearSolver LS_;

  /**
   * Constructor
   */
  kinsol_system_data(const F& f,
                     const Eigen::VectorXd& x,
                     const Eigen::VectorXd& y,
                     const std::vector<double>& dat,
                     const std::vector<int>& dat_int,
                     std::ostream* msgs)
    : f_(f), x_(x), y_(y), dat_(dat), dat_int_(dat_int),
      msgs_(msgs), N_(x.size()),
      nv_x_(N_VMake_Serial(N_, &to_array_1d(x_)[0])),  // FIX ME - wrap eigen directly
      J_(SUNDenseMatrix(N_, N_)),
      LS_(SUNLinSol_Dense(nv_x_, J_)) { }

  ~ kinsol_system_data() {
      N_VDestroy_Serial(nv_x_);  // FIX ME - will remove this. See above comment.
      SUNLinSolFree(LS_);
      SUNMatDestroy(J_);
  }

  /**
   * Implements the user-defined function passed to KINSOL.
   */
  static int kinsol_f_system (N_Vector x, N_Vector f, void *user_data) {
    const system_data* explicit_system
      = static_cast<const system_data*>(user_data);
    explicit_system->f_system(NV_DATA_S(x), NV_DATA_S(f));

    return 0;
  }

  /**
   * Implements the function of type CVDlsJacFn which is the user-defined
   * callbacks for KINSOL to calculate the jacobian of the root function.
   * The Jacobian is stored in column major format.
   * 
   * REMARK - tmp1 and tmp2 are pointers to memory allocated for variables
   * of type N_Vector which can be used by KINJacFN (the function which
   * computes the Jacobian) as temporary storage or work space.
   * See https://computation.llnl.gov/sites/default/files/public/kin_guide-dev.pdf,
   * page 55.
   */
  static int kinsol_jacobian (N_Vector x, N_Vector f,
                              SUNMatrix J, void *user_data,
                              N_Vector tmp1, N_Vector tmp2) {
    const system_data* explicit_system
      = static_cast<const system_data*>(user_data);
    return explicit_system->jacobian_states(NV_DATA_S(x), J);
  }

private:
  /**
   * Calculates the root function, using the user-supplied functor
   * for a given value of x.
   * 
   * CHECK -- do the above without converting types?
   */
  inline void f_system(const double x[], double f[]) const {
    const std::vector<double> x_vec(x, x + N_);
    const std::vector<double>& f_vec
      = to_array_1d(f_(to_vector(x_vec), y_, dat_, dat_int_, msgs_));
    std::move(f_vec.begin(), f_vec.end(), f);
  }
  
  /**
   * Calculate the Jacobian of the system function with respect to x.
   */
  inline int jacobian_states(const double x[], SUNMatrix J) const {
    const std::vector<double> x_vec(x, x + N_);
    
    // CHECK - is the best way of constructing the Jacobian
    system_functor<F, double, double, 1>
      system(f_, x_, y_, dat_, dat_int_, msgs_);
    Eigen::VectorXd fx;
    Eigen::MatrixXd Jac;
    jacobian(system, to_vector(x_vec), fx, Jac);

    std::vector<double> jacobian_x = std::vector<double>(N_ * N_);
    Eigen::Map<Eigen::MatrixXd>(&jacobian_x[0], N_, N_) = Jac;
    
    std::move(jacobian_x.begin(), jacobian_x.end(), SM_DATA_D(J));

    return 0;
  }
};

}  // namespace math
}  // namespace stan
#endif
