#ifndef STAN_MATH_REV_MAT_FUNCTOR_KINSOL_DATA_HPP
#define STAN_MATH_REV_MAT_FUNCTOR_KINSOL_DATA_HPP

#include <stan/math/prim/mat/fun/to_array_1d.hpp>
#include <stan/math/prim/mat/fun/to_vector.hpp>
#include <stan/math/rev/mat/functor/algebra_system.hpp>
#include <stan/math/rev/mat/functor/jacobian.hpp>

#include <kinsol/kinsol.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <nvector/nvector_serial.h>

#include <vector>

namespace stan {
namespace math {

/**
 * Default Jacobian builder using revser-mode autodiff.
 */
struct kinsol_J_f {
  template <typename F>
  inline int operator()(const F& f, const Eigen::VectorXd& x,
                        const Eigen::VectorXd& y,
                        const std::vector<double>& dat,
                        const std::vector<int>& dat_int, std::ostream* msgs,
                        const double x_sun[], SUNMatrix J) const {
    size_t N = x.size();
    const std::vector<double> x_vec(x_sun, x_sun + N);
    system_functor<F, double, double, 1> system(f, x, y, dat, dat_int, msgs);
    Eigen::VectorXd fx;
    Eigen::MatrixXd Jac;
    jacobian(system, to_vector(x_vec), fx, Jac);

    for (int i = 0; i < Jac.rows(); i++)
      for (int j = 0; j < Jac.cols(); j++)
        SM_ELEMENT_D(J, i, j) = Jac(i, j);

    return 0;
  }
};

/**
 * KINSOL algebraic system data holder.
 * Based on cvodes_ode_data.
 *
 * @tparam F1 functor type for system function.
 * @tparam F2 functor type for jacobian function. Default is 0.
 *         If 0, use rev mode autodiff to compute the Jacobian.
 */
template <typename F1, typename F2>
class kinsol_system_data {
  const F1& f_;
  const F2& J_f_;
  const Eigen::VectorXd& x_;
  const Eigen::VectorXd& y_;
  const size_t N_;
  const std::vector<double>& dat_;
  const std::vector<int>& dat_int_;
  std::ostream* msgs_;

  typedef kinsol_system_data<F1, F2> system_data;

 public:
  N_Vector nv_x_;
  SUNMatrix J_;
  SUNLinearSolver LS_;
  void* kinsol_memory_;

  /* Constructor */
  kinsol_system_data(const F1& f, const F2& J_f, const Eigen::VectorXd& x,
                     const Eigen::VectorXd& y, const std::vector<double>& dat,
                     const std::vector<int>& dat_int, std::ostream* msgs)
      : f_(f),
        J_f_(J_f),
        x_(x),
        y_(y),
        dat_(dat),
        dat_int_(dat_int),
        msgs_(msgs),
        N_(x.size()),
        nv_x_(N_VMake_Serial(N_, &to_array_1d(x_)[0])),
        J_(SUNDenseMatrix(N_, N_)),
        LS_(SUNLinSol_Dense(nv_x_, J_)),
        kinsol_memory_(KINCreate()) {}

  ~kinsol_system_data() {
    N_VDestroy_Serial(nv_x_);
    SUNLinSolFree(LS_);
    SUNMatDestroy(J_);
    KINFree(&kinsol_memory_);
  }

  /* Implements the user-defined function passed to KINSOL. */
  static int kinsol_f_system(N_Vector x, N_Vector f, void* user_data) {
    const system_data* explicit_system
        = static_cast<const system_data*>(user_data);

    Eigen::VectorXd x_eigen(
        Eigen::Map<Eigen::VectorXd>(NV_DATA_S(x), explicit_system->N_));

    Eigen::Map<Eigen::VectorXd>(N_VGetArrayPointer(f), explicit_system->N_)
        = explicit_system->f_(x_eigen, explicit_system->y_,
                              explicit_system->dat_, explicit_system->dat_int_,
                              explicit_system->msgs_);

    return 0;
  }

  /**
   * Implements the function of type CVDlsJacFn which is the user-defined
   * callbacks for KINSOL to calculate the jacobian of the system.
   * The Jacobian is stored in column major format.
   *
   * REMARK - tmp1 and tmp2 are pointers to memory allocated for variables
   * of type N_Vector which can be used by KINJacFN (the function which
   * computes the Jacobian) as temporary storage or work space.
   * See
   * https://computation.llnl.gov/sites/default/files/public/kin_guide-dev.pdf,
   * page 55.
   */
  static int kinsol_jacobian(N_Vector x, N_Vector f, SUNMatrix J,
                             void* user_data, N_Vector tmp1, N_Vector tmp2) {
    const system_data* explicit_system
        = static_cast<const system_data*>(user_data);
    return explicit_system->J_f_(explicit_system->f_, explicit_system->x_,
                                 explicit_system->y_, explicit_system->dat_,
                                 explicit_system->dat_int_,
                                 explicit_system->msgs_, NV_DATA_S(x), J);
  }
};

}  // namespace math
}  // namespace stan
#endif
