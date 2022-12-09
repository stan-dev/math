#ifndef STAN_MATH_REV_FUNCTOR_KINSOL_DATA_HPP
#define STAN_MATH_REV_FUNCTOR_KINSOL_DATA_HPP

#include <stan/math/rev/functor/algebra_system.hpp>
#include <stan/math/rev/functor/jacobian.hpp>
#include <stan/math/prim/fun/to_array_1d.hpp>
#include <stan/math/prim/fun/to_vector.hpp>
#include <stan/math/prim/functor/apply.hpp>
#include <sundials/sundials_context.h>
#include <kinsol/kinsol.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <nvector/nvector_serial.h>
#include <vector>
#include <tuple>

namespace stan {
namespace math {

/**
 * KINSOL algebraic system data holder.
 * Based on cvodes_ode_data.
 *
 * @tparam F1 functor type for system function.
 * @tparam F2 functor type for jacobian function. Default is 0.
 *         If 0, use rev mode autodiff to compute the Jacobian.
 */
template <typename F1, typename... Args>
class kinsol_system_data {
  const F1& f_;
  const Eigen::VectorXd& x_;
  const size_t N_;
  std::ostream* const msgs_;
  const std::tuple<const Args&...> args_tuple_;

  typedef kinsol_system_data<F1, Args...> system_data;

 public:
  sundials::Context sundials_context_;
  N_Vector nv_x_;
  SUNMatrix J_;
  SUNLinearSolver LS_;
  void* kinsol_memory_;

  /* Constructor */
  kinsol_system_data(const F1& f, const Eigen::VectorXd& x,
                     std::ostream* const msgs, const Args&... args)
      : f_(f),
        x_(x),
        N_(x.size()),
        msgs_(msgs),
        args_tuple_(args...),
        sundials_context_(),
        nv_x_(N_VMake_Serial(N_, &to_array_1d(x_)[0], sundials_context_)),
        J_(SUNDenseMatrix(N_, N_, sundials_context_)),
        LS_(SUNLinSol_Dense(nv_x_, J_, sundials_context_)),
        kinsol_memory_(KINCreate(sundials_context_)) {}

  ~kinsol_system_data() {
    N_VDestroy_Serial(nv_x_);
    SUNLinSolFree(LS_);
    SUNMatDestroy(J_);
    KINFree(&kinsol_memory_);
  }

  /* Implements the user-defined function passed to KINSOL. */
  static int kinsol_f_system(const N_Vector x, const N_Vector f_eval,
                             void* const user_data) {
    const system_data* explicit_system
        = static_cast<const system_data*>(user_data);

    Eigen::VectorXd x_eigen(
        Eigen::Map<Eigen::VectorXd>(NV_DATA_S(x), explicit_system->N_));

    Eigen::Map<Eigen::VectorXd> f_eval_map(N_VGetArrayPointer(f_eval),
                                           explicit_system->N_);
    auto result = math::apply(
        [&](const auto&... args) {
          return explicit_system->f_(x_eigen, explicit_system->msgs_, args...);
        },
        explicit_system->args_tuple_);
    check_matching_sizes("", "the algebraic system's output", result,
                         "the vector of unknowns, x,", f_eval_map);
    f_eval_map = result;
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
  static int kinsol_jacobian(const N_Vector x, const N_Vector f_eval,
                             const SUNMatrix J, void* const user_data,
                             const N_Vector tmp1, const N_Vector tmp2) {
    const system_data* explicit_system
        = static_cast<const system_data*>(user_data);

    Eigen::VectorXd x_eigen(
        Eigen::Map<Eigen::VectorXd>(NV_DATA_S(x), explicit_system->N_));
    Eigen::Map<Eigen::VectorXd> f_eval_map(N_VGetArrayPointer(f_eval),
                                           explicit_system->N_);

    auto f_wrt_x = [&](const auto& x) {
      return math::apply(
          [&](const auto&... args) {
            return explicit_system->f_(x, explicit_system->msgs_, args...);
          },
          explicit_system->args_tuple_);
    };

    Eigen::MatrixXd Jf_x;
    Eigen::VectorXd f_x;

    jacobian(f_wrt_x, x_eigen, f_x, Jf_x);

    f_eval_map = f_x;

    for (int j = 0; j < Jf_x.cols(); j++)
      for (int i = 0; i < Jf_x.rows(); i++)
        SM_ELEMENT_D(J, i, j) = Jf_x(i, j);

    return 0;
  }
};

}  // namespace math
}  // namespace stan
#endif
