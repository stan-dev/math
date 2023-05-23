#ifndef STAN_MATH_REV_FUNCTOR_IDAS_SERVICE_HPP
#define STAN_MATH_REV_FUNCTOR_IDAS_SERVICE_HPP

#include <stan/math/rev/core/recover_memory.hpp>
#include <stan/math/rev/meta/is_var.hpp>
#include <stan/math/prim/meta/require_helpers.hpp>
#include <stan/math/prim/err/check_flag_sundials.hpp>
#include <idas/idas.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sundials/sundials_context.h>
#include <ostream>
#include <vector>
#include <algorithm>
#include <type_traits>

namespace stan {
namespace math {

/**
 * For each type of Ode(with different rhs functor F and
 * senstivity parameters), we allocate mem and workspace for
 * idas. This service manages the
 * allocation/deallocation, so ODE systems only request
 * service by injection.
 * @tparam ode ode type
 * @tparam lmm_type IDAS solver type (BDF & ADAMS)
 * @tparam butcher_tab AKRODE Butcher table
 */
template <typename dae_type>
struct idas_service {
  sundials::Context sundials_context_;
  int ns;
  N_Vector nv_yy;
  N_Vector nv_yp;
  N_Vector* nv_yys;
  N_Vector* nv_yps;
  void* mem;
  SUNMatrix A;
  SUNLinearSolver LS;

  /**
   * Construct IDAS ODE mem & workspace
   *
   * @param t0 initial time
   * @param dae differential-algebraic system of equations
   */
  idas_service(double t0, dae_type& dae)
      : sundials_context_(),
        ns(dae.ns),
        nv_yy(N_VNew_Serial(dae.N, sundials_context_)),
        nv_yp(N_VNew_Serial(dae.N, sundials_context_)),
        nv_yys(nullptr),
        nv_yps(nullptr),
        mem(IDACreate(sundials_context_)),
        A(SUNDenseMatrix(dae.N, dae.N, sundials_context_)),
        LS(SUNLinSol_Dense(nv_yy, A, sundials_context_)) {
    const int n = dae.N;
    for (auto i = 0; i < n; ++i) {
      NV_Ith_S(nv_yy, i) = dae.dbl_yy[i];
      NV_Ith_S(nv_yp, i) = dae.dbl_yp[i];
    }

    CHECK_IDAS_CALL(IDAInit(mem, dae_type::idas_res, t0, nv_yy, nv_yp));
    CHECK_IDAS_CALL(IDASetUserData(mem, static_cast<void*>(&dae)));
    CHECK_IDAS_CALL(IDASetLinearSolver(mem, LS, A));

    idas_sens_init(nv_yys, nv_yps, ns, n);
  }

  ~idas_service() {
    SUNLinSolFree(LS);
    SUNMatDestroy(A);
    IDAFree(&mem);
    N_VDestroy(nv_yy);
    N_VDestroy(nv_yp);
    if (dae_type::use_fwd_sens) {
      N_VDestroyVectorArray(nv_yys, ns);
      N_VDestroyVectorArray(nv_yps, ns);
    }
  }

  template <typename dae_t = dae_type,
            std::enable_if_t<!dae_t::use_fwd_sens>* = nullptr>
  void idas_sens_init(N_Vector* yys, N_Vector* yps, int ns, int n) {}

  template <typename dae_t = dae_type,
            std::enable_if_t<dae_t::use_fwd_sens>* = nullptr>
  void idas_sens_init(N_Vector*& yys, N_Vector*& yps, int ns, int n) {
    yys = N_VCloneVectorArray(ns, nv_yy);
    yps = N_VCloneVectorArray(ns, nv_yp);
    for (size_t is = 0; is < ns; ++is) {
      N_VConst(RCONST(0.0), yys[is]);
      N_VConst(RCONST(0.0), yps[is]);
    }
    set_init_sens(yys, yps, n);
    CHECK_IDAS_CALL(
        IDASensInit(mem, ns, IDA_STAGGERED, dae_type::idas_sens_res, yys, yps));
  }

  template <typename dae_t = dae_type,
            std::enable_if_t<dae_t::is_var_yy0 && dae_t::is_var_yp0>* = nullptr>
  void set_init_sens(N_Vector*& yys, N_Vector*& yps, int n) {
    for (size_t i = 0; i < n; ++i) {
      NV_Ith_S(yys[i], i) = 1.0;
    }
    for (size_t i = 0; i < n; ++i) {
      NV_Ith_S(yps[i + n], i) = 1.0;
    }
  }

  template <
      typename dae_t = dae_type,
      std::enable_if_t<dae_t::is_var_yy0 && (!dae_t::is_var_yp0)>* = nullptr>
  void set_init_sens(N_Vector*& yys, N_Vector*& yps, int n) {
    for (size_t i = 0; i < n; ++i) {
      NV_Ith_S(yys[i], i) = 1.0;
    }
  }

  template <
      typename dae_t = dae_type,
      std::enable_if_t<(!dae_t::is_var_yy0) && dae_t::is_var_yp0>* = nullptr>
  void set_init_sens(N_Vector*& yys, N_Vector*& yps, int n) {
    for (size_t i = 0; i < n; ++i) {
      NV_Ith_S(yps[i], i) = 1.0;
    }
  }

  template <
      typename dae_t = dae_type,
      std::enable_if_t<(!dae_t::is_var_yy0) && (!dae_t::is_var_yp0)>* = nullptr>
  void set_init_sens(N_Vector*& yys, N_Vector*& yps, int n) {}
};
}  // namespace math
}  // namespace stan

#endif
