#ifndef STAN_MATH_REV_FUNCTOR_IDAS_INTEGRATOR_HPP
#define STAN_MATH_REV_FUNCTOR_IDAS_INTEGRATOR_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/functor/idas_service.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/err/check_flag_sundials.hpp>
#include <idas/idas.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <nvector/nvector_serial.h>
#include <ostream>
#include <vector>
#include <algorithm>

namespace stan {
namespace math {

/**
 * IDAS DAE integrator.
 */
class idas_integrator {
  sundials::Context sundials_context_;

  const double rtol_;
  const double atol_;
  const int64_t max_num_steps_;

 public:
  /**
   * constructor
   * @param[in] rtol relative tolerance
   * @param[in] atol absolute tolerance
   * @param[in] max_num_steps max nb. of times steps
   */
  idas_integrator(const double rtol, const double atol,
                  const int64_t max_num_steps)
      : rtol_(rtol), atol_(atol), max_num_steps_(max_num_steps) {}

  /**
   * Return the solutions for the specified DAE
   * given the specified initial state,
   * initial times, times of desired solution, and parameters and
   * data, writing error and warning messages to the specified
   * stream contained in the DAE system.
   *
   * @tparam DAE type of DAE system
   * @param[in] dae DAE system
   * @param[in] t0 initial time.
   * @param[in] ts times of the desired solutions, in strictly
   * increasing order, all greater than the initial time.
   * @return a vector of states, each state being a vector of the
   * same size as the state variable, corresponding to a time in ts.
   */
  template <typename dae_type>
  typename dae_type::return_t operator()(const char* func, dae_type& dae,
                                         double t0,
                                         const std::vector<double>& ts) {
    idas_service<dae_type> serv(t0, dae);

    void* mem = serv.mem;
    N_Vector& yy = serv.nv_yy;
    N_Vector& yp = serv.nv_yp;
    N_Vector* yys = serv.nv_yys;
    N_Vector* yps = serv.nv_yps;
    const size_t n = dae.N;

    CHECK_IDAS_CALL(IDASStolerances(mem, rtol_, atol_));
    CHECK_IDAS_CALL(IDASetMaxNumSteps(mem, max_num_steps_));

    if (dae_type::use_fwd_sens) {
      CHECK_IDAS_CALL(IDASensEEtolerances(mem));
      CHECK_IDAS_CALL(IDAGetSensConsistentIC(mem, yys, yps));
    }

    size_t nt = ts.size();
    typename dae_type::return_t res_yy(
        nt, Eigen::Matrix<typename dae_type::scalar_t, -1, 1>::Zero(n));
    double t1 = t0;
    for (size_t i = 0; i < nt; ++i) {
      CHECK_IDAS_CALL(IDASolve(mem, ts[i], &t1, yy, yp, IDA_NORMAL));
      if (dae_type::use_fwd_sens) {
        CHECK_IDAS_CALL(IDAGetSens(mem, &t1, yys));
      }
      collect(yy, yys, dae, res_yy[i]);
    }

    return res_yy;
  }

  /**
   *
   * Solve DAE system, no sensitivity
   *
   * @tparam F DAE functor type
   * @param[out] dae DAE system
   * @param[in] t0 initial time
   * @param[in] ts times of the desired solutions
   * @param[out] res_yy DAE solutions
   */
  template <typename dae_type>
  void collect(N_Vector const& yy, N_Vector const* yys, dae_type& dae,
               Eigen::Matrix<double, -1, 1>& y) {
    for (int i = 0; i < dae.N; ++i) {
      y.coeffRef(i) = NV_Ith_S(yy, i);
    }
  }

  template <typename dae_type>
  void collect(N_Vector const& yy, N_Vector const* yys, dae_type& dae,
               Eigen::Matrix<stan::math::var, -1, 1>& y) {
    std::vector<double> g(dae.ns);
    for (size_t i = 0; i < dae.N; ++i) {
      for (size_t j = 0; j < dae.ns; ++j) {
        g[j] = NV_Ith_S(yys[j], i);
      }
      y[i] = precomputed_gradients(NV_Ith_S(yy, i), dae.all_vars, g);
    }
  }
};

}  //   namespace math
}  //   namespace stan

#endif
