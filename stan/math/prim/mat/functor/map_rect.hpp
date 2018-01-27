#ifndef STAN_MATH_PRIM_MAT_FUNCTOR_MAP_RECT_HPP
#define STAN_MATH_PRIM_MAT_FUNCTOR_MAP_RECT_HPP

#include <stan/math/prim/mat/fun/typedefs.hpp>

#define STAN_REGISTER_MAP_RECT(CALLID, FUNCTOR)

#ifdef STAN_HAS_MPI
#error "MPI not yet supported"
#else
#include <stan/math/prim/mat/functor/map_rect_serial.hpp>
#endif

#include <vector>

namespace stan {
namespace math {

template <int call_id, typename F, typename T_shared_param,
          typename T_job_param>
Eigen::Matrix<typename stan::return_type<T_shared_param, T_job_param>::type,
              Eigen::Dynamic, 1>
map_rect(const Eigen::Matrix<T_shared_param, Eigen::Dynamic, 1>& shared_params,
         const std::vector<Eigen::Matrix<T_job_param, Eigen::Dynamic, 1>>&
             job_params,
         const std::vector<std::vector<double>>& x_r,
         const std::vector<std::vector<int>>& x_i, std::ostream* msgs = 0) {
#ifdef STAN_HAS_MPI
#error "MPI not yet supported"
#else
  return map_rect_serial<call_id, F, T_shared_param, T_job_param>(
      shared_params, job_params, x_r, x_i, msgs);
#endif
}

}  // namespace math
}  // namespace stan

#endif
