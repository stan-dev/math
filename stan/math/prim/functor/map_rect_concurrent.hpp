#ifndef STAN_MATH_PRIM_FUNCTOR_MAP_RECT_CONCURRENT_HPP
#define STAN_MATH_PRIM_FUNCTOR_MAP_RECT_CONCURRENT_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <vector>

namespace stan {
namespace math {
namespace internal {

template <int call_id, typename F, typename T_shared_param,
          typename T_job_param,
          require_eigen_col_vector_t<T_shared_param>* = nullptr>
Eigen::Matrix<return_type_t<T_shared_param, T_job_param>, Eigen::Dynamic, 1>
map_rect_concurrent(
    const T_shared_param& shared_params,
    const std::vector<Eigen::Matrix<T_job_param, Eigen::Dynamic, 1>>&
        job_params,
    const std::vector<std::vector<double>>& x_r,
    const std::vector<std::vector<int>>& x_i, std::ostream* msgs);


template <int call_id, typename F, typename T_shared_param,
          typename T_job_param,
          require_eigen_col_vector_t<T_shared_param>* = nullptr>
Eigen::Matrix<return_type_t<T_shared_param, T_job_param>, Eigen::Dynamic, 1>
map_rect_concurrent(
    const T_shared_param& shared_params,
    const std::vector<Eigen::Matrix<T_job_param, Eigen::Dynamic, 1>>&
        job_params,
    const std::vector<std::vector<double>>& x_r,
    const std::vector<std::vector<int>>& x_i) {
    return map_rect_concurrent<call_id, F>(
        shared_params, job_params, x_r, x_i, nullptr);
    }
}  // namespace internal
}  // namespace math
}  // namespace stan

#endif
