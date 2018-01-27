#ifndef STAN_MATH_PRIM_MAT_FUNCTOR_MAP_RECT_COMBINE_HPP
#define STAN_MATH_PRIM_MAT_FUNCTOR_MAP_RECT_COMBINE_HPP

#include <stan/math/prim/scal/meta/operands_and_partials.hpp>

namespace stan {
namespace math {
namespace internal {

template <typename F, typename T_shared_param, typename T_job_param>
class map_rect_combine {
  typedef operands_and_partials<
      Eigen::Matrix<T_shared_param, Eigen::Dynamic, 1>,
      Eigen::Matrix<T_job_param, Eigen::Dynamic, 1>>
      ops_partials_t;
  std::vector<ops_partials_t> ops_partials_;

  const std::size_t num_shared_operands_;
  const std::size_t num_job_operands_;

 public:
  typedef Eigen::Matrix<
      typename stan::return_type<T_shared_param, T_job_param>::type,
      Eigen::Dynamic, 1>
      result_type;

  map_rect_combine()
      : ops_partials_(), num_shared_operands_(0), num_job_operands_(0) {}
  map_rect_combine(
      const Eigen::Matrix<T_shared_param, Eigen::Dynamic, 1>& shared_params,
      const std::vector<Eigen::Matrix<T_job_param, Eigen::Dynamic, 1>>&
          job_params)
      : ops_partials_(),
        num_shared_operands_(shared_params.rows()),
        num_job_operands_(dims(job_params)[1]) {
    ops_partials_.reserve(job_params.size());
    for (std::size_t i = 0; i < job_params.size(); i++) {
      ops_partials_.push_back(ops_partials_t(shared_params, job_params[i]));
    }
  }

  result_type operator()(const matrix_d& world_result,
                         const std::vector<int>& world_f_out) {
    const std::size_t num_jobs = world_f_out.size();
    const std::size_t offset_job_params
        = is_constant_struct<T_shared_param>::value ? 1
                                                    : 1 + num_shared_operands_;
    const std::size_t size_world_f_out = sum(world_f_out);

    result_type out(size_world_f_out);

    for (std::size_t i = 0, ij = 0; i != num_jobs; ++i) {
      for (std::size_t j = 0; j != world_f_out[i]; ++j, ++ij) {
        if (!is_constant_struct<T_shared_param>::value)
          ops_partials_[i].edge1_.partials_
              = world_result.block(1, ij, num_shared_operands_, 1);

        if (!is_constant_struct<T_job_param>::value)
          ops_partials_[i].edge2_.partials_
              = world_result.block(offset_job_params, ij, num_job_operands_, 1);

        out(ij) = ops_partials_[i].build(world_result(0, ij));
      }
    }

    return out;
  }
};

}  // namespace internal
}  // namespace math
}  // namespace stan

#endif
