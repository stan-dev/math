#ifndef STAN_MATH_REV_FUNCTOR_MAP_RECT_REDUCE_HPP
#define STAN_MATH_REV_FUNCTOR_MAP_RECT_REDUCE_HPP

#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/core/typedefs.hpp>
#include <stan/math/rev/fun/to_var.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <stan/math/prim/functor/map_rect_reduce.hpp>
#include <vector>

namespace stan {
namespace math {
namespace internal {

template <typename F>
struct map_rect_reduce<F, var, var> {
  matrix_d operator()(const vector_d& shared_params,
                      const vector_d& job_specific_params,
                      const std::vector<double>& x_r,
                      const std::vector<int>& x_i,
                      std::ostream* msgs = nullptr) const {
    const size_type num_shared_params = shared_params.rows();
    const size_type num_job_specific_params = job_specific_params.rows();
    matrix_d out(1 + num_shared_params + num_job_specific_params, 0);

    // Run nested autodiff in this scope
    nested_rev_autodiff nested;

    vector_v shared_params_v = to_var(shared_params);
    vector_v job_specific_params_v = to_var(job_specific_params);

    vector_v fx_v = F()(shared_params_v, job_specific_params_v, x_r, x_i, msgs);

    const size_type size_f = fx_v.rows();

    out.resize(Eigen::NoChange, size_f);

    for (size_type i = 0; i < size_f; ++i) {
      out(0, i) = fx_v(i).val();
      nested.set_zero_all_adjoints();
      fx_v(i).grad();
      for (size_type j = 0; j < num_shared_params; ++j) {
        out(1 + j, i) = shared_params_v(j).vi_->adj_;
      }
      for (size_type j = 0; j < num_job_specific_params; ++j) {
        out(1 + num_shared_params + j, i) = job_specific_params_v(j).vi_->adj_;
      }
    }
    return out;
  }
};

template <typename F>
struct map_rect_reduce<F, double, var> {
  matrix_d operator()(const vector_d& shared_params,
                      const vector_d& job_specific_params,
                      const std::vector<double>& x_r,
                      const std::vector<int>& x_i,
                      std::ostream* msgs = nullptr) const {
    const size_type num_job_specific_params = job_specific_params.rows();
    matrix_d out(1 + num_job_specific_params, 0);

    // Run nested autodiff in this scope
    nested_rev_autodiff nested;

    vector_v job_specific_params_v = to_var(job_specific_params);

    vector_v fx_v = F()(shared_params, job_specific_params_v, x_r, x_i, msgs);

    const size_type size_f = fx_v.rows();

    out.resize(Eigen::NoChange, size_f);

    for (size_type i = 0; i < size_f; ++i) {
      out(0, i) = fx_v(i).val();
      nested.set_zero_all_adjoints();
      fx_v(i).grad();
      for (size_type j = 0; j < num_job_specific_params; ++j) {
        out(1 + j, i) = job_specific_params_v(j).vi_->adj_;
      }
    }
    return out;
  }
};

template <typename F>
struct map_rect_reduce<F, var, double> {
  matrix_d operator()(const vector_d& shared_params,
                      const vector_d& job_specific_params,
                      const std::vector<double>& x_r,
                      const std::vector<int>& x_i,
                      std::ostream* msgs = nullptr) const {
    const size_type num_shared_params = shared_params.rows();
    matrix_d out(1 + num_shared_params, 0);

    // Run nested autodiff in this scope
    nested_rev_autodiff nested;

    vector_v shared_params_v = to_var(shared_params);

    vector_v fx_v = F()(shared_params_v, job_specific_params, x_r, x_i, msgs);

    const size_type size_f = fx_v.rows();

    out.resize(Eigen::NoChange, size_f);

    for (size_type i = 0; i < size_f; ++i) {
      out(0, i) = fx_v(i).val();
      nested.set_zero_all_adjoints();
      fx_v(i).grad();
      for (size_type j = 0; j < num_shared_params; ++j) {
        out(1 + j, i) = shared_params_v(j).vi_->adj_;
      }
    }

    return out;
  }
};

}  // namespace internal
}  // namespace math
}  // namespace stan

#ifdef STAN_REGISTER_MPI_MAP_RECT_ALL

#undef STAN_REGISTER_MPI_MAP_RECT_ALL

#define STAN_REGISTER_MPI_MAP_RECT_ALL(CALLID, FUNCTOR)       \
  STAN_REGISTER_MPI_MAP_RECT(CALLID, FUNCTOR, double, double) \
  STAN_REGISTER_MPI_MAP_RECT(CALLID, FUNCTOR, double, var)    \
  STAN_REGISTER_MPI_MAP_RECT(CALLID, FUNCTOR, var, double)    \
  STAN_REGISTER_MPI_MAP_RECT(CALLID, FUNCTOR, var, var)

#endif

#endif
