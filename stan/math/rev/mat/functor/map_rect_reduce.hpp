#ifndef STAN_MATH_REV_MAT_FUNCTOR_MAP_RECT_REDUCE_HPP
#define STAN_MATH_REV_MAT_FUNCTOR_MAP_RECT_REDUCE_HPP

namespace stan {
  namespace math {

    template <typename F>
    struct map_rect_reduce<F, var, var> {
      matrix_d operator()(const vector_d& shared_params, const vector_d& job_specific_params,
                          const std::vector<double>& x_r, const std::vector<int>& x_i) const {
        const size_type num_shared_params = shared_params.rows();
        const size_type num_job_specific_params = job_specific_params.rows();
        const size_type num_params = num_shared_params + num_job_specific_params;
        matrix_d out(1+num_params,0);

        try {
          start_nested();
          vector_v shared_params_v(num_shared_params);
          vector_v job_specific_params_v(num_job_specific_params);

          for(std::size_t i=0; i < num_shared_params; ++i)
            shared_params_v(i) = shared_params(i);
          for(std::size_t i=0; i < num_job_specific_params; ++i)
            job_specific_params_v(i) = job_specific_params(i);

          std::vector<var> z_vars(num_params);

          for(size_type i = 0; i < num_shared_params; ++i)
            z_vars[i] = shared_params_v(i);
          for(size_type i = 0; i < num_job_specific_params; ++i)
            z_vars[num_shared_params + i] = job_specific_params_v(i);

          const F f;
          vector_v fx_v = f(shared_params_v, job_specific_params_v, x_r, x_i, 0);

          const size_t size_f = fx_v.rows();

          out.resize(Eigen::NoChange, size_f);

          for(size_type i = 0; i < size_f; ++i) {
            out(0,i) = fx_v(i).val();            
            set_zero_all_adjoints_nested();
            fx_v(i).grad();
            for(std::size_t j = 0; j < num_params; ++j)
              out(1+j,i) = z_vars[j].vi_->adj_;
          }
          recover_memory_nested();
        } catch (const std::exception& e) {
          recover_memory_nested();
          throw;
        }
        return out;
      }
    };

    template <typename F>
    struct map_rect_reduce<F, double, var> {
      matrix_d operator()(const vector_d& shared_params, const vector_d& job_specific_params,
                          const std::vector<double>& x_r, const std::vector<int>& x_i) const {
        const size_type num_job_specific_params = job_specific_params.rows();
        matrix_d out(1+num_job_specific_params,0);

        try {
          start_nested();
          vector_v job_specific_params_v(num_job_specific_params);

          for(std::size_t i=0; i < num_job_specific_params; ++i)
            job_specific_params_v(i) = job_specific_params(i);

          const F f;
          vector_v fx_v = f(shared_params, job_specific_params_v, x_r, x_i, 0);

          const size_t size_f = fx_v.rows();

          out.resize(Eigen::NoChange, size_f);

          for(size_type i = 0; i < size_f; ++i) {
            out(0,i) = fx_v(i).val();            
            set_zero_all_adjoints_nested();
            fx_v(i).grad();
            for (std::size_t j = 0; j < num_job_specific_params; ++j)
              out(1+j,i) = job_specific_params_v(j).vi_->adj_;
          }
          recover_memory_nested();
        } catch (const std::exception& e) {
          recover_memory_nested();
          throw;
        }
        return out;
      }
    };

    template <typename F>
    struct map_rect_reduce<F, var, double> {
      matrix_d operator()(const vector_d& shared_params, const vector_d& job_specific_params,
                          const std::vector<double>& x_r, const std::vector<int>& x_i) const {
        const size_type num_shared_params = shared_params.rows();
        matrix_d out(1+num_shared_params,0);

        try {
          start_nested();
          vector_v shared_params_v(num_shared_params);

          for (std::size_t i=0; i < num_shared_params; ++i)
            shared_params_v(i) = shared_params(i);

          const F f;
          vector_v fx_v = f(shared_params_v, job_specific_params, x_r, x_i, 0);

          const size_t size_f = fx_v.rows();

          out.resize(Eigen::NoChange, size_f);

          for (size_type i = 0; i < size_f; ++i) {
            out(0,i) = fx_v(i).val();            
            set_zero_all_adjoints_nested();
            fx_v(i).grad();
            for (std::size_t j = 0; j < num_shared_params; ++j)
              out(1+j,i) = shared_params_v(j).vi_->adj_;
          }
          recover_memory_nested();
        } catch (const std::exception& e) {
          recover_memory_nested();
          throw;
        }
        return out;
      }
    };

  }
}

#ifdef STAN_REGISTER_MPI_MAP_RECT_ALL

#undef STAN_REGISTER_MPI_MAP_RECT_ALL

#define STAN_REGISTER_MPI_MAP_RECT_ALL(CALLID, FUNCTOR)         \
  STAN_REGISTER_MPI_MAP_RECT(CALLID, FUNCTOR, double, double) \
  STAN_REGISTER_MPI_MAP_RECT(CALLID, FUNCTOR, double, var) \
  STAN_REGISTER_MPI_MAP_RECT(CALLID, FUNCTOR, var, double) \
  STAN_REGISTER_MPI_MAP_RECT(CALLID, FUNCTOR, var, var)

#endif

#endif
