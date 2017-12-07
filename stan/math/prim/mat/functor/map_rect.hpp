#pragma once

#ifdef STAN_HAS_MPI
#include <stan/math/prim/mat/functor/map_rect_mpi.hpp>
#endif

namespace stan {
  namespace math {

    template <int call_id, typename F, typename T_shared_param, typename T_job_param>
    Eigen::Matrix<typename stan::return_type<T_shared_param, T_job_param>::type, Eigen::Dynamic, 1>
    map_rect_serial(const Eigen::Matrix<T_shared_param, Eigen::Dynamic, 1>& shared_params,
                    const std::vector<Eigen::Matrix<T_job_param, Eigen::Dynamic, 1> >& job_params,
                    const std::vector<std::vector<double> >& x_r,
                    const std::vector<std::vector<int> >& x_i,
                    std::ostream* msgs = 0) {
      typedef typename stan::return_type<T_shared_param, T_job_param>::type result_type;
      Eigen::Matrix<result_type, Eigen::Dynamic, 1> out;
      const std::size_t num_jobs = job_params.size();
      int out_size = 0;
      const F f;
      
      for(std::size_t i = 0; i != num_jobs; ++i) {
        const Eigen::Matrix<result_type, Eigen::Dynamic, 1> fx = f(shared_params, job_params[i], x_r[i], x_i[i], msgs);
        const int fx_size = fx.rows();
        out_size += fx_size;
        if(i == 0) out.resize(num_jobs * fx_size);
        if(out.rows() < out_size)
          out.conservativeResize(2*out_size);
        out.segment(out_size - fx_size, fx_size) = fx;
      }
      out.conservativeResize(out_size);
      return(out);
    }

    template <int call_id, typename F, typename T_shared_param, typename T_job_param>
    Eigen::Matrix<typename stan::return_type<T_shared_param, T_job_param>::type, Eigen::Dynamic, 1>
    map_rect(const Eigen::Matrix<T_shared_param, Eigen::Dynamic, 1>& shared_params,
             const std::vector<Eigen::Matrix<T_job_param, Eigen::Dynamic, 1> >& job_params,
             const std::vector<std::vector<double> >& x_r,
             const std::vector<std::vector<int> >& x_i,
             std::ostream* msgs = 0) {
#ifdef STAN_HAS_MPI
      return(map_rect_mpi<call_id,F,T_shared_param,T_job_param>(shared_params, job_params, x_r, x_i, msgs));
#else
      return(map_rect_serial<call_id,F,T_shared_param,T_job_param>(shared_params, job_params, x_r, x_i, msgs));
#endif
    }

    
  }
}

