#pragma once

#ifdef STAN_HAS_MPI
#include <stan/math/prim/mat/functor/map_rect_mpi.hpp>
#endif

namespace stan {
  namespace math {

    /*
    template <typename T>
    std::vector<std::vector<T> > transpose(const std::vector<std::vector<T> >& in) {
      const std::size_t num_cols = in.size();
      if(num_cols == 0) return(in);
      const std::size_t num_rows = in[0].size();
      std::vector<std::vector<T> > out(num_rows, std::vector<T>(num_cols));

      for(std::size_t i=0; i < num_cols; ++i)
        for(std::size_t j=0; j < num_rows; ++j)
          out[j][i] = in[i][j];
      
      return(out);
    }
    */

    template <typename F, typename T_shared_param, typename T_job_param>
    Eigen::Matrix<typename stan::return_type<T_shared_param, T_job_param>::type, Eigen::Dynamic, 1>
    map_rect(const Eigen::Matrix<T_shared_param, Eigen::Dynamic, 1>& shared_params,
             const std::vector<Eigen::Matrix<T_job_param, Eigen::Dynamic, 1> >& job_params,
             const std::vector<std::vector<double> >& x_r,
             const std::vector<std::vector<int> >& x_i,
             const int callsite_id) {
#ifdef STAN_HAS_MPI
      return(map_rect_mpi<F>(eta, theta, x_r, x_i, callsite_id));
#else
      typedef typename stan::return_type<T_shared_param, T_job_param>::type result_type;
      Eigen::Matrix<result_type, Eigen::Dynamic, 1> out(job_params.size());
      const std::size_t num_jobs = job_params.size();
      int out_size = 0;

      for(std::size_t i = 0; i != num_jobs; ++i) {
        const Eigen::Matrix<result_type, Eigen::Dynamic, 1> f = F::apply(shared_params, job_params[i], x_r[i], x_i[i]);
        const int f_size = f.rows();
        out_size += f_size;
        if(out.rows() < out_size)
          out.conservativeResize(2*out_size, 1);
        out.segment(out_size - f_size, f_size) = f;
      }
      out.conservativeResize(out_size, 1);
      return(out);
#endif
    }
  }
}

