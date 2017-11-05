#pragma once

#ifdef STAN_HAS_MPI
#include <stan/math/prim/mat/functor/map_rect_mpi.hpp>
#endif

namespace stan {
  namespace math {

    template <typename F, typename T_shared_param, typename T_param>
    Eigen::Matrix<typename stan::return_type<T_shared_param, T_param>::type, Eigen::Dynamic, 1>
    map_rect(const Eigen::Matrix<T_shared_param, Eigen::Dynamic, 1>& eta,
             const Eigen::Matrix<T_param, Eigen::Dynamic, Eigen::Dynamic>& theta,
             const Eigen::MatrixXd& x_r,
             const std::vector<std::vector<int> >& x_i,
             const int callsite_id) {
#ifdef STAN_HAS_MPI
      return(map_rect_mpi<F>(eta, theta, x_r, x_i, callsite_id));
#else
      typedef typename stan::return_type<T_shared_param, T_param>::type T_out;
      Eigen::Matrix<T_out, Eigen::Dynamic, 1> res;
      const std::size_t num_jobs = theta.cols();

      // TODO: transpose x_i

      for(std::size_t i = 0; i != num_jobs; ++i) {
        const Eigen::Matrix<T_out, Eigen::Dynamic, 1> f = F::apply(eta, theta.col(i), x_r.col(i), x_i[i]);
        const std::size_t rows = res.rows();
        const std::size_t f_size = f.rows();
        res.conservativeResize(rows + f_size, 1);
        res.bottomRows(f_size) = f;
      }
      return(res);
#endif
    }
  }
}
