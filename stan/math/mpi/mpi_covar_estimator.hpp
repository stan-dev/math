#ifndef STAN_MATH_MPI_COVAR_ESTIMATOR_HPP
#define STAN_MATH_MPI_COVAR_ESTIMATOR_HPP

#ifdef STAN_LANG_MPI

#include <stan/math/mpi/envionment.hpp>
#include <stan/math/prim/mat/fun/welford_covar_estimator.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/sum.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/statistics/count.hpp>
#include <boost/mpi.hpp>
#include <algorithm>
#include <vector>

namespace stan {
namespace math {
namespace mpi {

  /*
   * brute force by using welford. Assuming all chains have
   * same nb. of draws.
   */
class mpi_covar_estimator {
  boost::accumulators::accumulator_set<int, boost::accumulators::features<boost::accumulators::tag::count>> count_acc_;
 public:
  size_t num_iterations;

protected:
  Eigen::MatrixXd draws_;

public:
  explicit mpi_covar_estimator(size_t num_params,
                               size_t num_iterations0)
    : num_iterations(num_iterations0),
      draws_(num_params, num_iterations)
  {}

  void restart() { count_acc_ = {}; }

  void restart(size_t num_params) { count_acc_ = {}; }

  inline void add_sample(const Eigen::VectorXd& q) {
    int i = boost::accumulators::count(count_acc_);
    draws_.col(i) = q;
    count_acc_(0);
  }

  inline int num_samples(const Communicator& comm) {
    int n = boost::accumulators::count(count_acc_);
    int n_all = n * comm.size();
    return n_all;
  }

  inline void sample_mean(Eigen::VectorXd& var) {
    // todo
  }

  inline void sample_covariance(Eigen::MatrixXd& covar, const Communicator& comm) {
    sample_covariance(covar, 0, boost::accumulators::count(count_acc_), comm);
  }

  inline void sample_covariance(Eigen::MatrixXd& covar,
                                int col_begin, int n_samples,
                                const Communicator& comm) {                                
    int num_params = draws_.rows();
    Eigen::MatrixXd all_draws(num_params, n_samples * comm.size());
    MPI_Allgather(draws_.data() + col_begin * num_params,
                  n_samples * num_params, MPI_DOUBLE, all_draws.data(),
                  n_samples * num_params, MPI_DOUBLE, comm.comm());
    welford_covar_estimator est(num_params);
    for (int j = 0; j < all_draws.cols(); ++j) {    
      est.add_sample(all_draws.col(j));
    }
    est.sample_covariance(covar);
  }

  inline size_t num_params() {
    return draws_.rows();
  }
};

}
}  // namespace math
}  // namespace stan

#endif

#endif
