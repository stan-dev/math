#ifndef STAN_MATH_MPI_VAR_ESTIMATOR_HPP
#define STAN_MATH_MPI_VAR_ESTIMATOR_HPP

#ifdef STAN_LANG_MPI

#include <stan/math/mpi/envionment.hpp>
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

class mpi_var_estimator {
  using acc_param = boost::accumulators::accumulator_set<double,
                                                         boost::accumulators::stats<
                                                         boost::accumulators::tag::sum,
                                                         boost::accumulators::tag::mean,
                                                         boost::accumulators::tag::variance>>;
  using acc_count = boost::accumulators::accumulator_set<int,
                                                         boost::accumulators::features<boost::accumulators::tag::count>>;

 public:
  explicit mpi_var_estimator(size_t num_params)
    : acc_(num_params)
  {}

  void restart() {
    for (size_t i = 0; i < acc_.size(); ++i) {
      acc_[i] = {};
    }
    count_acc_ = {};
  }

  void restart(size_t num_params) {
    acc_.resize(num_params);
    restart();
  }

  inline void add_sample(const Eigen::VectorXd& q) {
    for (size_t i = 0; i < acc_.size(); ++i) {
      acc_[i](q(i));
    }
    count_acc_(0);
  }

  inline std::vector<int> num_samples(const Communicator& comm) {
    int n = boost::accumulators::count(count_acc_);
    int n_all;
    MPI_Allreduce(&n, &n_all, 1, MPI_INT, MPI_SUM, comm.comm());    
    return {n, n_all};
  }

  inline void sample_mean(Eigen::VectorXd& var, int n_all, const Communicator& comm) {
    int num_params = acc_.size();
    Eigen::ArrayXd work(num_params);
    var.resize(num_params);
    for (int i = 0; i < num_params; ++i) {    
      work(i) = boost::accumulators::sum(acc_[i]);
    }
    MPI_Allreduce(work.data(), var.data(), num_params, MPI_DOUBLE, MPI_SUM, comm.comm());
    var /= n_all;
  }

  inline void sample_mean(Eigen::VectorXd& var, const Communicator& comm) {
    int n_all = num_samples(comm)[1];
    sample_mean(var, n_all, comm);
  }

  inline int sample_variance(Eigen::VectorXd& var, const Communicator& comm) {
    std::vector<int> n_samples(num_samples(comm));
    int n = n_samples[0];
    int n_all = n_samples[1];
    int m = comm.size();
    int num_params = acc_.size();
    Eigen::VectorXd work(num_params);
    sample_mean(work, n_all, comm);
    for (int i = 0; i < num_params; ++i) {    
      double d = boost::accumulators::mean(acc_[i]) - work(i);
      work(i) = n * boost::accumulators::variance(acc_[i]) + n * d * d;
    }
    var.resize(num_params);
    MPI_Allreduce(work.data(), var.data(), num_params, MPI_DOUBLE, MPI_SUM, comm.comm());
    var /= (n_all - 1.0);
    return n_all;
  }

  inline size_t num_params() {
    return acc_.size();
  }

 protected:
  std::vector<acc_param> acc_;
  acc_count count_acc_;
};

  }
}  // namespace math
}  // namespace stan

#endif

#endif
