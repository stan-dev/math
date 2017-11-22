#include <stan/math.hpp>
//#include <gtest/gtest.h>

#include <test/unit/math/prim/mat/functor/hard_work.hpp>

#include <iostream>

#ifdef STAN_HAS_MPI

typedef stan::math::map_rect_reduce<hard_work, double, double> hard_work_reducer_dd;
typedef stan::math::map_rect_combine<hard_work, double, double> hard_work_combiner_dd;
typedef stan::math::mpi_parallel_call<hard_work_reducer_dd,hard_work_combiner_dd> hard_work_parallel_call;
BOOST_CLASS_EXPORT(stan::math::mpi_distributed_apply<hard_work_parallel_call>);
BOOST_CLASS_TRACKING(stan::math::mpi_distributed_apply<hard_work_parallel_call>,track_never);

#endif

int main(int argc, const char* argv[]) {
#ifdef STAN_HAS_MPI
  boost::mpi::environment env;
  
  // on non-root processes this makes the workers listen to commands
  // send from the root
  stan::math::mpi_cluster cluster;

#endif

  Eigen::VectorXd shared_params_d(2);
  shared_params_d << 2, 0;
  std::vector<Eigen::VectorXd> job_params_d;

  const std::size_t N = 10;

  for(std::size_t n = 0; n != N; ++n) {
    Eigen::VectorXd job_d(2);
    job_d << 0, n * n;
    job_params_d.push_back(job_d);
  }
  
  std::vector<std::vector<double> > x_r(N, std::vector<double>(1,1.0));
  std::vector<std::vector<int> > x_i(N, std::vector<int>(0));

  Eigen::VectorXd result = stan::math::map_rect<hard_work>(shared_params_d, job_params_d, x_r, x_i, 0);

  std::cout << "Root process ends." << std::endl;
  
  return(0);
}


