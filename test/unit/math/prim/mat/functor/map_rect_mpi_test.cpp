#include <stan/math.hpp>
#include <gtest/gtest.h>

#include <test/unit/math/prim/mat/functor/mpi_test_env.hpp>

#include <test/unit/math/prim/mat/functor/hard_work.hpp>
#include <test/unit/math/prim/mat/functor/faulty_functor.hpp>

#include <iostream>

STAN_REGISTER_MPI_MAP_RECT(0, hard_work, double, double)
STAN_REGISTER_MPI_MAP_RECT(1, faulty_functor, double, double)
STAN_REGISTER_MPI_MAP_RECT(2, faulty_functor, double, double)

struct MpiJob : public ::testing::Test {
  Eigen::VectorXd shared_params_d;
  std::vector<Eigen::VectorXd> job_params_d;
  const std::size_t N = 10;
  std::vector<std::vector<double> > x_r = std::vector<std::vector<double>>(N, std::vector<double>(1,1.0));
  std::vector<std::vector<int> > x_i = std::vector<std::vector<int>>(N, std::vector<int>(1,0));

   virtual void SetUp() {
     shared_params_d.resize(2);
     shared_params_d << 2, 0;
     
     for(std::size_t n = 0; n != N; ++n) {
       x_i[n][0] = n;
       Eigen::VectorXd job_d(2);
       job_d << n+1.0, n * n;
       job_params_d.push_back(job_d);
     }
   }
  
};

TEST_F(MpiJob, hard_work_dd) {
  if(rank != 0) return;
  
  Eigen::VectorXd result_mpi = stan::math::map_rect_mpi<0,hard_work>(shared_params_d, job_params_d, x_r, x_i);
  Eigen::VectorXd result_serial = stan::math::map_rect_serial<0,hard_work>(shared_params_d, job_params_d, x_r, x_i);

  EXPECT_EQ(result_mpi.rows(), result_serial.rows() );

  for(int i = 0; i < result_mpi.rows(); ++i) {
    EXPECT_DOUBLE_EQ(result_mpi(i), result_serial(i));
  }

} 

TEST_F(MpiJob, always_faulty_functor) {
  if(rank != 0) return;

  Eigen::VectorXd result;

  EXPECT_NO_THROW((result = stan::math::map_rect<1,faulty_functor>(shared_params_d, job_params_d, x_r, x_i)));

  // faulty functor throws on theta(0) being -1.0
  // throwing during the first evaluation is quite severe and will
  // lead to a respective runtime error
  job_params_d[0](0) = -1;

  // upon the second evaluation throwing is handled internally different
  EXPECT_ANY_THROW((result = stan::math::map_rect<1,faulty_functor>(shared_params_d, job_params_d, x_r, x_i)));

  // thorwing on the very first evaluation
  EXPECT_ANY_THROW((result = stan::math::map_rect<2,faulty_functor>(shared_params_d,job_params_d, x_r, x_i)));
}
