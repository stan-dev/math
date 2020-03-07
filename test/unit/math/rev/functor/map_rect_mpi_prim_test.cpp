// these tests can only be compiled and executed with availability of
// MPI
#ifdef STAN_MPI

#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>

#include <test/unit/math/prim/functor/hard_work.hpp>
#include <test/unit/math/prim/functor/faulty_functor.hpp>

#include <iostream>
#include <vector>

STAN_REGISTER_MAP_RECT(0, hard_work)
STAN_REGISTER_MAP_RECT(1, hard_work)
STAN_REGISTER_MAP_RECT(2, faulty_functor)
STAN_REGISTER_MAP_RECT(3, faulty_functor)

struct MpiJob : public ::testing::Test {
  Eigen::VectorXd shared_params_d;
  std::vector<Eigen::VectorXd> job_params_d;
  const std::size_t N = 10;
  std::vector<std::vector<double>> x_r
      = std::vector<std::vector<double>>(N, std::vector<double>(1, 1.0));
  std::vector<std::vector<int>> x_i
      = std::vector<std::vector<int>>(N, std::vector<int>(1, 0));

  virtual void SetUp() {
    shared_params_d.resize(2);
    shared_params_d << 2, 0;

    for (std::size_t n = 0; n != N; ++n) {
      x_i[n][0] = n;
      Eigen::VectorXd job_d(2);
      job_d << n + 1.0, n * n;
      job_params_d.push_back(job_d);
    }
  }
};

struct MpiJobSmallWorld : public ::testing::Test {
  Eigen::VectorXd shared_params_d;
  std::vector<Eigen::VectorXd> job_params_d;
  boost::mpi::communicator world_;
  const std::size_t world_size_ = world_.size();
  const std::size_t N = world_size_ > 1 ? world_size_ - 1 : 1;
  std::vector<std::vector<double>> x_r
      = std::vector<std::vector<double>>(N, std::vector<double>(1, 1.0));
  std::vector<std::vector<int>> x_i
      = std::vector<std::vector<int>>(N, std::vector<int>(1, 0));

  virtual void SetUp() {
    shared_params_d.resize(2);
    shared_params_d << 2, 0;

    for (std::size_t n = 0; n != N; ++n) {
      x_i[n][0] = n;
      Eigen::VectorXd job_d(2);
      job_d << n + 1.0, n * n;
      job_params_d.push_back(job_d);
    }
  }
};

TEST_F(MpiJob, hard_work_dd) {
  Eigen::VectorXd result_mpi = stan::math::map_rect<0, hard_work>(
      shared_params_d, job_params_d, x_r, x_i);
  Eigen::VectorXd result_concurrent
      = stan::math::internal::map_rect_concurrent<0, hard_work>(
          shared_params_d, job_params_d, x_r, x_i);

  EXPECT_EQ(result_mpi.rows(), result_concurrent.rows());

  for (int i = 0; i < result_mpi.rows(); ++i) {
    EXPECT_DOUBLE_EQ(result_mpi(i), result_concurrent(i));
  }
}

// execute fewer jobs than the cluster size
TEST_F(MpiJobSmallWorld, hard_work_dd) {
  Eigen::VectorXd result_mpi = stan::math::map_rect<1, hard_work>(
      shared_params_d, job_params_d, x_r, x_i);
  Eigen::VectorXd result_concurrent
      = stan::math::internal::map_rect_concurrent<1, hard_work>(
          shared_params_d, job_params_d, x_r, x_i);

  EXPECT_EQ(result_mpi.rows(), result_concurrent.rows());

  for (int i = 0; i < result_mpi.rows(); ++i) {
    EXPECT_DOUBLE_EQ(result_mpi(i), result_concurrent(i));
  }
}

TEST_F(MpiJob, always_faulty_functor) {
  Eigen::VectorXd result;

  EXPECT_NO_THROW((result = stan::math::map_rect<2, faulty_functor>(
                       shared_params_d, job_params_d, x_r, x_i)));

  // faulty functor throws on theta(0) being -1.0
  // throwing during the first evaluation is quite severe and will
  // lead to a respective runtime error
  job_params_d[0](0) = -1;

  // upon the second evaluation throwing is handled internally
  // different
  EXPECT_THROW_MSG((result = stan::math::map_rect<2, faulty_functor>(
                        shared_params_d, job_params_d, x_r, x_i)),
                   std::domain_error, "Error during MPI evaluation.");

  // throwing on the very first evaluation
  EXPECT_THROW_MSG((result = stan::math::map_rect<3, faulty_functor>(
                        shared_params_d, job_params_d, x_r, x_i)),
                   std::domain_error, "MPI error on first evaluation.");
  // things don't get better if we repeat
  EXPECT_THROW_MSG((result = stan::math::map_rect<3, faulty_functor>(
                        shared_params_d, job_params_d, x_r, x_i)),
                   std::domain_error, "MPI error on first evaluation.");
}

#endif
