// these tests can only be compiled and executed with availability of
// MPI
#ifdef STAN_MPI

#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>

#include <test/unit/math/prim/mat/functor/faulty_functor.hpp>

#include <iostream>
#include <vector>

using stan::math::matrix_d;
using stan::math::vector_d;

struct mock_reduce {
  matrix_d operator()(const vector_d& shared_params,
                      const vector_d& job_specific_params,
                      const std::vector<double>& x_r,
                      const std::vector<int>& x_i,
                      std::ostream* msgs = nullptr) const {
    boost::mpi::communicator world;
    if (stan::math::abs(job_specific_params(0) + 1.0)
        < 1E-7) {  // check for param being 1.0
      throw std::domain_error("Illegal parameter!");
    }
    return stan::math::rep_matrix(world.rank(), 3, world.rank());
  }
};

template <typename F, typename T_shared_param, typename T_job_param>
struct mock_combine {
 public:
  typedef matrix_d result_t;

  mock_combine() {}
  mock_combine(
      const Eigen::Matrix<T_shared_param, Eigen::Dynamic, 1>& shared_params,
      const std::vector<Eigen::Matrix<T_job_param, Eigen::Dynamic, 1>>&
          job_params) {}

  result_t operator()(const matrix_d& local_result,
                      const std::vector<int>& world_f_out) {
    return local_result;
  }
};

typedef mock_combine<faulty_functor, double, double> mock_combine_dd;

typedef stan::math::mpi_parallel_call<0, mock_reduce, mock_combine_dd>
    mock0_call_t;
STAN_REGISTER_MPI_DISTRIBUTED_APPLY(mock0_call_t)

typedef stan::math::mpi_parallel_call<1, mock_reduce, mock_combine_dd>
    mock1_call_t;
STAN_REGISTER_MPI_DISTRIBUTED_APPLY(mock1_call_t)

typedef stan::math::mpi_parallel_call<2, mock_reduce, mock_combine_dd>
    mock_call_t;
STAN_REGISTER_MPI_DISTRIBUTED_APPLY(mock_call_t)

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

TEST_F(MpiJob, no_job_input_ok_dd) {
  job_params_d.resize(0);
  x_i.resize(0);
  x_r.resize(0);

  std::shared_ptr<mock0_call_t> call;
  EXPECT_NO_THROW((call = std::shared_ptr<mock0_call_t>(new mock0_call_t(
                       shared_params_d, job_params_d, x_r, x_i))));

  matrix_d res = call->reduce_combine();

  EXPECT_EQ(res.size(), 0);
}

TEST_F(MpiJob, one_job_input_ok_dd) {
  job_params_d.resize(1);
  x_i.resize(1);
  x_r.resize(1);

  std::shared_ptr<mock1_call_t> call;
  EXPECT_NO_THROW((call = std::shared_ptr<mock1_call_t>(new mock1_call_t(
                       shared_params_d, job_params_d, x_r, x_i))));

  matrix_d res = call->reduce_combine();

  EXPECT_EQ(res.size(), 0);
  EXPECT_EQ(res.cols(), 0);
  EXPECT_EQ(res.rows(), 3);
}

TEST_F(MpiJob, size_mismatch_job_params_dd) {
  job_params_d.pop_back();

  std::shared_ptr<mock_call_t> call;
  EXPECT_THROW((call = std::shared_ptr<mock_call_t>(
                    new mock_call_t(shared_params_d, job_params_d, x_r, x_i))),
               std::invalid_argument);
}

TEST_F(MpiJob, size_mismatch_real_data_dd) {
  x_r.pop_back();

  std::shared_ptr<mock_call_t> call;
  EXPECT_THROW((call = std::shared_ptr<mock_call_t>(
                    new mock_call_t(shared_params_d, job_params_d, x_r, x_i))),
               std::invalid_argument);
}

TEST_F(MpiJob, size_mismatch_int_data_dd) {
  x_i.pop_back();

  std::shared_ptr<mock_call_t> call;
  EXPECT_THROW((call = std::shared_ptr<mock_call_t>(
                    new mock_call_t(shared_params_d, job_params_d, x_r, x_i))),
               std::invalid_argument);
}

TEST_F(MpiJob, recover_on_first_evaluation_dd) {
  // trigger exception on first evaluation
  job_params_d[0](0) = -1.0;

  std::shared_ptr<mock_call_t> call;
  EXPECT_NO_THROW((call = std::shared_ptr<mock_call_t>(new mock_call_t(
                       shared_params_d, job_params_d, x_r, x_i))));

  EXPECT_THROW(call->reduce_combine(), std::domain_error);

  // get rid of call object to free MPI resource
  call.reset();

  // if things are again ok, then this works
  job_params_d[0](0) = 1.0;
  EXPECT_NO_THROW((call = std::shared_ptr<mock_call_t>(new mock_call_t(
                       shared_params_d, job_params_d, x_r, x_i))));

  EXPECT_NO_THROW(call->reduce_combine());

  // note: tests below have data already cached for mock_call_t !
}

TEST_F(MpiJob, MPI_busy) {
  std::shared_ptr<mock_call_t> call1;
  EXPECT_NO_THROW((call1 = std::shared_ptr<mock_call_t>(new mock_call_t(
                       shared_params_d, job_params_d, x_r, x_i))));

  EXPECT_NO_THROW(call1->reduce_combine());

  // we still hold a reference to mpi_parallel_call which blocks the
  // MPI resource until its deallocated
  std::shared_ptr<mock_call_t> call2;
  EXPECT_THROW((call2 = std::shared_ptr<mock_call_t>(
                    new mock_call_t(shared_params_d, job_params_d, x_r, x_i))),
               stan::math::mpi_is_in_use);

  // free it up
  call1.reset();

  // now it works
  EXPECT_NO_THROW((call2 = std::shared_ptr<mock_call_t>(new mock_call_t(
                       shared_params_d, job_params_d, x_r, x_i))));
  EXPECT_NO_THROW(call2->reduce_combine());
}

TEST_F(MpiJob, size_mismatch_cached_jobs_dd) {
  // run call once to cache data
  std::shared_ptr<mock_call_t> call;
  EXPECT_NO_THROW((call = std::shared_ptr<mock_call_t>(new mock_call_t(
                       shared_params_d, job_params_d, x_r, x_i))));

  call->reduce_combine();

  job_params_d.pop_back();

  // the # of jobs must not change when calling again
  EXPECT_THROW((call = std::shared_ptr<mock_call_t>(
                    new mock_call_t(shared_params_d, job_params_d, x_r, x_i))),
               std::invalid_argument);
}

TEST_F(MpiJob, root_not_confused_dd) {
  // the root must not call the distributed_apply ever
  EXPECT_THROW_MSG(mock_call_t::distributed_apply(), std::runtime_error,
                   "problem sizes must be defined on the root.");
}

#endif
