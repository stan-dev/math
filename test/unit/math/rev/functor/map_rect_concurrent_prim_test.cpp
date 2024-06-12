// the tests here check that map_rect_concurrent works correct as such we
// enforce that STAN_MPI is NOT defined

#ifndef STAN_MPI

#include <gtest/gtest.h>
#include <stan/math/rev.hpp>

#include <test/unit/math/prim/functor/hard_work.hpp>
#include <test/unit/math/prim/functor/utils_threads.hpp>

#include <iostream>
#include <vector>
#include <string>

STAN_REGISTER_MAP_RECT(0, hard_work)
STAN_REGISTER_MAP_RECT(1, hard_work)
STAN_REGISTER_MAP_RECT(2, hard_work)
STAN_REGISTER_MAP_RECT(3, hard_work)
STAN_REGISTER_MAP_RECT(4, hard_work)
STAN_REGISTER_MAP_RECT(5, hard_work)

struct map_rect_con_prim : public ::testing::Test {
  Eigen::VectorXd shared_params_d;
  std::vector<Eigen::VectorXd> job_params_d;
  std::vector<std::vector<double> > x_r;
  std::vector<std::vector<int> > x_i;
  const int N = 100;

  virtual void SetUp() {
    set_n_threads(4);
    shared_params_d.resize(2);
    shared_params_d << 2, 0;

    for (int n = 0; n != N; ++n) {
      Eigen::VectorXd job_d(2);
      job_d << 0, n * n;
      job_params_d.push_back(job_d);
    }

    x_r = std::vector<std::vector<double> >(N, std::vector<double>(1, 1.0));
    x_i = std::vector<std::vector<int> >(N, std::vector<int>(1, 0));
  }
};

TEST_F(map_rect_con_prim, concurrent_ragged_return_size_dd) {
  Eigen::VectorXd res1 = stan::math::map_rect<0, hard_work>(
      shared_params_d, job_params_d, x_r, x_i);

  EXPECT_EQ(res1.size(), 2 * N);

  x_i[1] = std::vector<int>(1, 1);

  Eigen::VectorXd res2 = stan::math::map_rect<1, hard_work>(
      shared_params_d, job_params_d, x_r, x_i);

  EXPECT_EQ(res2.size(), 2 * N + 1);
}

TEST_F(map_rect_con_prim, concurrent_eval_ok_dd) {
  Eigen::VectorXd res1 = stan::math::map_rect<0, hard_work>(
      shared_params_d, job_params_d, x_r, x_i);
  for (int i = 0, j = 0; i < N; i++) {
    j = 2 * i;
    EXPECT_FLOAT_EQ(
        res1(j), job_params_d[i](0) * job_params_d[i](0) + shared_params_d(0));
    EXPECT_FLOAT_EQ(res1(j + 1),
                    x_r[i][0] * job_params_d[i](1) * job_params_d[i](0)
                        + 2 * shared_params_d(0) + shared_params_d(1));
  }
}

// the tests here check that map_rect_concurrent works correct as such
// we enforce that STAN_MPI is NOT defined

struct map_rect_con : public ::testing::Test {
  Eigen::VectorXd shared_params_d;
  std::vector<Eigen::VectorXd> job_params_d;
  std::vector<std::vector<double> > x_r;
  std::vector<std::vector<int> > x_i;
  const std::size_t N = 100;

  virtual void SetUp() {
    set_n_threads(4);
    shared_params_d.resize(2);
    shared_params_d << 2, 0;

    for (std::size_t n = 0; n != N; ++n) {
      Eigen::VectorXd job_d(2);
      job_d << 0, n * n;
      job_params_d.push_back(job_d);
    }

    x_r = std::vector<std::vector<double> >(N, std::vector<double>(1, 1.0));
    x_i = std::vector<std::vector<int> >(N, std::vector<int>(1, 0));
  }
};

TEST_F(map_rect_con, concurrent_eval_ok_vd) {
  stan::math::vector_v shared_params_v = stan::math::to_var(shared_params_d);
  stan::math::vector_v res1 = stan::math::map_rect<0, hard_work>(
      shared_params_v, job_params_d, x_r, x_i);

  for (std::size_t i = 0, j = 0; i < N; i++) {
    j = 2 * i;
    EXPECT_FLOAT_EQ(
        stan::math::value_of(res1(j)),
        job_params_d[i](0) * job_params_d[i](0) + shared_params_d(0));
    EXPECT_FLOAT_EQ(stan::math::value_of(res1(j + 1)),
                    x_r[i][0] * job_params_d[i](1) * job_params_d[i](0)
                        + 2 * shared_params_d(0) + shared_params_d(1));

    stan::math::set_zero_all_adjoints();
    res1(j).grad();
    EXPECT_FLOAT_EQ(shared_params_v(0).vi_->adj_, 1.0);
    EXPECT_FLOAT_EQ(shared_params_v(1).vi_->adj_, 0.0);

    stan::math::set_zero_all_adjoints();
    res1(j + 1).grad();
    EXPECT_FLOAT_EQ(shared_params_v(0).vi_->adj_, 2.0);
    EXPECT_FLOAT_EQ(shared_params_v(1).vi_->adj_, 1.0);
  }
}

TEST_F(map_rect_con, concurrent_eval_ok_dv) {
  std::vector<stan::math::vector_v> job_params_v;

  for (std::size_t i = 0; i < N; i++)
    job_params_v.push_back(stan::math::to_var(job_params_d[i]));

  stan::math::vector_v res1 = stan::math::map_rect<0, hard_work>(
      shared_params_d, job_params_v, x_r, x_i);

  for (std::size_t i = 0, j = 0; i < N; i++) {
    j = 2 * i;
    EXPECT_FLOAT_EQ(
        stan::math::value_of(res1(j)),
        job_params_d[i](0) * job_params_d[i](0) + shared_params_d(0));
    EXPECT_FLOAT_EQ(stan::math::value_of(res1(j + 1)),
                    x_r[i][0] * job_params_d[i](1) * job_params_d[i](0)
                        + 2 * shared_params_d(0) + shared_params_d(1));

    stan::math::set_zero_all_adjoints();
    res1(j).grad();
    for (std::size_t k = 0; k < N; k++) {
      if (k == i) {
        EXPECT_FLOAT_EQ(job_params_v[i](0).vi_->adj_, 2.0 * job_params_d[i](0));
        EXPECT_FLOAT_EQ(job_params_v[i](1).vi_->adj_, 0.0);
      } else {
        EXPECT_FLOAT_EQ(job_params_v[k](0).vi_->adj_, 0.0);
        EXPECT_FLOAT_EQ(job_params_v[k](1).vi_->adj_, 0.0);
      }
    }

    stan::math::set_zero_all_adjoints();
    res1(j + 1).grad();
    for (std::size_t k = 0; k < N; k++) {
      if (k == i) {
        EXPECT_FLOAT_EQ(job_params_v[i](0).vi_->adj_,
                        x_r[i][0] * job_params_d[i](1));
        EXPECT_FLOAT_EQ(job_params_v[i](1).vi_->adj_,
                        x_r[i][0] * job_params_d[i](0));
      } else {
        EXPECT_FLOAT_EQ(job_params_v[k](0).vi_->adj_, 0.0);
        EXPECT_FLOAT_EQ(job_params_v[k](1).vi_->adj_, 0.0);
      }
    }
  }
}

TEST_F(map_rect_con, concurrent_eval_ok_vv) {
  stan::math::vector_v shared_params_v = stan::math::to_var(shared_params_d);
  std::vector<stan::math::vector_v> job_params_v;

  for (std::size_t i = 0; i < N; i++)
    job_params_v.push_back(stan::math::to_var(job_params_d[i]));

  stan::math::vector_v res1 = stan::math::map_rect<0, hard_work>(
      shared_params_v, job_params_v, x_r, x_i);

  for (std::size_t i = 0, j = 0; i < N; i++) {
    j = 2 * i;
    EXPECT_FLOAT_EQ(
        stan::math::value_of(res1(j)),
        job_params_d[i](0) * job_params_d[i](0) + shared_params_d(0));
    EXPECT_FLOAT_EQ(stan::math::value_of(res1(j + 1)),
                    x_r[i][0] * job_params_d[i](1) * job_params_d[i](0)
                        + 2 * shared_params_d(0) + shared_params_d(1));

    stan::math::set_zero_all_adjoints();
    res1(j).grad();
    EXPECT_FLOAT_EQ(shared_params_v(0).vi_->adj_, 1.0);
    EXPECT_FLOAT_EQ(shared_params_v(1).vi_->adj_, 0.0);

    for (std::size_t k = 0; k < N; k++) {
      if (k == i) {
        EXPECT_FLOAT_EQ(job_params_v[i](0).vi_->adj_, 2.0 * job_params_d[i](0));
        EXPECT_FLOAT_EQ(job_params_v[i](1).vi_->adj_, 0.0);
      } else {
        EXPECT_FLOAT_EQ(job_params_v[k](0).vi_->adj_, 0.0);
        EXPECT_FLOAT_EQ(job_params_v[k](1).vi_->adj_, 0.0);
      }
    }

    stan::math::set_zero_all_adjoints();
    res1(j + 1).grad();
    EXPECT_FLOAT_EQ(shared_params_v(0).vi_->adj_, 2.0);
    EXPECT_FLOAT_EQ(shared_params_v(1).vi_->adj_, 1.0);

    for (std::size_t k = 0; k < N; k++) {
      if (k == i) {
        EXPECT_FLOAT_EQ(job_params_v[i](0).vi_->adj_,
                        x_r[i][0] * job_params_d[i](1));
        EXPECT_FLOAT_EQ(job_params_v[i](1).vi_->adj_,
                        x_r[i][0] * job_params_d[i](0));
      } else {
        EXPECT_FLOAT_EQ(job_params_v[k](0).vi_->adj_, 0.0);
        EXPECT_FLOAT_EQ(job_params_v[k](1).vi_->adj_, 0.0);
      }
    }
  }
}

// the tests here check that map_rect_concurrent works correct as such we
// enforce that STAN_MPI is NOT defined and these tests only run if
// STAN_THREADS is defined

#ifdef STAN_THREADS

inline void setup_job(int N, Eigen::VectorXd& shared_params_d,
                      std::vector<Eigen::VectorXd>& job_params_d,
                      std::vector<std::vector<double> >& x_r,
                      std::vector<std::vector<int> >& x_i) {
  shared_params_d.resize(2);
  shared_params_d << 2, 0;

  job_params_d.clear();
  for (int n = 0; n != N; ++n) {
    Eigen::VectorXd job_d(2);
    job_d << 1.1, n * n;
    job_params_d.push_back(job_d);
  }

  x_r.clear();
  x_i.clear();

  if (N > 0) {
    x_r = std::vector<std::vector<double> >(N, std::vector<double>(1, 1.0));
    x_i = std::vector<std::vector<int> >(N, std::vector<int>(1, 0));
  }
}

struct map_rect_con_threads : public ::testing::Test {
  Eigen::VectorXd shared_params_d;
  std::vector<Eigen::VectorXd> job_params_d;
  std::vector<std::vector<double> > x_r;
  std::vector<std::vector<int> > x_i;
  std::vector<int> N_test{0, 1, 2, 3, 4, 5, 7, 10, 13, 67, 100};

  virtual void SetUp() {}
};

TEST_F(map_rect_con_threads, concurrent_varying_num_threads_ragged_dd) {
  set_n_threads(-1);

  for (std::size_t i = 1; i < 20; ++i) {
    for (std::size_t n = 0; n < N_test.size(); ++n) {
      const int N = N_test[n];

      setup_job(N, shared_params_d, job_params_d, x_r, x_i);

      Eigen::VectorXd res1 = stan::math::map_rect<0, hard_work>(
          shared_params_d, job_params_d, x_r, x_i);
      EXPECT_EQ(res1.size(), 2 * N);
      for (int i = 0, j = 0; i < N; i++) {
        j = 2 * i;
        EXPECT_FLOAT_EQ(res1(j), job_params_d[i](0) * job_params_d[i](0)
                                     + shared_params_d(0));
        EXPECT_FLOAT_EQ(res1(j + 1),
                        x_r[i][0] * job_params_d[i](1) * job_params_d[i](0)
                            + 2 * shared_params_d(0) + shared_params_d(1));
      }
    }
    set_n_threads(i);
  }

  set_n_threads(-1);
  for (std::size_t i = 1; i < 20; ++i) {
    for (std::size_t n = 1; n < N_test.size(); ++n) {
      const int N = N_test[n];
      setup_job(N, shared_params_d, job_params_d, x_r, x_i);

      x_i[0] = std::vector<int>(1, 1);

      Eigen::VectorXd res2 = stan::math::map_rect<0, hard_work>(
          shared_params_d, job_params_d, x_r, x_i);
      EXPECT_EQ(res2.size(), 2 * N + 1);
    }
    set_n_threads(i);
  }
}

TEST_F(map_rect_con_threads, concurrent_varying_num_threads_eval_ok_dd) {
  set_n_threads(-1);
  for (std::size_t i = 1; i < 20; ++i) {
    for (std::size_t n = 0; n < N_test.size(); ++n) {
      const int N = N_test[n];
      setup_job(N, shared_params_d, job_params_d, x_r, x_i);

      Eigen::VectorXd res1 = stan::math::map_rect<0, hard_work>(
          shared_params_d, job_params_d, x_r, x_i);
      EXPECT_EQ(res1.size(), 2 * N);
      for (int i = 0, j = 0; i < N; i++) {
        j = 2 * i;
        EXPECT_FLOAT_EQ(res1(j), job_params_d[i](0) * job_params_d[i](0)
                                     + shared_params_d(0));
        EXPECT_FLOAT_EQ(res1(j + 1),
                        x_r[i][0] * job_params_d[i](1) * job_params_d[i](0)
                            + 2 * shared_params_d(0) + shared_params_d(1));
      }
    }
    set_n_threads(i);
  }
}
#endif

// the tests here check that map_rect refuses mal-formatted input as
// such it does not matter if STAN_MPI is defined or not

struct map_rect_prim : public ::testing::Test {
  Eigen::VectorXd shared_params_d;
  std::vector<Eigen::VectorXd> job_params_d;
  std::vector<std::vector<double> > x_r;
  std::vector<std::vector<int> > x_i;
  const std::size_t N = 10;

  virtual void SetUp() {
    shared_params_d.resize(2);
    shared_params_d << 2, 0;

    for (std::size_t n = 0; n != N; ++n) {
      Eigen::VectorXd job_d(2);
      job_d << 0, n * n;
      job_params_d.push_back(job_d);
    }

    x_r = std::vector<std::vector<double> >(N, std::vector<double>(1, 1.0));
    x_i = std::vector<std::vector<int> >(N, std::vector<int>(1, 0));
  }
};

TEST_F(map_rect_prim, nompi_no_job_input_ok_dd) {
  job_params_d.resize(0);
  x_i.resize(0);
  x_r.resize(0);

  EXPECT_NO_THROW((stan::math::map_rect<0, hard_work>(shared_params_d,
                                                      job_params_d, x_r, x_i)));
}

TEST_F(map_rect_prim, nompi_size_mismatch_job_params_dd) {
  job_params_d.pop_back();

  EXPECT_THROW((stan::math::map_rect<1, hard_work>(shared_params_d,
                                                   job_params_d, x_r, x_i)),
               std::invalid_argument);
}

TEST_F(map_rect_prim, nompi_size_mismatch_real_data_dd) {
  x_r.pop_back();

  EXPECT_THROW((stan::math::map_rect<2, hard_work>(shared_params_d,
                                                   job_params_d, x_r, x_i)),
               std::invalid_argument);
}

TEST_F(map_rect_prim, nompi_size_mismatch_int_data_dd) {
  x_i.pop_back();

  EXPECT_THROW((stan::math::map_rect<3, hard_work>(shared_params_d,
                                                   job_params_d, x_r, x_i)),
               std::invalid_argument);
}

TEST_F(map_rect_prim, nompi_wrong_size_job_params_dd) {
  job_params_d[1].resize(5);

  EXPECT_THROW((stan::math::map_rect<1, hard_work>(shared_params_d,
                                                   job_params_d, x_r, x_i)),
               std::invalid_argument);
}

TEST_F(map_rect_prim, nompi_wrong_size_real_data_dd) {
  x_r[1] = std::vector<double>(5, 1);

  EXPECT_THROW((stan::math::map_rect<4, hard_work>(shared_params_d,
                                                   job_params_d, x_r, x_i)),
               std::invalid_argument);
}

TEST_F(map_rect_prim, nompi_wrong_size_int_data_dd) {
  x_i[1] = std::vector<int>(5, 1);

  EXPECT_THROW((stan::math::map_rect<5, hard_work>(shared_params_d,
                                                   job_params_d, x_r, x_i)),
               std::invalid_argument);
}
#endif
