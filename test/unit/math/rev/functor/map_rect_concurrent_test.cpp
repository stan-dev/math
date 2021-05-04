// the tests here check that map_rect_concurrent works correct as such
// we enforce that STAN_MPI is NOT defined

#ifdef STAN_MPI
#undef STAN_MPI
#endif

#include <gtest/gtest.h>
#include <stan/math/rev.hpp>

#include <test/unit/math/prim/functor/hard_work.hpp>
#include <test/unit/math/prim/functor/utils_threads.hpp>

#include <iostream>
#include <vector>

STAN_REGISTER_MAP_RECT(0, hard_work)

struct map_rect : public ::testing::Test {
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

TEST_F(map_rect, concurrent_eval_ok_vd) {
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

TEST_F(map_rect, concurrent_eval_ok_dv) {
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

TEST_F(map_rect, concurrent_eval_ok_vv) {
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
