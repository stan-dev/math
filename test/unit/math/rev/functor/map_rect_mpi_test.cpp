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
STAN_REGISTER_MAP_RECT(1, faulty_functor)
STAN_REGISTER_MAP_RECT(2, faulty_functor)

struct MpiJob : public ::testing::Test {
  stan::math::vector_d shared_params_d;
  std::vector<stan::math::vector_d> job_params_d;
  stan::math::vector_v shared_params_v;
  std::vector<stan::math::vector_v> job_params_v;
  stan::math::vector_v shared_params_v2;
  std::vector<stan::math::vector_v> job_params_v2;
  const std::size_t N = 10;
  std::vector<std::vector<double>> x_r
      = std::vector<std::vector<double>>(N, std::vector<double>(1, 1.0));
  std::vector<std::vector<int>> x_i
      = std::vector<std::vector<int>>(N, std::vector<int>(1, 0));

  virtual void SetUp() {
    shared_params_v.resize(2);
    shared_params_v << 2, 0;
    shared_params_d = stan::math::value_of(shared_params_v);

    shared_params_v2.resize(2);
    shared_params_v2 << 2, 0;

    for (std::size_t n = 0; n != N; ++n) {
      x_i[n][0] = n;
      stan::math::vector_v job_v(2);
      job_v << n + 1.0, n * n;
      job_params_v.push_back(job_v);

      job_params_d.push_back(stan::math::value_of(job_v));

      stan::math::vector_v job_v2(2);
      job_v2 << n + 1.0, n * n;
      job_params_v2.push_back(job_v2);
    }
  }
};

TEST_F(MpiJob, hard_work_vv) {
  std::vector<stan::math::var> shared_params_v_vec
      = stan::math::to_array_1d(shared_params_v);
  std::vector<stan::math::var> shared_params_v2_vec
      = stan::math::to_array_1d(shared_params_v2);

  std::vector<std::vector<stan::math::var>> job_params_v_vec;
  std::vector<std::vector<stan::math::var>> job_params_v2_vec;

  for (std::size_t i = 0; i < N; ++i) {
    job_params_v_vec.push_back(stan::math::to_array_1d(job_params_v[i]));
    job_params_v2_vec.push_back(stan::math::to_array_1d(job_params_v2[i]));
  }

  stan::math::vector_v result_mpi = stan::math::map_rect<0, hard_work>(
      shared_params_v, job_params_v, x_r, x_i, 0);

  stan::math::vector_v result_concurrent
      = stan::math::internal::map_rect_concurrent<0, hard_work>(
          shared_params_v2, job_params_v2, x_r, x_i, 0);

  std::vector<double> z_grad1;
  std::vector<double> z_grad2;

  EXPECT_EQ(result_mpi.rows(), result_concurrent.rows());

  for (std::size_t i = 0, ij = 0; i < job_params_v_vec.size(); ++i) {
    for (std::size_t j = 0; j < 2; ++j, ++ij) {
      EXPECT_DOUBLE_EQ(stan::math::value_of(result_mpi(ij)),
                       stan::math::value_of(result_concurrent(ij)));

      std::vector<stan::math::var> z_var1, z_var2;

      z_var1.insert(z_var1.end(), shared_params_v_vec.begin(),
                    shared_params_v_vec.end());
      z_var1.insert(z_var1.end(), job_params_v_vec[i].begin(),
                    job_params_v_vec[i].end());

      z_var2.insert(z_var2.end(), shared_params_v2_vec.begin(),
                    shared_params_v2_vec.end());
      z_var2.insert(z_var2.end(), job_params_v2_vec[i].begin(),
                    job_params_v2_vec[i].end());

      stan::math::set_zero_all_adjoints();

      result_mpi(ij).grad(z_var1, z_grad1);
      result_concurrent(ij).grad(z_var2, z_grad2);

      for (std::size_t k = 0; k < z_grad1.size(); ++k) {
        EXPECT_DOUBLE_EQ(z_grad1[k], z_grad2[k]);
      }
    }
  }
}

TEST_F(MpiJob, always_faulty_functor_vv) {
  stan::math::vector_v result;

  EXPECT_NO_THROW((result = stan::math::map_rect<1, faulty_functor>(
                       shared_params_v, job_params_v, x_r, x_i)));

  // faulty functor throws on theta(0) being -1.0
  // throwing during the first evaluation is quite severe and will
  // lead to a respective runtime error
  job_params_v[0](0) = -1;

  // upon the second evaluation throwing is handled internally different
  EXPECT_THROW_MSG((result = stan::math::map_rect<1, faulty_functor>(
                        shared_params_v, job_params_v, x_r, x_i)),
                   std::domain_error, "Error during MPI evaluation.");

  // throwing on the very first evaluation
  EXPECT_THROW_MSG((result = stan::math::map_rect<2, faulty_functor>(
                        shared_params_v, job_params_v, x_r, x_i)),
                   std::domain_error, "MPI error on first evaluation.");
}

TEST_F(MpiJob, always_faulty_functor_vd) {
  stan::math::vector_v result;

  EXPECT_NO_THROW((result = stan::math::map_rect<1, faulty_functor>(
                       shared_params_v, job_params_d, x_r, x_i)));

  // faulty functor throws on theta(0) being -1.0
  // throwing during the first evaluation is quite severe and will
  // lead to a respective runtime error
  job_params_d[0](0) = -1;

  // upon the second evaluation throwing is handled internally different
  EXPECT_THROW_MSG((result = stan::math::map_rect<1, faulty_functor>(
                        shared_params_v, job_params_d, x_r, x_i)),
                   std::domain_error, "Error during MPI evaluation.");

  // throwing on the very first evaluation
  EXPECT_THROW_MSG((result = stan::math::map_rect<2, faulty_functor>(
                        shared_params_v, job_params_d, x_r, x_i)),
                   std::domain_error, "MPI error on first evaluation.");
}

TEST_F(MpiJob, always_faulty_functor_dv) {
  stan::math::vector_v result;

  EXPECT_NO_THROW((result = stan::math::map_rect<1, faulty_functor>(
                       shared_params_d, job_params_v, x_r, x_i)));

  // faulty functor throws on theta(0) being -1.0
  // throwing during the first evaluation is quite severe and will
  // lead to a respective runtime error
  job_params_v[0](0) = -1;

  // upon the second evaluation throwing is handled internally different
  EXPECT_THROW_MSG((result = stan::math::map_rect<1, faulty_functor>(
                        shared_params_d, job_params_v, x_r, x_i)),
                   std::domain_error, "Error during MPI evaluation.");

  // throwing on the very first evaluation
  EXPECT_THROW_MSG((result = stan::math::map_rect<2, faulty_functor>(
                        shared_params_d, job_params_v, x_r, x_i)),
                   std::domain_error, "MPI error on first evaluation.");
}

#endif
