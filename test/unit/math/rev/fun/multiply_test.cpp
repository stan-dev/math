#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/rev/mat/util.hpp>

#ifdef STAN_OPENCL
#include <boost/random/mersenne_twister.hpp>
#define EXPECT_MATRIX_NEAR(A, B, DELTA) \
  for (int i = 0; i < A.size(); i++)    \
    EXPECT_NEAR(stan::math::value_of(A(i)), stan::math::value_of(B(i)), DELTA);

boost::random::mt19937 rng;
#define MULTIPLY_OPENCL_OVERRIDE 0
#define MULTIPLY_CPU_OVERRIDE INT_MAX
TEST(AgradRevMatrix, multiply_val_vv_cl) {
  int temp = stan::math::opencl_context.tuning_opts()
                 .multiply_dim_prod_worth_transfer;
  using stan::math::matrix_d;
  using stan::math::matrix_v;
  using stan::math::multiply;
  int mat_size = 234;
  matrix_v Av(mat_size, mat_size);
  matrix_v Bv(mat_size, mat_size);
  matrix_v C, C_cl;
  for (int i = 0; i < mat_size; i++) {
    for (int j = 0; j < mat_size; j++) {
      Av(i, j) = stan::math::uniform_rng(-5, 5, rng);
      Bv(i, j) = stan::math::uniform_rng(-5, 5, rng);
    }
  }
  stan::math::opencl_context.tuning_opts().multiply_dim_prod_worth_transfer
      = MULTIPLY_OPENCL_OVERRIDE;
  C_cl = multiply(Av, Bv);
  C_cl(0, 0).grad();
  stan::math::opencl_context.tuning_opts().multiply_dim_prod_worth_transfer
      = MULTIPLY_CPU_OVERRIDE;
  C = multiply(Av, Bv);
  C(0, 0).grad();
  EXPECT_MATRIX_NEAR(C, C_cl, 1.0E-12);
  EXPECT_MATRIX_NEAR(C.adj(), C_cl.adj(), 1.0E-12);
  stan::math::recover_memory();
  stan::math::opencl_context.tuning_opts().multiply_dim_prod_worth_transfer
      = temp;
}

TEST(AgradRevMatrix, multiply_val_vd_cl) {
  int temp = stan::math::opencl_context.tuning_opts()
                 .multiply_dim_prod_worth_transfer;
  using stan::math::matrix_d;
  using stan::math::matrix_v;
  using stan::math::multiply;
  int mat_size = 256;
  matrix_v Av(mat_size, mat_size);
  matrix_v Bd(mat_size, mat_size);
  matrix_v C, C_cl;
  for (int i = 0; i < mat_size; i++) {
    for (int j = 0; j < mat_size; j++) {
      Av(i, j) = stan::math::uniform_rng(-5, 5, rng);
      Bd(i, j) = stan::math::uniform_rng(-5, 5, rng);
    }
  }

  stan::math::opencl_context.tuning_opts().multiply_dim_prod_worth_transfer
      = MULTIPLY_OPENCL_OVERRIDE;
  C_cl = multiply(Av, Bd);
  C_cl(0, 0).grad();
  stan::math::opencl_context.tuning_opts().multiply_dim_prod_worth_transfer
      = MULTIPLY_CPU_OVERRIDE;
  C = multiply(Av, Bd);
  C(0, 0).grad();
  EXPECT_MATRIX_NEAR(C, C_cl, 1.0E-12);
  EXPECT_MATRIX_NEAR(C.adj(), C_cl.adj(), 1.0E-12);
  stan::math::recover_memory();
  stan::math::opencl_context.tuning_opts().multiply_dim_prod_worth_transfer
      = temp;
}

TEST(AgradRevMatrix, multiply_val_dv_cl) {
  int temp = stan::math::opencl_context.tuning_opts()
                 .multiply_dim_prod_worth_transfer;
  using stan::math::matrix_d;
  using stan::math::matrix_v;
  using stan::math::multiply;
  int mat_size = 321;
  matrix_v Ad(mat_size, mat_size);
  matrix_v Bv(mat_size, mat_size);
  matrix_v C, C_cl;
  for (int i = 0; i < mat_size; i++) {
    for (int j = 0; j < mat_size; j++) {
      Ad(i, j) = stan::math::uniform_rng(-5, 5, rng);
      Bv(i, j) = stan::math::uniform_rng(-5, 5, rng);
    }
  }
  stan::math::opencl_context.tuning_opts().multiply_dim_prod_worth_transfer
      = MULTIPLY_OPENCL_OVERRIDE;
  C_cl = multiply(Ad, Bv);
  C_cl(0, 0).grad();
  stan::math::opencl_context.tuning_opts().multiply_dim_prod_worth_transfer
      = MULTIPLY_CPU_OVERRIDE;
  C = multiply(Ad, Bv);
  C(0, 0).grad();
  EXPECT_MATRIX_NEAR(C, C_cl, 1.0E-12);
  EXPECT_MATRIX_NEAR(C.adj(), C_cl.adj(), 1.0E-12);
  stan::math::recover_memory();
  stan::math::opencl_context.tuning_opts().multiply_dim_prod_worth_transfer
      = temp;
}

#endif
