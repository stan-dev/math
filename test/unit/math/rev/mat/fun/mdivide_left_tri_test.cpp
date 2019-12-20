#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/rev/mat/util.hpp>

// TODO(carpenter): move this to test framework;  should be able to put
// all this GPU config into the functor
// TODO(carpenter): safer cleanup on resetting tuning_opts

#ifdef STAN_OPENCL
#include <boost/random/mersenne_twister.hpp>

#define EXPECT_MATRIX_NEAR(A, B, DELTA) \
  for (int i = 0; i < A.size(); i++)    \
    EXPECT_NEAR(stan::math::value_of(A(i)), stan::math::value_of(B(i)), DELTA);

boost::random::mt19937 rng;
#define MDIVIDE_OPENCL_OVERRIDE 0
#define MDIVIDE_CPU_OVERRIDE INT_MAX
TEST(AgradRevMatrix, mdivide_left_tri_val_cl) {
  int temp = stan::math::opencl_context.tuning_opts()
                 .tri_inverse_size_worth_transfer;
  using stan::math::matrix_d;
  using stan::math::matrix_v;
  using stan::math::mdivide_left_tri;
  int mat_size = 111;
  matrix_v Av(mat_size, mat_size);
  matrix_d Ad(mat_size, mat_size);
  matrix_v I, I_cl;
  for (int i = 0; i < mat_size; i++) {
    for (int j = 0; j < i; j++) {
      Ad(i, j) = stan::math::uniform_rng(-5, 5, rng);
      Av(i, j) = Ad(i, j);
    }
    Ad(i, i) = 20.0;
    Av(i, i) = Ad(i, i);
    for (int j = i + 1; j < mat_size; j++) {
      Ad(i, j) = 0.0;
      Av(i, j) = Ad(i, j);
    }
  }
  stan::math::opencl_context.tuning_opts().tri_inverse_size_worth_transfer
      = MDIVIDE_OPENCL_OVERRIDE;
  I_cl = mdivide_left_tri<Eigen::Lower>(Av, Av);
  stan::math::opencl_context.tuning_opts().tri_inverse_size_worth_transfer
      = MDIVIDE_CPU_OVERRIDE;
  I = mdivide_left_tri<Eigen::Lower>(Av, Av);
  EXPECT_MATRIX_NEAR(I, I_cl, 1.0E-12);

  stan::math::opencl_context.tuning_opts().tri_inverse_size_worth_transfer
      = MDIVIDE_OPENCL_OVERRIDE;
  I_cl = mdivide_left_tri<Eigen::Lower>(Av, Ad);
  stan::math::opencl_context.tuning_opts().tri_inverse_size_worth_transfer
      = MDIVIDE_CPU_OVERRIDE;
  I = mdivide_left_tri<Eigen::Lower>(Av, Ad);
  EXPECT_MATRIX_NEAR(I, I_cl, 1.0E-12);

  stan::math::opencl_context.tuning_opts().tri_inverse_size_worth_transfer
      = MDIVIDE_OPENCL_OVERRIDE;
  I_cl = mdivide_left_tri<Eigen::Lower>(Ad, Av);
  stan::math::opencl_context.tuning_opts().tri_inverse_size_worth_transfer
      = MDIVIDE_CPU_OVERRIDE;
  I = mdivide_left_tri<Eigen::Lower>(Ad, Av);
  EXPECT_MATRIX_NEAR(I, I_cl, 1.0E-12);

  stan::math::opencl_context.tuning_opts().tri_inverse_size_worth_transfer
      = MDIVIDE_OPENCL_OVERRIDE;
  matrix_v Av_inv_cl = mdivide_left_tri<Eigen::Lower>(Av);
  I_cl = Av * Av_inv_cl;
  stan::math::opencl_context.tuning_opts().tri_inverse_size_worth_transfer
      = MDIVIDE_CPU_OVERRIDE;
  matrix_v Av_inv = mdivide_left_tri<Eigen::Lower>(Av);
  I = Av * Av_inv;
  EXPECT_MATRIX_NEAR(I, I_cl, 1.0E-12);

  for (int i = 0; i < mat_size; i++) {
    for (int j = 0; j < i; j++) {
      Ad(i, j) = 0.0;
      Av(i, j) = Ad(i, j);
    }
    Ad(i, i) = 20.0;
    Av(i, i) = Ad(i, i);
    for (int j = i + 1; j < mat_size; j++) {
      Ad(i, j) = stan::math::uniform_rng(-5, 5, rng);
      Av(i, j) = Ad(i, j);
    }
  }

  stan::math::opencl_context.tuning_opts().tri_inverse_size_worth_transfer
      = MDIVIDE_OPENCL_OVERRIDE;
  I_cl = mdivide_left_tri<Eigen::Upper>(Av, Av);
  stan::math::opencl_context.tuning_opts().tri_inverse_size_worth_transfer
      = MDIVIDE_CPU_OVERRIDE;
  I = mdivide_left_tri<Eigen::Upper>(Av, Av);
  EXPECT_MATRIX_NEAR(I, I_cl, 1.0E-12);

  stan::math::opencl_context.tuning_opts().tri_inverse_size_worth_transfer
      = MDIVIDE_OPENCL_OVERRIDE;
  I_cl = mdivide_left_tri<Eigen::Upper>(Av, Ad);
  stan::math::opencl_context.tuning_opts().tri_inverse_size_worth_transfer
      = MDIVIDE_CPU_OVERRIDE;
  I = mdivide_left_tri<Eigen::Upper>(Av, Ad);
  EXPECT_MATRIX_NEAR(I, I_cl, 1.0E-12);

  stan::math::opencl_context.tuning_opts().tri_inverse_size_worth_transfer
      = MDIVIDE_OPENCL_OVERRIDE;
  I_cl = mdivide_left_tri<Eigen::Upper>(Ad, Av);
  stan::math::opencl_context.tuning_opts().tri_inverse_size_worth_transfer
      = MDIVIDE_CPU_OVERRIDE;
  I = mdivide_left_tri<Eigen::Upper>(Ad, Av);
  EXPECT_MATRIX_NEAR(I, I_cl, 1.0E-12);

  // restore default value
  stan::math::opencl_context.tuning_opts().tri_inverse_size_worth_transfer
      = temp;
}

TEST(AgradRevMatrix, mdivide_left_tri_lower_grad_vv_cl) {
  int temp = stan::math::opencl_context.tuning_opts()
                 .tri_inverse_size_worth_transfer;

  using stan::math::matrix_d;
  using stan::math::matrix_v;
  using stan::math::mdivide_left_tri;
  using stan::math::multiply;
  int mat_size = 256;

  matrix_v Av(mat_size, mat_size);
  matrix_v Bv(mat_size, mat_size);
  matrix_v Cv(mat_size, mat_size);
  matrix_v Cv_cl(mat_size, mat_size);
  for (int i = 0; i < mat_size; i++) {
    for (int j = 0; j < i; j++) {
      Av(i, j) = stan::math::uniform_rng(-5, 5, rng);
    }
    Av(i, i) = 20.0;
    for (int j = i + 1; j < mat_size; j++) {
      Av(i, j) = 0.0;
    }
    for (int j = 0; j < mat_size; j++) {
      Bv(i, j) = stan::math::uniform_rng(-5, 5, rng);
    }
  }
  stan::math::opencl_context.tuning_opts().tri_inverse_size_worth_transfer
      = MDIVIDE_OPENCL_OVERRIDE;
  Cv_cl = mdivide_left_tri<Eigen::Lower>(Av, Bv);
  Cv_cl(0, 0).grad();
  stan::math::opencl_context.tuning_opts().tri_inverse_size_worth_transfer
      = MDIVIDE_CPU_OVERRIDE;
  Cv = mdivide_left_tri<Eigen::Lower>(Av, Bv);
  Cv(0, 0).grad();
  EXPECT_MATRIX_NEAR(Cv, Cv_cl, 1.0E-12);
  EXPECT_MATRIX_NEAR(Cv.adj(), Cv_cl.adj(), 1.0E-12);
  stan::math::recover_memory();
  stan::math::opencl_context.tuning_opts().tri_inverse_size_worth_transfer
      = temp;
}

TEST(AgradRevMatrix, mdivide_left_tri_lower_grad_dv_cl) {
  int temp = stan::math::opencl_context.tuning_opts()
                 .tri_inverse_size_worth_transfer;

  using stan::math::matrix_d;
  using stan::math::matrix_v;
  using stan::math::mdivide_left_tri;
  using stan::math::multiply;
  int mat_size = 213;

  matrix_d Ad(mat_size, mat_size);
  matrix_v Bv(mat_size, mat_size);
  matrix_v Cv(mat_size, mat_size);
  matrix_v Cv_cl(mat_size, mat_size);
  for (int i = 0; i < mat_size; i++) {
    for (int j = 0; j < i; j++) {
      Ad(i, j) = stan::math::uniform_rng(-5, 5, rng);
    }
    Ad(i, i) = 20.0;
    for (int j = i + 1; j < mat_size; j++) {
      Ad(i, j) = 0.0;
    }
    for (int j = 0; j < mat_size; j++) {
      Bv(i, j) = stan::math::uniform_rng(-5, 5, rng);
    }
  }
  stan::math::opencl_context.tuning_opts().tri_inverse_size_worth_transfer
      = MDIVIDE_OPENCL_OVERRIDE;
  Cv_cl = mdivide_left_tri<Eigen::Lower>(Ad, Bv);
  Cv_cl(0, 0).grad();
  stan::math::opencl_context.tuning_opts().tri_inverse_size_worth_transfer
      = MDIVIDE_CPU_OVERRIDE;
  Cv = mdivide_left_tri<Eigen::Lower>(Ad, Bv);
  Cv(0, 0).grad();
  EXPECT_MATRIX_NEAR(Cv, Cv_cl, 1.0E-12);
  EXPECT_MATRIX_NEAR(Cv.adj(), Cv_cl.adj(), 1.0E-12);
  stan::math::recover_memory();
  stan::math::opencl_context.tuning_opts().tri_inverse_size_worth_transfer
      = temp;
}

TEST(AgradRevMatrix, mdivide_left_tri_lower_grad_vd_cl) {
  int temp = stan::math::opencl_context.tuning_opts()
                 .tri_inverse_size_worth_transfer;

  using stan::math::matrix_d;
  using stan::math::matrix_v;
  using stan::math::mdivide_left_tri;
  using stan::math::multiply;
  int mat_size = 250;

  matrix_d Av(mat_size, mat_size);
  matrix_v Bd(mat_size, mat_size);
  matrix_v Cv(mat_size, mat_size);
  matrix_v Cv_cl(mat_size, mat_size);
  for (int i = 0; i < mat_size; i++) {
    for (int j = 0; j < i; j++) {
      Av(i, j) = stan::math::uniform_rng(-5, 5, rng);
    }
    Av(i, i) = 20.0;
    for (int j = i + 1; j < mat_size; j++) {
      Av(i, j) = 0.0;
    }
    for (int j = 0; j < mat_size; j++) {
      Bd(i, j) = stan::math::uniform_rng(-5, 5, rng);
    }
  }
  stan::math::opencl_context.tuning_opts().tri_inverse_size_worth_transfer
      = MDIVIDE_OPENCL_OVERRIDE;
  Cv_cl = mdivide_left_tri<Eigen::Lower>(Av, Bd);
  Cv_cl(0, 0).grad();
  stan::math::opencl_context.tuning_opts().tri_inverse_size_worth_transfer
      = MDIVIDE_CPU_OVERRIDE;
  Cv = mdivide_left_tri<Eigen::Lower>(Av, Bd);
  Cv(0, 0).grad();
  EXPECT_MATRIX_NEAR(Cv, Cv_cl, 1.0E-12);
  EXPECT_MATRIX_NEAR(Cv.adj(), Cv_cl.adj(), 1.0E-12);
  stan::math::recover_memory();
  stan::math::opencl_context.tuning_opts().tri_inverse_size_worth_transfer
      = temp;
}

TEST(AgradRevMatrix, mdivide_left_tri_upper_grad_vv_cl) {
  int temp = stan::math::opencl_context.tuning_opts()
                 .tri_inverse_size_worth_transfer;

  using stan::math::matrix_d;
  using stan::math::matrix_v;
  using stan::math::mdivide_left_tri;
  using stan::math::multiply;
  int mat_size = 299;

  matrix_v Av(mat_size, mat_size);
  matrix_v Bv(mat_size, mat_size);
  matrix_v Cv(mat_size, mat_size);
  matrix_v Cv_cl(mat_size, mat_size);
  for (int i = 0; i < mat_size; i++) {
    for (int j = 0; j < i; j++) {
      Av(i, j) = 0.0;
    }
    Av(i, i) = 20.0;
    for (int j = i + 1; j < mat_size; j++) {
      Av(i, j) = stan::math::uniform_rng(-5, 5, rng);
    }
    for (int j = 0; j < mat_size; j++) {
      Bv(i, j) = stan::math::uniform_rng(-5, 5, rng);
    }
  }
  stan::math::opencl_context.tuning_opts().tri_inverse_size_worth_transfer
      = MDIVIDE_OPENCL_OVERRIDE;
  Cv_cl = mdivide_left_tri<Eigen::Upper>(Av, Bv);
  Cv_cl(0, 0).grad();
  stan::math::opencl_context.tuning_opts().tri_inverse_size_worth_transfer
      = MDIVIDE_CPU_OVERRIDE;
  Cv = mdivide_left_tri<Eigen::Upper>(Av, Bv);
  Cv(0, 0).grad();
  EXPECT_MATRIX_NEAR(Cv, Cv_cl, 1.0E-12);
  EXPECT_MATRIX_NEAR(Cv.adj(), Cv_cl.adj(), 1.0E-12);
  stan::math::recover_memory();
  stan::math::opencl_context.tuning_opts().tri_inverse_size_worth_transfer
      = temp;
}

TEST(AgradRevMatrix, mdivide_left_tri_upper_grad_dv_cl) {
  int temp = stan::math::opencl_context.tuning_opts()
                 .tri_inverse_size_worth_transfer;

  using stan::math::matrix_d;
  using stan::math::matrix_v;
  using stan::math::mdivide_left_tri;
  using stan::math::multiply;
  int mat_size = 301;

  matrix_d Ad(mat_size, mat_size);
  matrix_v Bv(mat_size, mat_size);
  matrix_v Cv(mat_size, mat_size);
  matrix_v Cv_cl(mat_size, mat_size);
  for (int i = 0; i < mat_size; i++) {
    for (int j = 0; j < i; j++) {
      Ad(i, j) = 0.0;
    }
    Ad(i, i) = 20.0;
    for (int j = i + 1; j < mat_size; j++) {
      Ad(i, j) = stan::math::uniform_rng(-5, 5, rng);
    }
    for (int j = 0; j < mat_size; j++) {
      Bv(i, j) = stan::math::uniform_rng(-5, 5, rng);
    }
  }
  stan::math::opencl_context.tuning_opts().tri_inverse_size_worth_transfer
      = MDIVIDE_OPENCL_OVERRIDE;
  Cv_cl = mdivide_left_tri<Eigen::Upper>(Ad, Bv);
  Cv_cl(0, 0).grad();
  stan::math::opencl_context.tuning_opts().tri_inverse_size_worth_transfer
      = MDIVIDE_CPU_OVERRIDE;
  Cv = mdivide_left_tri<Eigen::Upper>(Ad, Bv);
  Cv(0, 0).grad();
  EXPECT_MATRIX_NEAR(Cv, Cv_cl, 1.0E-12);
  EXPECT_MATRIX_NEAR(Cv.adj(), Cv_cl.adj(), 1.0E-12);
  stan::math::recover_memory();
  stan::math::opencl_context.tuning_opts().tri_inverse_size_worth_transfer
      = temp;
}

TEST(AgradRevMatrix, mdivide_left_tri_upper_grad_vd_cl) {
  int temp = stan::math::opencl_context.tuning_opts()
                 .tri_inverse_size_worth_transfer;

  using stan::math::matrix_d;
  using stan::math::matrix_v;
  using stan::math::mdivide_left_tri;
  using stan::math::multiply;
  int mat_size = 300;

  matrix_d Av(mat_size, mat_size);
  matrix_v Bd(mat_size, mat_size);
  matrix_v Cv(mat_size, mat_size);
  matrix_v Cv_cl(mat_size, mat_size);
  for (int i = 0; i < mat_size; i++) {
    for (int j = 0; j < i; j++) {
      Av(i, j) = 0.0;
    }
    Av(i, i) = 20.0;
    for (int j = i + 1; j < mat_size; j++) {
      Av(i, j) = stan::math::uniform_rng(-5, 5, rng);
    }
    for (int j = 0; j < mat_size; j++) {
      Bd(i, j) = stan::math::uniform_rng(-5, 5, rng);
    }
  }
  stan::math::opencl_context.tuning_opts().tri_inverse_size_worth_transfer
      = MDIVIDE_OPENCL_OVERRIDE;
  Cv_cl = mdivide_left_tri<Eigen::Upper>(Av, Bd);
  Cv_cl(0, 0).grad();
  stan::math::opencl_context.tuning_opts().tri_inverse_size_worth_transfer
      = MDIVIDE_CPU_OVERRIDE;
  Cv = mdivide_left_tri<Eigen::Upper>(Av, Bd);
  Cv(0, 0).grad();
  EXPECT_MATRIX_NEAR(Cv, Cv_cl, 1.0E-12);
  EXPECT_MATRIX_NEAR(Cv.adj(), Cv_cl.adj(), 1.0E-12);
  stan::math::recover_memory();
  stan::math::opencl_context.tuning_opts().tri_inverse_size_worth_transfer
      = temp;
}
#endif
