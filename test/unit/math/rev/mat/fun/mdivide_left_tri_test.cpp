#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/rev/mat/util.hpp>

#ifdef STAN_OPENCL
#include <boost/random/mersenne_twister.hpp>
#endif

TEST(AgradRevMatrix, mdivide_left_tri_val) {
  using stan::math::matrix_d;
  using stan::math::matrix_v;
  using stan::math::mdivide_left_tri;

  matrix_v Av(2, 2);
  matrix_v Av_inv(2, 2);
  matrix_d Ad(2, 2);
  matrix_v I;

  Av << 2.0, 0.0, 5.0, 7.0;
  Ad << 2.0, 0.0, 5.0, 7.0;

  I = mdivide_left_tri<Eigen::Lower>(Av, Av);
  EXPECT_NEAR(1.0, I(0, 0).val(), 1.0E-12);
  EXPECT_NEAR(0.0, I(0, 1).val(), 1.0E-12);
  EXPECT_NEAR(0.0, I(1, 0).val(), 1.0E-12);
  EXPECT_NEAR(1.0, I(1, 1).val(), 1.0e-12);

  I = mdivide_left_tri<Eigen::Lower>(Av, Ad);
  EXPECT_NEAR(1.0, I(0, 0).val(), 1.0E-12);
  EXPECT_NEAR(0.0, I(0, 1).val(), 1.0E-12);
  EXPECT_NEAR(0.0, I(1, 0).val(), 1.0E-12);
  EXPECT_NEAR(1.0, I(1, 1).val(), 1.0e-12);

  I = mdivide_left_tri<Eigen::Lower>(Ad, Av);
  EXPECT_NEAR(1.0, I(0, 0).val(), 1.0E-12);
  EXPECT_NEAR(0.0, I(0, 1).val(), 1.0E-12);
  EXPECT_NEAR(0.0, I(1, 0).val(), 1.0E-12);
  EXPECT_NEAR(1.0, I(1, 1).val(), 1.0e-12);

  Av_inv = mdivide_left_tri<Eigen::Lower>(Av);
  I = Av * Av_inv;
  EXPECT_NEAR(1.0, I(0, 0).val(), 1.0E-12);
  EXPECT_NEAR(0.0, I(0, 1).val(), 1.0E-12);
  EXPECT_NEAR(0.0, I(1, 0).val(), 1.0E-12);
  EXPECT_NEAR(1.0, I(1, 1).val(), 1.0e-12);

  Av << 2.0, 3.0, 0.0, 7.0;
  Ad << 2.0, 3.0, 0.0, 7.0;

  I = mdivide_left_tri<Eigen::Upper>(Av, Av);
  EXPECT_NEAR(1.0, I(0, 0).val(), 1.0E-12);
  EXPECT_NEAR(0.0, I(0, 1).val(), 1.0E-12);
  EXPECT_NEAR(0.0, I(1, 0).val(), 1.0E-12);
  EXPECT_NEAR(1.0, I(1, 1).val(), 1.0e-12);

  I = mdivide_left_tri<Eigen::Upper>(Av, Ad);
  EXPECT_NEAR(1.0, I(0, 0).val(), 1.0E-12);
  EXPECT_NEAR(0.0, I(0, 1).val(), 1.0E-12);
  EXPECT_NEAR(0.0, I(1, 0).val(), 1.0E-12);
  EXPECT_NEAR(1.0, I(1, 1).val(), 1.0e-12);

  I = mdivide_left_tri<Eigen::Upper>(Ad, Av);
  EXPECT_NEAR(1.0, I(0, 0).val(), 1.0E-12);
  EXPECT_NEAR(0.0, I(0, 1).val(), 1.0E-12);
  EXPECT_NEAR(0.0, I(1, 0).val(), 1.0E-12);
  EXPECT_NEAR(1.0, I(1, 1).val(), 1.0e-12);
}

TEST(AgradRevMatrix, mdivide_left_tri_lower_grad_vv) {
  using stan::math::matrix_d;
  using stan::math::matrix_v;
  using stan::math::mdivide_left_tri;
  using stan::math::multiply;

  matrix_d Ad(2, 2), Ad_tmp(2, 2);
  matrix_d Bd(2, 2), Bd_tmp(2, 2);
  matrix_d Cd(2, 2);

  Ad << 2.0, 0.0, 5.0, 7.0;
  Bd << 12.0, 13.0, 15.0, 17.0;

  size_type i, j, k, l;
  for (i = 0; i < Bd.rows(); i++) {
    for (j = 0; j < Bd.cols(); j++) {
      matrix_v A(2, 2);
      matrix_v B(2, 2);
      matrix_v C;
      VEC g;

      for (k = 0; k < 4; k++) {
        A(k) = Ad(k);
        B(k) = Bd(k);
      }

      C = mdivide_left_tri<Eigen::Lower>(A, B);
      AVEC x = createAVEC(A(0, 0), A(1, 0), A(0, 1), A(1, 1), B(0, 0), B(1, 0),
                          B(0, 1), B(1, 1));
      C(i, j).grad(x, g);

      for (l = 0; l < Ad_tmp.cols(); l++) {
        for (k = 0; k < Ad_tmp.rows(); k++) {
          if (k >= l) {
            Ad_tmp.setZero();
            Ad_tmp(k, l) = 1.0;
            Cd = -mdivide_left_tri<Eigen::Lower>(
                Ad, multiply(Ad_tmp, mdivide_left_tri<Eigen::Lower>(Ad, Bd)));
            EXPECT_NEAR(Cd(i, j), g[k + l * Ad_tmp.rows()], 1.0E-12);
          } else {
            EXPECT_NEAR(0.0, g[k + l * Ad_tmp.rows()], 1.0E-12);
          }
        }
      }
      for (k = 0; k < 4; k++) {
        Bd_tmp.setZero();
        Bd_tmp(k) = 1.0;
        Cd = mdivide_left_tri<Eigen::Lower>(Ad, Bd_tmp);
        EXPECT_NEAR(Cd(i, j), g[4 + k], 1.0E-12);
      }
    }
  }
}

TEST(AgradRevMatrix, mdivide_left_tri_lower_grad_dv) {
  using stan::math::matrix_d;
  using stan::math::matrix_v;
  using stan::math::mdivide_left_tri;
  using stan::math::multiply;

  matrix_d Ad(2, 2), Ad_tmp(2, 2);
  matrix_d Bd(2, 2), Bd_tmp(2, 2);
  matrix_d Cd(2, 2);

  Ad << 2.0, 0.0, 5.0, 7.0;
  Bd << 12.0, 13.0, 15.0, 17.0;

  size_type i, j, k;
  for (i = 0; i < Bd.rows(); i++) {
    for (j = 0; j < Bd.cols(); j++) {
      matrix_v B(2, 2);
      matrix_v C;
      VEC g;

      for (k = 0; k < 4; k++) {
        B(k) = Bd(k);
      }

      C = mdivide_left_tri<Eigen::Lower>(Ad, B);
      AVEC x = createAVEC(B(0, 0), B(1, 0), B(0, 1), B(1, 1));
      C(i, j).grad(x, g);

      for (k = 0; k < 4; k++) {
        Bd_tmp.setZero();
        Bd_tmp(k) = 1.0;
        Cd = mdivide_left_tri<Eigen::Lower>(Ad, Bd_tmp);
        EXPECT_NEAR(Cd(i, j), g[k], 1.0E-12);
      }
    }
  }
}

TEST(AgradRevMatrix, mdivide_left_tri_lower_grad_vd) {
  using stan::math::matrix_d;
  using stan::math::matrix_v;
  using stan::math::mdivide_left_tri;
  using stan::math::multiply;

  matrix_d Ad(2, 2), Ad_tmp(2, 2);
  matrix_d Bd(2, 2), Bd_tmp(2, 2);
  matrix_d Cd(2, 2);

  Ad << 2.0, 0.0, 5.0, 7.0;
  Bd << 12.0, 13.0, 15.0, 17.0;

  size_type i, j, k, l;
  for (i = 0; i < Bd.rows(); i++) {
    for (j = 0; j < Bd.cols(); j++) {
      matrix_v A(2, 2);
      matrix_v C;
      VEC g;

      for (k = 0; k < 4; k++) {
        A(k) = Ad(k);
      }

      C = mdivide_left_tri<Eigen::Lower>(A, Bd);
      AVEC x = createAVEC(A(0, 0), A(1, 0), A(0, 1), A(1, 1));
      C(i, j).grad(x, g);

      for (l = 0; l < Ad_tmp.cols(); l++) {
        for (k = 0; k < Ad_tmp.rows(); k++) {
          if (k >= l) {
            Ad_tmp.setZero();
            Ad_tmp(k, l) = 1.0;
            Cd = -mdivide_left_tri<Eigen::Lower>(
                Ad, multiply(Ad_tmp, mdivide_left_tri<Eigen::Lower>(Ad, Bd)));
            EXPECT_NEAR(Cd(i, j), g[k + l * Ad_tmp.rows()], 1.0E-12);
          } else {
            EXPECT_NEAR(0.0, g[k + l * Ad_tmp.rows()], 1.0E-12);
          }
        }
      }
    }
  }
}
TEST(AgradRevMatrix, mdivide_left_tri_upper_grad_vv) {
  using stan::math::matrix_d;
  using stan::math::matrix_v;
  using stan::math::mdivide_left_tri;
  using stan::math::multiply;

  matrix_d Ad(2, 2), Ad_tmp(2, 2);
  matrix_d Bd(2, 2), Bd_tmp(2, 2);
  matrix_d Cd(2, 2);

  Ad << 2.0, 3.0, 0.0, 7.0;
  Bd << 12.0, 13.0, 15.0, 17.0;

  size_type i, j, k, l;
  for (i = 0; i < Bd.rows(); i++) {
    for (j = 0; j < Bd.cols(); j++) {
      matrix_v A(2, 2);
      matrix_v B(2, 2);
      matrix_v C;
      VEC g;

      for (k = 0; k < 4; k++) {
        A(k) = Ad(k);
        B(k) = Bd(k);
      }

      C = mdivide_left_tri<Eigen::Upper>(A, B);
      AVEC x = createAVEC(A(0, 0), A(1, 0), A(0, 1), A(1, 1), B(0, 0), B(1, 0),
                          B(0, 1), B(1, 1));
      C(i, j).grad(x, g);

      for (l = 0; l < Ad_tmp.cols(); l++) {
        for (k = 0; k < Ad_tmp.rows(); k++) {
          if (k <= l) {
            Ad_tmp.setZero();
            Ad_tmp(k, l) = 1.0;
            Cd = -mdivide_left_tri<Eigen::Upper>(
                Ad, multiply(Ad_tmp, mdivide_left_tri<Eigen::Upper>(Ad, Bd)));
            EXPECT_NEAR(Cd(i, j), g[k + l * Ad_tmp.rows()], 1.0E-12);
          } else {
            EXPECT_NEAR(0.0, g[k + l * Ad_tmp.rows()], 1.0E-12);
          }
        }
      }
      for (k = 0; k < 4; k++) {
        Bd_tmp.setZero();
        Bd_tmp(k) = 1.0;
        Cd = mdivide_left_tri<Eigen::Upper>(Ad, Bd_tmp);
        EXPECT_NEAR(Cd(i, j), g[4 + k], 1.0E-12);
      }
    }
  }
}

TEST(AgradRevMatrix, mdivide_left_tri_upper_grad_dv) {
  using stan::math::matrix_d;
  using stan::math::matrix_v;
  using stan::math::mdivide_left_tri;
  using stan::math::multiply;

  matrix_d Ad(2, 2), Ad_tmp(2, 2);
  matrix_d Bd(2, 2), Bd_tmp(2, 2);
  matrix_d Cd(2, 2);

  Ad << 2.0, 3.0, 0.0, 7.0;
  Bd << 12.0, 13.0, 15.0, 17.0;

  size_type i, j, k;
  for (i = 0; i < Bd.rows(); i++) {
    for (j = 0; j < Bd.cols(); j++) {
      matrix_v B(2, 2);
      matrix_v C;
      VEC g;

      for (k = 0; k < 4; k++) {
        B(k) = Bd(k);
      }

      C = mdivide_left_tri<Eigen::Upper>(Ad, B);
      AVEC x = createAVEC(B(0, 0), B(1, 0), B(0, 1), B(1, 1));
      C(i, j).grad(x, g);

      for (k = 0; k < 4; k++) {
        Bd_tmp.setZero();
        Bd_tmp(k) = 1.0;
        Cd = mdivide_left_tri<Eigen::Upper>(Ad, Bd_tmp);
        EXPECT_NEAR(Cd(i, j), g[k], 1.0E-12);
      }
    }
  }
}

TEST(AgradRevMatrix, mdivide_left_tri_upper_grad_vd) {
  using stan::math::matrix_d;
  using stan::math::matrix_v;
  using stan::math::mdivide_left_tri;
  using stan::math::multiply;

  matrix_d Ad(2, 2), Ad_tmp(2, 2);
  matrix_d Bd(2, 2), Bd_tmp(2, 2);
  matrix_d Cd(2, 2);

  Ad << 2.0, 3.0, 0.0, 7.0;
  Bd << 12.0, 13.0, 15.0, 17.0;

  size_type i, j, k, l;
  for (i = 0; i < Bd.rows(); i++) {
    for (j = 0; j < Bd.cols(); j++) {
      matrix_v A(2, 2);
      matrix_v C;
      VEC g;

      for (k = 0; k < 4; k++) {
        A(k) = Ad(k);
      }

      C = mdivide_left_tri<Eigen::Upper>(A, Bd);
      AVEC x = createAVEC(A(0, 0), A(1, 0), A(0, 1), A(1, 1));
      C(i, j).grad(x, g);

      for (l = 0; l < Ad_tmp.cols(); l++) {
        for (k = 0; k < Ad_tmp.rows(); k++) {
          if (k <= l) {
            Ad_tmp.setZero();
            Ad_tmp(k, l) = 1.0;
            Cd = -mdivide_left_tri<Eigen::Upper>(
                Ad, multiply(Ad_tmp, mdivide_left_tri<Eigen::Upper>(Ad, Bd)));
            EXPECT_NEAR(Cd(i, j), g[k + l * Ad_tmp.rows()], 1.0E-12);
          } else {
            EXPECT_NEAR(0.0, g[k + l * Ad_tmp.rows()], 1.0E-12);
          }
        }
      }
    }
  }
}

TEST(AgradRevMatrix, check_varis_on_stack) {
  using stan::math::value_of;
  stan::math::matrix_v A(2, 2);
  A << 2.0, 0.0, 5.0, 7.0;

  test::check_varis_on_stack(stan::math::mdivide_left_tri<Eigen::Lower>(A, A));
  test::check_varis_on_stack(
      stan::math::mdivide_left_tri<Eigen::Lower>(A, value_of(A)));
  test::check_varis_on_stack(
      stan::math::mdivide_left_tri<Eigen::Lower>(value_of(A), A));
}

#ifdef STAN_OPENCL
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
  int size = 111;
  matrix_v Av(size, size);
  matrix_d Ad(size, size);
  matrix_v I, I_cl;
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < i; j++) {
      Ad(i, j) = stan::math::uniform_rng(-5, 5, rng);
      Av(i, j) = Ad(i, j);
    }
    Ad(i, i) = 20.0;
    Av(i, i) = Ad(i, i);
    for (int j = i + 1; j < size; j++) {
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

  for (int i = 0; i < size; i++) {
    for (int j = 0; j < i; j++) {
      Ad(i, j) = 0.0;
      Av(i, j) = Ad(i, j);
    }
    Ad(i, i) = 20.0;
    Av(i, i) = Ad(i, i);
    for (int j = i + 1; j < size; j++) {
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
  int size = 256;

  matrix_v Av(size, size);
  matrix_v Bv(size, size);
  matrix_v Cv(size, size);
  matrix_v Cv_cl(size, size);
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < i; j++) {
      Av(i, j) = stan::math::uniform_rng(-5, 5, rng);
    }
    Av(i, i) = 20.0;
    for (int j = i + 1; j < size; j++) {
      Av(i, j) = 0.0;
    }
    for (int j = 0; j < size; j++) {
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
  int size = 213;

  matrix_d Ad(size, size);
  matrix_v Bv(size, size);
  matrix_v Cv(size, size);
  matrix_v Cv_cl(size, size);
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < i; j++) {
      Ad(i, j) = stan::math::uniform_rng(-5, 5, rng);
    }
    Ad(i, i) = 20.0;
    for (int j = i + 1; j < size; j++) {
      Ad(i, j) = 0.0;
    }
    for (int j = 0; j < size; j++) {
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
  int size = 250;

  matrix_d Av(size, size);
  matrix_v Bd(size, size);
  matrix_v Cv(size, size);
  matrix_v Cv_cl(size, size);
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < i; j++) {
      Av(i, j) = stan::math::uniform_rng(-5, 5, rng);
    }
    Av(i, i) = 20.0;
    for (int j = i + 1; j < size; j++) {
      Av(i, j) = 0.0;
    }
    for (int j = 0; j < size; j++) {
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
  int size = 299;

  matrix_v Av(size, size);
  matrix_v Bv(size, size);
  matrix_v Cv(size, size);
  matrix_v Cv_cl(size, size);
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < i; j++) {
      Av(i, j) = 0.0;
    }
    Av(i, i) = 20.0;
    for (int j = i + 1; j < size; j++) {
      Av(i, j) = stan::math::uniform_rng(-5, 5, rng);
    }
    for (int j = 0; j < size; j++) {
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
  int size = 301;

  matrix_d Ad(size, size);
  matrix_v Bv(size, size);
  matrix_v Cv(size, size);
  matrix_v Cv_cl(size, size);
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < i; j++) {
      Ad(i, j) = 0.0;
    }
    Ad(i, i) = 20.0;
    for (int j = i + 1; j < size; j++) {
      Ad(i, j) = stan::math::uniform_rng(-5, 5, rng);
    }
    for (int j = 0; j < size; j++) {
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
  int size = 300;

  matrix_d Av(size, size);
  matrix_v Bd(size, size);
  matrix_v Cv(size, size);
  matrix_v Cv_cl(size, size);
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < i; j++) {
      Av(i, j) = 0.0;
    }
    Av(i, i) = 20.0;
    for (int j = i + 1; j < size; j++) {
      Av(i, j) = stan::math::uniform_rng(-5, 5, rng);
    }
    for (int j = 0; j < size; j++) {
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
