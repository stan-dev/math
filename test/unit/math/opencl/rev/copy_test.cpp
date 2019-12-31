#ifdef STAN_OPENCL
#include <stan/math/opencl/rev/opencl.hpp>
#include <gtest/gtest.h>
#include <CL/cl2.hpp>
#include <algorithm>
#include <vector>

TEST(MathMatrixRevCL, matrix_cl_vector_copy) {
  using stan::math::var;
  stan::math::vector_v d1_cpu(3);
  stan::math::vector_d d1_a_cpu(3);
  stan::math::vector_d d1_b_cpu(3);
  d1_cpu << 1, 2, 3;
  // vector
  stan::math::matrix_cl<var> d11_cl(3, 1);
  stan::math::matrix_cl<var> d111_cl(3, 1);
  EXPECT_NO_THROW(d11_cl = stan::math::to_matrix_cl(d1_cpu));
  EXPECT_NO_THROW(d11_cl.adj() = stan::math::copy_cl(d11_cl.val()));
  EXPECT_NO_THROW(d111_cl = stan::math::copy_cl(d11_cl));
  EXPECT_NO_THROW(d1_a_cpu = stan::math::from_matrix_cl(d11_cl));
  EXPECT_NO_THROW(d1_b_cpu = stan::math::from_matrix_cl(d111_cl));
  EXPECT_EQ(1, d1_a_cpu(0));
  EXPECT_EQ(2, d1_a_cpu(1));
  EXPECT_EQ(3, d1_a_cpu(2));
  EXPECT_EQ(1, d1_b_cpu(0));
  EXPECT_EQ(2, d1_b_cpu(1));
  EXPECT_EQ(3, d1_b_cpu(2));
}

TEST(MathMatrixRevCL, matrix_cl_matrix_copy) {
  using stan::math::var;
  stan::math::matrix_v d2_cpu(2, 3);
  stan::math::matrix_d d2_a_cpu(2, 3);
  stan::math::matrix_d d2_b_cpu(2, 3);
  stan::math::matrix_v d0_cpu(0, 0);
  d2_cpu << 1, 2, 3, 4, 5, 6;
  // matrix
  stan::math::matrix_cl<var> d00_cl;
  stan::math::matrix_cl<var> d000_cl;
  stan::math::matrix_cl<var> d22_cl(2, 3);
  stan::math::matrix_cl<var> d222_cl(2, 3);
  EXPECT_NO_THROW(d22_cl = stan::math::to_matrix_cl(d2_cpu));
  EXPECT_NO_THROW(d222_cl = stan::math::copy_cl(d22_cl));
  EXPECT_NO_THROW(d2_a_cpu = stan::math::from_matrix_cl(d22_cl.val()));
  EXPECT_NO_THROW(d2_b_cpu = stan::math::from_matrix_cl(d222_cl.val()));
  EXPECT_EQ(1, d2_a_cpu(0, 0));
  EXPECT_EQ(2, d2_a_cpu(0, 1));
  EXPECT_EQ(3, d2_a_cpu(0, 2));
  EXPECT_EQ(4, d2_a_cpu(1, 0));
  EXPECT_EQ(5, d2_a_cpu(1, 1));
  EXPECT_EQ(6, d2_a_cpu(1, 2));
  EXPECT_EQ(1, d2_b_cpu(0, 0));
  EXPECT_EQ(2, d2_b_cpu(0, 1));
  EXPECT_EQ(3, d2_b_cpu(0, 2));
  EXPECT_EQ(4, d2_b_cpu(1, 0));
  EXPECT_EQ(5, d2_b_cpu(1, 1));
  EXPECT_EQ(6, d2_b_cpu(1, 2));
  d00_cl = stan::math::to_matrix_cl(d0_cpu);
  EXPECT_NO_THROW(d00_cl = stan::math::to_matrix_cl(d0_cpu));
  EXPECT_NO_THROW(d0_cpu = stan::math::from_matrix_cl(d00_cl.val()));
  EXPECT_NO_THROW(d000_cl = stan::math::copy_cl(d00_cl));
}

TEST(MathMatrixRevCL, matrix_cl_pack_unpack_copy_lower) {
  using stan::math::var;
  using stan::math::vari;
  int size = 42;
  int packed_size = size * (size + 1) / 2;
  vari** packed_mat
      = stan::math::ChainableStack::instance_->memalloc_.alloc_array<vari*>(
          packed_size);
  std::vector<double> packed_mat_dst(packed_size);
  for (size_t i = 0; i < packed_size; i++) {
    packed_mat[i] = new stan::math::vari(i);
  }
  stan::math::matrix_d m_flat_cpu(size, size);
  auto m_cl = stan::math::packed_copy<stan::math::matrix_cl_view::Lower>(
      packed_mat, size);
  m_flat_cpu = stan::math::from_matrix_cl(m_cl.val());
  size_t pos = 0;
  for (size_t j = 0; j < size; ++j) {
    for (size_t i = 0; i < j; i++) {
      EXPECT_EQ(m_flat_cpu(i, j), 0.0);
    }
    for (size_t i = j; i < size; ++i) {
      EXPECT_EQ(m_flat_cpu(i, j), packed_mat[pos]->val_);
      pos++;
    }
  }
  packed_mat_dst = stan::math::packed_copy(m_cl);
  for (size_t i = 0; i < packed_size; i++) {
    EXPECT_EQ(packed_mat[i]->adj_, packed_mat_dst[i]);
  }
  packed_mat_dst = stan::math::packed_copy(m_cl.val());
  for (size_t i = 0; i < packed_size; i++) {
    EXPECT_EQ(packed_mat[i]->val_, packed_mat_dst[i]);
  }
}

TEST(MathMatrixRevCL, matrix_cl_pack_unpack_copy_upper) {
  using stan::math::var;
  using stan::math::vari;
  int size = 51;
  int packed_size = size * (size + 1) / 2;
  vari** packed_mat
      = stan::math::ChainableStack::instance_->memalloc_.alloc_array<vari*>(
          packed_size);
  std::vector<double> packed_mat_dst(packed_size);
  for (size_t i = 0; i < packed_size; i++) {
    packed_mat[i] = new stan::math::vari(i);
  }
  stan::math::matrix_d m_flat_cpu(size, size);
  auto m_cl = stan::math::packed_copy<stan::math::matrix_cl_view::Upper>(
      packed_mat, size);
  m_flat_cpu = stan::math::from_matrix_cl(m_cl.val());
  size_t pos = 0;
  for (size_t j = 0; j < size; ++j) {
    for (size_t i = 0; i <= j; i++) {
      EXPECT_EQ(m_flat_cpu(i, j), packed_mat[pos]->val_);
      pos++;
    }
    for (size_t i = j + 1; i < size; ++i) {
      EXPECT_EQ(m_flat_cpu(i, j), 0.0);
    }
  }
  packed_mat_dst = stan::math::packed_copy(m_cl);
  for (size_t i = 0; i < packed_size; i++) {
    EXPECT_EQ(packed_mat[i]->adj_, packed_mat_dst[i]);
  }
  packed_mat_dst = stan::math::packed_copy(m_cl.val());
  for (size_t i = 0; i < packed_size; i++) {
    EXPECT_EQ(packed_mat[i]->val_, packed_mat_dst[i]);
  }
}

TEST(MathMatrixRevCL, matrix_cl_pack_unpack_copy_exception) {
  using stan::math::var;
  using stan::math::vari;
  vari** packed_mat;
  stan::math::matrix_cl<var> m_cl_zero;
  EXPECT_NO_THROW(stan::math::packed_copy<stan::math::matrix_cl_view::Upper>(
      packed_mat, 0));
  m_cl_zero.view(stan::math::matrix_cl_view::Upper);
  EXPECT_NO_THROW(stan::math::packed_copy(m_cl_zero));
  m_cl_zero.view(stan::math::matrix_cl_view::Entire);
  EXPECT_THROW(stan::math::packed_copy(m_cl_zero), std::invalid_argument);
  /* With vari** we can't check to throw here.
  EXPECT_THROW(
      stan::math::packed_copy<stan::math::matrix_cl_view::Upper>(packed_mat, 1),
      std::invalid_argument);
  */
}
#endif
