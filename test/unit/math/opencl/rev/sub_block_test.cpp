#ifdef STAN_OPENCL
#include <stan/math/prim.hpp>
#include <stan/math/opencl/rev/copy.hpp>
#include <stan/math/opencl/rev/matrix_cl.hpp>
#include <stan/math/opencl/rev/sub_block.hpp>
#include <stan/math/opencl/zeros.hpp>
#include <gtest/gtest.h>
#include <algorithm>

TEST(MathMatrixRevCL, sub_block_pass_vari) {
  using stan::math::matrix_cl;
  using stan::math::matrix_d;
  using stan::math::matrix_v;
  using stan::math::matrix_vi;
  using stan::math::var;
  using stan::math::vari;
  vari** d1_vals(
      stan::math::ChainableStack::instance_->memalloc_.alloc_array<vari*>(9));
  vari** d2_vals(
      stan::math::ChainableStack::instance_->memalloc_.alloc_array<vari*>(16));

  for (int i = 0; i < 9; i++) {
    d1_vals[i] = new vari(i);
  }
  for (int i = 0; i < 16; i++) {
    d2_vals[i] = new vari(15 - i);
  }
  const matrix_vi d1 = Eigen::Map<matrix_vi>(d1_vals, 3, 3);
  const matrix_vi d2 = Eigen::Map<matrix_vi>(d2_vals, 4, 4);
  matrix_cl<var> d11(d1);
  matrix_cl<var> d22(d2);
  d22.sub_block(d11, 0, 0, 0, 0, 2, 2);
  matrix_d d3 = stan::math::from_matrix_cl(d22.val());
  EXPECT_EQ(0, d3(0, 0));
  EXPECT_EQ(3, d3(0, 1));
  EXPECT_EQ(1, d3(1, 0));
  EXPECT_EQ(4, d3(1, 1));
  stan::math::recover_memory();
}

TEST(MathMatrixRevCL, sub_block_exception) {
  using stan::math::var;
  using stan::math::vari;
  stan::math::matrix_v d1(3, 3);
  stan::math::matrix_v d2(4, 4);
  for (int i = 0; i < 9; i++) {
    d1(i) = new vari(i);
  }
  for (int i = 0; i < 16; i++) {
    d2(i) = new vari(15 - i);
  }

  stan::math::matrix_cl<var> d11(d1);
  stan::math::matrix_cl<var> d22(d2);
  EXPECT_THROW(d22.sub_block(d11, 1, 1, 0, 0, 4, 4), std::domain_error);
  EXPECT_THROW(d22.sub_block(d11, 4, 4, 0, 0, 2, 2), std::domain_error);
}

TEST(MathMatrixRevCL, sub_block_triangular) {
  using stan::math::var;
  Eigen::Matrix<var, -1, -1> a = Eigen::Matrix<var, -1, -1>::Zero(3, 3);
  Eigen::Matrix<var, -1, -1> b = Eigen::Matrix<var, -1, -1>::Ones(3, 3);
  stan::math::matrix_cl<var> a_cl(3, 3);
  stan::math::matrix_cl<var> b_cl(b);
  Eigen::Matrix<var, -1, -1> c;

  a_cl = stan::math::to_matrix_cl(a);
  a_cl.view(stan::math::matrix_cl_view::Lower);
  b_cl.view(stan::math::matrix_cl_view::Lower);
  a_cl.sub_block(b_cl, 0, 1, 0, 1, 2, 2);
  EXPECT_EQ(a_cl.val().view(), stan::math::matrix_cl_view::Lower);
  c = stan::math::from_matrix_cl(a_cl.val());
  EXPECT_EQ(c(0, 0), 0);
  EXPECT_EQ(c(0, 1), 0);
  EXPECT_EQ(c(0, 2), 0);
  EXPECT_EQ(c(1, 0), 0);
  EXPECT_EQ(c(1, 1), 1);
  EXPECT_EQ(c(1, 2), 0);
  EXPECT_EQ(c(2, 0), 0);
  EXPECT_EQ(c(2, 1), 0);
  EXPECT_EQ(c(2, 2), 0);

  a_cl = stan::math::to_matrix_cl(a);
  a_cl.view(stan::math::matrix_cl_view::Lower);
  b_cl.view(stan::math::matrix_cl_view::Lower);
  a_cl.sub_block(b_cl, 1, 0, 1, 1, 2, 2);
  EXPECT_EQ(a_cl.val().view(), stan::math::matrix_cl_view::Entire);
  c = stan::math::from_matrix_cl(a_cl.val());
  EXPECT_EQ(c(0, 0), 0);
  EXPECT_EQ(c(0, 1), 0);
  EXPECT_EQ(c(0, 2), 0);
  EXPECT_EQ(c(1, 0), 0);
  EXPECT_EQ(c(1, 1), 1);
  EXPECT_EQ(c(1, 2), 1);
  EXPECT_EQ(c(2, 0), 0);
  EXPECT_EQ(c(2, 1), 1);
  EXPECT_EQ(c(2, 2), 1);

  a_cl = stan::math::to_matrix_cl(a);
  a_cl.view(stan::math::matrix_cl_view::Lower);
  b_cl.view(stan::math::matrix_cl_view::Upper);
  a_cl.sub_block(b_cl, 0, 0, 1, 0, 2, 2);
  EXPECT_EQ(a_cl.val().view(), stan::math::matrix_cl_view::Lower);
  c = stan::math::from_matrix_cl(a_cl.val());
  EXPECT_EQ(c(0, 0), 0);
  EXPECT_EQ(c(0, 1), 0);
  EXPECT_EQ(c(0, 2), 0);
  EXPECT_EQ(c(1, 0), 1);
  EXPECT_EQ(c(1, 1), 1);
  EXPECT_EQ(c(1, 2), 0);
  EXPECT_EQ(c(2, 0), 0);
  EXPECT_EQ(c(2, 1), 1);
  EXPECT_EQ(c(2, 2), 0);

  a_cl = stan::math::to_matrix_cl(a);
  a_cl.view(stan::math::matrix_cl_view::Lower);
  b_cl.view(stan::math::matrix_cl_view::Upper);
  a_cl.sub_block(b_cl, 0, 0, 0, 0, 2, 2);
  EXPECT_EQ(a_cl.val().view(), stan::math::matrix_cl_view::Entire);
  c = stan::math::from_matrix_cl(a_cl.val());
  EXPECT_EQ(c(0, 0), 1);
  EXPECT_EQ(c(0, 1), 1);
  EXPECT_EQ(c(0, 2), 0);
  EXPECT_EQ(c(1, 0), 0);
  EXPECT_EQ(c(1, 1), 1);
  EXPECT_EQ(c(1, 2), 0);
  EXPECT_EQ(c(2, 0), 0);
  EXPECT_EQ(c(2, 1), 0);
  EXPECT_EQ(c(2, 2), 0);

  a_cl = stan::math::to_matrix_cl(a);
  a_cl.view(stan::math::matrix_cl_view::Lower);
  b_cl.view(stan::math::matrix_cl_view::Upper);
  a_cl.sub_block(b_cl, 1, 0, 1, 0, 2, 2);
  EXPECT_EQ(a_cl.val().view(), stan::math::matrix_cl_view::Diagonal);
  c = stan::math::from_matrix_cl(a_cl.val());
  EXPECT_EQ(c(0, 0), 0);
  EXPECT_EQ(c(0, 1), 0);
  EXPECT_EQ(c(0, 2), 0);
  EXPECT_EQ(c(1, 0), 0);
  EXPECT_EQ(c(1, 1), 1);
  EXPECT_EQ(c(1, 2), 0);
  EXPECT_EQ(c(2, 0), 0);
  EXPECT_EQ(c(2, 1), 0);
  EXPECT_EQ(c(2, 2), 0);

  a_cl = stan::math::to_matrix_cl(a);
  a_cl.view(stan::math::matrix_cl_view::Lower);
  b_cl.view(stan::math::matrix_cl_view::Upper);
  a_cl.sub_block(b_cl, 1, 0, 1, 0, 2, 3);
  EXPECT_EQ(a_cl.val().view(), stan::math::matrix_cl_view::Upper);
  c = stan::math::from_matrix_cl(a_cl.val());
  EXPECT_EQ(c(0, 0), 0);
  EXPECT_EQ(c(0, 1), 0);
  EXPECT_EQ(c(0, 2), 0);
  EXPECT_EQ(c(1, 0), 0);
  EXPECT_EQ(c(1, 1), 1);
  EXPECT_EQ(c(1, 2), 1);
  EXPECT_EQ(c(2, 0), 0);
  EXPECT_EQ(c(2, 1), 0);
  EXPECT_EQ(c(2, 2), 1);

  a_cl = stan::math::to_matrix_cl(a);
  a_cl.view(stan::math::matrix_cl_view::Upper);
  b_cl.view(stan::math::matrix_cl_view::Upper);
  a_cl.sub_block(b_cl, 1, 0, 1, 0, 2, 2);
  EXPECT_EQ(a_cl.val().view(), stan::math::matrix_cl_view::Upper);
  c = stan::math::from_matrix_cl(a_cl.val());
  EXPECT_EQ(c(0, 0), 0);
  EXPECT_EQ(c(0, 1), 0);
  EXPECT_EQ(c(0, 2), 0);
  EXPECT_EQ(c(1, 0), 0);
  EXPECT_EQ(c(1, 1), 1);
  EXPECT_EQ(c(1, 2), 0);
  EXPECT_EQ(c(2, 0), 0);
  EXPECT_EQ(c(2, 1), 0);
  EXPECT_EQ(c(2, 2), 0);

  a_cl = stan::math::to_matrix_cl(a);
  a_cl.view(stan::math::matrix_cl_view::Upper);
  b_cl.view(stan::math::matrix_cl_view::Upper);
  a_cl.sub_block(b_cl, 0, 1, 1, 1, 2, 2);
  EXPECT_EQ(a_cl.val().view(), stan::math::matrix_cl_view::Entire);
  c = stan::math::from_matrix_cl(a_cl.val());
  EXPECT_EQ(c(0, 0), 0);
  EXPECT_EQ(c(0, 1), 0);
  EXPECT_EQ(c(0, 2), 0);
  EXPECT_EQ(c(1, 0), 0);
  EXPECT_EQ(c(1, 1), 1);
  EXPECT_EQ(c(1, 2), 1);
  EXPECT_EQ(c(2, 0), 0);
  EXPECT_EQ(c(2, 1), 1);
  EXPECT_EQ(c(2, 2), 1);

  a_cl = stan::math::to_matrix_cl(a);
  a_cl.view(stan::math::matrix_cl_view::Upper);
  b_cl.view(stan::math::matrix_cl_view::Lower);
  a_cl.sub_block(b_cl, 0, 0, 0, 1, 2, 2);
  EXPECT_EQ(a_cl.val().view(), stan::math::matrix_cl_view::Upper);
  c = stan::math::from_matrix_cl(a_cl.val());
  EXPECT_EQ(c(0, 0), 0);
  EXPECT_EQ(c(0, 1), 1);
  EXPECT_EQ(c(0, 2), 0);
  EXPECT_EQ(c(1, 0), 0);
  EXPECT_EQ(c(1, 1), 1);
  EXPECT_EQ(c(1, 2), 1);
  EXPECT_EQ(c(2, 0), 0);
  EXPECT_EQ(c(2, 1), 0);
  EXPECT_EQ(c(2, 2), 0);

  a_cl = stan::math::to_matrix_cl(a);
  a_cl.view(stan::math::matrix_cl_view::Upper);
  b_cl.view(stan::math::matrix_cl_view::Lower);
  a_cl.sub_block(b_cl, 0, 0, 0, 0, 2, 2);
  EXPECT_EQ(a_cl.val().view(), stan::math::matrix_cl_view::Entire);
  c = stan::math::from_matrix_cl(a_cl.val());
  EXPECT_EQ(c(0, 0), 1);
  EXPECT_EQ(c(0, 1), 0);
  EXPECT_EQ(c(0, 2), 0);
  EXPECT_EQ(c(1, 0), 1);
  EXPECT_EQ(c(1, 1), 1);
  EXPECT_EQ(c(1, 2), 0);
  EXPECT_EQ(c(2, 0), 0);
  EXPECT_EQ(c(2, 1), 0);
  EXPECT_EQ(c(2, 2), 0);

  a_cl = stan::math::to_matrix_cl(a);
  a_cl.view(stan::math::matrix_cl_view::Upper);
  b_cl.view(stan::math::matrix_cl_view::Lower);
  a_cl.sub_block(b_cl, 0, 1, 0, 1, 2, 2);
  EXPECT_EQ(a_cl.val().view(), stan::math::matrix_cl_view::Diagonal);
  c = stan::math::from_matrix_cl(a_cl.val());
  EXPECT_EQ(c(0, 0), 0);
  EXPECT_EQ(c(0, 1), 0);
  EXPECT_EQ(c(0, 2), 0);
  EXPECT_EQ(c(1, 0), 0);
  EXPECT_EQ(c(1, 1), 1);
  EXPECT_EQ(c(1, 2), 0);
  EXPECT_EQ(c(2, 0), 0);
  EXPECT_EQ(c(2, 1), 0);
  EXPECT_EQ(c(2, 2), 0);

  a_cl = stan::math::to_matrix_cl(a);
  a_cl.view(stan::math::matrix_cl_view::Upper);
  b_cl.view(stan::math::matrix_cl_view::Lower);
  a_cl.sub_block(b_cl, 0, 1, 0, 1, 3, 2);
  EXPECT_EQ(a_cl.val().view(), stan::math::matrix_cl_view::Lower);
  c = stan::math::from_matrix_cl(a_cl.val());
  EXPECT_EQ(c(0, 0), 0);
  EXPECT_EQ(c(0, 1), 0);
  EXPECT_EQ(c(0, 2), 0);
  EXPECT_EQ(c(1, 0), 0);
  EXPECT_EQ(c(1, 1), 1);
  EXPECT_EQ(c(1, 2), 0);
  EXPECT_EQ(c(2, 0), 0);
  EXPECT_EQ(c(2, 1), 1);
  EXPECT_EQ(c(2, 2), 1);

  a_cl = stan::math::to_matrix_cl(a);
  a_cl.view(stan::math::matrix_cl_view::Upper);
  b_cl.view(stan::math::matrix_cl_view::Lower);
  a_cl.sub_block(b_cl, 0, 0, 0, 0, 3, 3);
  EXPECT_EQ(a_cl.val().view(), stan::math::matrix_cl_view::Lower);
  c = stan::math::from_matrix_cl(a_cl.val());
  EXPECT_EQ(c(0, 0), 1);
  EXPECT_EQ(c(0, 1), 0);
  EXPECT_EQ(c(0, 2), 0);
  EXPECT_EQ(c(1, 0), 1);
  EXPECT_EQ(c(1, 1), 1);
  EXPECT_EQ(c(1, 2), 0);
  EXPECT_EQ(c(2, 0), 1);
  EXPECT_EQ(c(2, 1), 1);
  EXPECT_EQ(c(2, 2), 1);
}

#endif
