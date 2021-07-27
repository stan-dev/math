#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(AgradRev, arena_matrix_cl_shallow_copies) {
  stan::math::arena_matrix_cl<double> a(3, 2);
  stan::math::arena_matrix_cl<double> b(a);
  stan::math::arena_matrix_cl<double> c;
  c = a;
  EXPECT_EQ(a.buffer()(), b.buffer()());
  EXPECT_EQ(a.buffer()(), c.buffer()());
}

TEST(AgradRev, arena_matrix_cl_to_matrix_cl_conversion) {
  stan::math::arena_matrix_cl<double> a(3, 2);
  const stan::math::matrix_cl<double>& b(a);
  EXPECT_EQ(a.buffer()(), b.buffer()());
}

TEST(AgradRev, arena_matrix_cl_to_matrix_cl_move_construction) {
  stan::math::arena_matrix_cl<double> a(3, 2);
  cl::Buffer a_buf = a.buffer();
  stan::math::matrix_cl<double> b(std::move(a));
  EXPECT_EQ(a_buf(), a.buffer()());
  EXPECT_EQ(a_buf(), b.buffer()());
}

TEST(AgradRev, arena_matrix_cl_to_matrix_cl_move_assignment) {
  stan::math::arena_matrix_cl<double> a(3, 2);
  cl::Buffer a_buf = a.buffer();
  stan::math::matrix_cl<double> b;
  b = std::move(a);
  EXPECT_EQ(a_buf(), a.buffer()());
  EXPECT_EQ(a_buf(), b.buffer()());
}

#endif
