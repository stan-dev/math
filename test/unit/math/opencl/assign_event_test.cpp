#ifdef STAN_OPENCL

#include <stan/math/opencl/prim.hpp>
#include <stan/math/opencl/kernel_cl.hpp>
#include <gtest/gtest.h>

using stan::math::matrix_cl;
using stan::math::opencl_kernels::in_buffer;
using stan::math::opencl_kernels::in_out_buffer;
using stan::math::opencl_kernels::out_buffer;
using stan::math::opencl_kernels::internal::assign_events;

TEST(assign_event, correct_vectors) {
  matrix_cl<double> m;
  // pointers not set up to work yet; TBD if needed
  // matrix_cl *mp = &m;
  cl::Event e;
  assign_events<in_buffer>(e, m);
  EXPECT_EQ(m.read_events().size(), 1);
  EXPECT_EQ(m.write_events().size(), 0);

  assign_events<out_buffer>(e, m);
  EXPECT_EQ(m.read_events().size(), 1);
  EXPECT_EQ(m.write_events().size(), 1);

  assign_events<in_out_buffer>(e, m);
  EXPECT_EQ(m.read_events().size(), 2);
  EXPECT_EQ(m.write_events().size(), 2);

  m.clear_read_write_events();
}
#endif
