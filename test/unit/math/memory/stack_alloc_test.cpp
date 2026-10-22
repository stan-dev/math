#include <gtest/gtest.h>
#include <stan/math/memory/stack_alloc.hpp>
#include <stdlib.h>
#include <cmath>
#include <utility>
#include <vector>

TEST(MemoryStackAlloc, allocArray) {
  // just an example to show how alloc_array is used
  stan::math::stack_alloc allocator;
  double* x = allocator.alloc_array<double>(10U);
  for (int i = 0; i < 10; ++i)
    x[i] = 3.0;
  for (int i = 0; i < 10; ++i)
    EXPECT_EQ(3.0, x[i]);
}

struct biggy {
  double r[10];
};

TEST(MemoryStackAlloc, allocArrayBigger) {
  size_t N = 1000;
  size_t K = 10;
  stan::math::stack_alloc allocator;
  biggy* x = allocator.alloc_array<biggy>(N);
  for (size_t i = 0; i < N; ++i)
    for (size_t k = 0; k < K; ++k)
      x[i].r[k] = k * i;
  for (size_t i = 0; i < N; ++i)
    for (size_t k = 0; k < K; ++k)
      EXPECT_FLOAT_EQ(k * i, x[i].r[k]);
}
TEST(stack_alloc, bytes_allocated) {
  stan::math::stack_alloc allocator;
  EXPECT_LE(0L, allocator.bytes_allocated());
  for (size_t n = 1; n <= 10000; ++n) {
    allocator.alloc(n);
    size_t bytes_requested = (n * (n + 1)) / 2;
    size_t bytes_allocated = allocator.bytes_allocated();
    EXPECT_TRUE(bytes_requested <= bytes_allocated)
        << "bytes_requested: " << bytes_requested << std::endl
        << "bytes_allocated: " << bytes_allocated;
    // 1 << 16 is initial allocation;  *3 is to account for slop at end
    EXPECT_TRUE(bytes_allocated < ((1 << 16) + bytes_requested * 3));
  }
}

TEST(stack_alloc, is_aligned) {
  char* ptr = static_cast<char*>(malloc(1024));
  EXPECT_TRUE(stan::math::is_aligned(ptr, 1U));
  EXPECT_TRUE(stan::math::is_aligned(ptr, 2U));
  EXPECT_TRUE(stan::math::is_aligned(ptr, 4U));
  EXPECT_TRUE(stan::math::is_aligned(ptr, 8U));

  EXPECT_FALSE(stan::math::is_aligned(ptr + 1, 8U));
  // not very safe, but just a test
  free(ptr);
}

TEST(stack_alloc, alloc) {
  std::vector<double*> ds;
  std::vector<int*> is;
  std::vector<char*> cs;

  stan::math::stack_alloc allocator;

  for (int i = 0; i < 100000; ++i) {
    allocator.alloc(1317);
    double* foo = reinterpret_cast<double*>(allocator.alloc(sizeof(double)));
    *foo = 9.0;
    ds.push_back(foo);
    int* bar = reinterpret_cast<int*>(allocator.alloc(sizeof(int)));
    *bar = 17;
    is.push_back(bar);
    char* baz = reinterpret_cast<char*>(allocator.alloc(sizeof(char)));
    *baz = 3;
    cs.push_back(baz);
    allocator.alloc(13);

    EXPECT_FLOAT_EQ(9.0, *foo);
    EXPECT_EQ(17, *bar);
    EXPECT_EQ(3, *baz);
  }
  for (int i = 0; i < 10000; ++i) {
    EXPECT_FLOAT_EQ(9.0, *ds[i]);
    EXPECT_EQ(17, *is[i]);
    EXPECT_EQ(3, *cs[i]);
  }

  allocator.recover_all();
}

TEST(stack_alloc, alloc_aligned) {
  stan::math::stack_alloc allocator;
  int* x = allocator.alloc_array<int>(3);

  double* y = allocator.alloc_array<double>(4);
  EXPECT_TRUE(stan::math::is_aligned(y, 8));
  allocator.recover_all();
}

TEST(stack_alloc, alloc_over_aligned) {
  stan::math::stack_alloc allocator;
  // leave the next free byte 8 but not 64 byte aligned
  allocator.alloc(8);
  for (size_t alignment : {16U, 32U, 64U}) {
    double* x
        = static_cast<double*>(allocator.alloc(4 * sizeof(double), alignment));
    EXPECT_TRUE(stan::math::is_aligned(x, alignment)) << alignment;
    for (int i = 0; i < 4; ++i) {
      x[i] = i;
    }
    EXPECT_FLOAT_EQ(3.0, x[3]);
    // plain allocations after an aligned one stay 8 byte aligned
    EXPECT_TRUE(stan::math::is_aligned(allocator.alloc(sizeof(double)), 8U));
  }
}

TEST(stack_alloc, alloc_over_aligned_new_blocks) {
  stan::math::stack_alloc allocator;
  // larger than the initial block, so every allocation needs a new block
  const size_t big = 1 << 17;
  for (int i = 0; i < 4; ++i) {
    allocator.alloc(8);
    char* x = static_cast<char*>(allocator.alloc(big, 64));
    EXPECT_TRUE(stan::math::is_aligned(x, 64U));
    x[0] = 1;
    x[big - 1] = 2;
    EXPECT_TRUE(allocator.in_stack(x));
  }
  // reused blocks after recovering memory are aligned too
  allocator.recover_all();
  for (int i = 0; i < 4; ++i) {
    allocator.alloc(8);
    EXPECT_TRUE(stan::math::is_aligned(allocator.alloc(big, 64), 64U));
  }
}

TEST(stack_alloc, in_stack) {
  stan::math::stack_alloc allocator;

  double* x = allocator.alloc_array<double>(10U);
  double y = 0;

  EXPECT_TRUE(allocator.in_stack(x));
  EXPECT_FALSE(allocator.in_stack(&y));
}

TEST(stack_alloc, in_stack_second_block) {
  stan::math::stack_alloc allocator;

  char* x = allocator.alloc_array<char>(
      stan::math::internal::DEFAULT_INITIAL_NBYTES);
  EXPECT_TRUE(allocator.in_stack(x));
  EXPECT_FALSE(
      allocator.in_stack(x + stan::math::internal::DEFAULT_INITIAL_NBYTES));

  char* y = allocator.alloc_array<char>(1);
  EXPECT_TRUE(allocator.in_stack(x));
  EXPECT_TRUE(allocator.in_stack(y));
  // effect of padding
  EXPECT_TRUE(allocator.in_stack(y + 1));
  EXPECT_TRUE(allocator.in_stack(y + 7));

  EXPECT_FALSE(allocator.in_stack(y + 8));

  allocator.recover_all();
  EXPECT_FALSE(allocator.in_stack(x));
  EXPECT_FALSE(allocator.in_stack(y));
}

TEST(stack_alloc, aligned_alloc_and_free) {
  for (std::size_t alignment : {8, 16, 32, 64, 128}) {
    void* ptr = STAN_ALIGNED_ALLOC(alignment, 4 * alignment);
    ASSERT_NE(ptr, nullptr) << alignment;
    EXPECT_TRUE(stan::math::is_aligned(ptr, alignment)) << alignment;
    STAN_ALIGNED_FREE(ptr);
  }
}

TEST(stack_alloc, first_block_aligned) {
  stan::math::stack_alloc allocator;
  EXPECT_TRUE(stan::math::is_aligned(allocator.alloc(1), 8U));
}
