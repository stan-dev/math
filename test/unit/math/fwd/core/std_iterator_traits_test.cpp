#include <stan/math/mix.hpp>
#include <gtest/gtest.h>

template <typename T>
void expect_iterator_traits() {
  using traits = std::iterator_traits<T>;
  T a;
  T b;
  T* a_ptr = &a;
  T* b_ptr = &b;
  typename traits::difference_type diff = a_ptr - b_ptr;
  EXPECT_EQ(a_ptr - b_ptr, diff);

  typename traits::value_type c = a;
  EXPECT_EQ(a.val(), c.val());
  EXPECT_EQ(a.tangent(), c.tangent());

  typename traits::pointer a_ptr_copy = a_ptr;
  EXPECT_EQ(a_ptr, a_ptr_copy);

  typename traits::reference d_ref = a;
  EXPECT_EQ(a.val(), d_ref.val());
  EXPECT_EQ(a.tangent(), d_ref.tangent());
}

TEST(fwdCore, stdIteratorTraits) {
  using stan::math::fvar;
  expect_iterator_traits<fvar<double>>();
  expect_iterator_traits<fvar<fvar<double>>>();
}
