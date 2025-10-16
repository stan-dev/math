#include <stan/math.hpp>
#include <gtest/gtest.h>

TEST(revCore, stdIteratorTraits) {
  using stan::math::var;
  using traits = std::iterator_traits<stan::math::var>;
  var a;
  var b;
  var* a_ptr = &a;
  var* b_ptr = &b;
  traits::difference_type diff = a_ptr - b_ptr;
  EXPECT_EQ(a_ptr - b_ptr, diff);

  traits::value_type c = a;
  EXPECT_EQ(a.vi_, c.vi_);

  traits::pointer a_ptr_copy = a_ptr;
  EXPECT_EQ(a_ptr, a_ptr_copy);

  traits::reference d_ref = a;
  EXPECT_EQ(d_ref.vi_, a.vi_);
}
