#include <stan/math/rev/core.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <stan/math/rev/scal/fun/value_of_rec.hpp>
#include <stan/math/prim/arr/fun/value_of_rec.hpp>
#include <gtest/gtest.h>

template <typename T>
void fill(const std::vector<double>& contents,
          std::vector<T>& x){
  x.assign(contents.size(), T());
  for (size_t i = 0; i < contents.size(); ++i)
    x[i] = T(contents[i]);
}

TEST(MathMatrix,value_of_rec) {
  using stan::math::value_of_rec;
  using std::vector;
  using stan::math::var;

  vector<double> a_vals;

  for (size_t i = 0; i < 10; ++i)
    a_vals.push_back(i + 1);

  vector<double> b_vals;

  for (size_t i = 10; i < 15; ++i)
    b_vals.push_back(i + 1);

  vector<var> a;
  fill(a_vals, a);
  vector<var> b;
  fill(b_vals, b);

  vector<double> d_a = value_of_rec(a);
  vector<double> d_b = value_of_rec(b);

  for (int i = 0; i < 5; ++i)
    EXPECT_FLOAT_EQ(b[i].val(), d_b[i]);

  for (int i = 0; i < 10; ++i)
    EXPECT_FLOAT_EQ(a[i].val(), d_a[i]);
}
