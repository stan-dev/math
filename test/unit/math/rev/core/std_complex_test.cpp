#include <stan/math/rev/scal.hpp>
#include <stan/math/rev/core/std_complex.hpp>
#include <gtest/gtest.h>

// For a definition of the spec:
// https://en.cppreference.com/w/cpp/numeric/complex
class MathRev : public testing::Test {
 public:
  void SetUp() { stan::math::recover_memory(); }
};

TEST_F(MathRev, complex_constructor) {
  std::complex<stan::math::var> x;
  EXPECT_EQ(2, stan::math::ChainableStack::instance().var_stack_.size());
  EXPECT_FLOAT_EQ(0, x.real().val());
  EXPECT_FLOAT_EQ(0, x.imag().val());
  stan::math::recover_memory();

  std::complex<stan::math::var> y{0};
  EXPECT_EQ(2, stan::math::ChainableStack::instance().var_stack_.size());
  stan::math::recover_memory();

  std::complex<stan::math::var> z{0, 0};
  EXPECT_EQ(2, stan::math::ChainableStack::instance().var_stack_.size());
  stan::math::recover_memory();
}

TEST_F(MathRev, assignment_double) {
  double rhs = 1;

  std::complex<stan::math::var> lhs;
  lhs = rhs;
  EXPECT_FLOAT_EQ(1, lhs.real().val());
  EXPECT_FLOAT_EQ(0, lhs.imag().val());
  EXPECT_EQ(4, stan::math::ChainableStack::instance().var_stack_.size());
}

TEST_F(MathRev, assignment_var) {
  stan::math::var rhs = 2;

  std::complex<stan::math::var> lhs;
  lhs = rhs;
  EXPECT_FLOAT_EQ(2, lhs.real().val());
  EXPECT_FLOAT_EQ(0, lhs.imag().val());
  EXPECT_EQ(4, stan::math::ChainableStack::instance().var_stack_.size());
  EXPECT_EQ(rhs.vi_, lhs.real().vi_);
}

TEST_F(MathRev, assignment_complex_real_only) {
  std::complex<stan::math::var> rhs{3};

  std::complex<stan::math::var> lhs;
  lhs = rhs;
  EXPECT_FLOAT_EQ(3, lhs.real().val());
  EXPECT_FLOAT_EQ(0, lhs.imag().val());
  EXPECT_EQ(4, stan::math::ChainableStack::instance().var_stack_.size());
  EXPECT_EQ(rhs.real().vi_, lhs.real().vi_);
  EXPECT_EQ(rhs.imag().vi_, lhs.imag().vi_);
}

TEST_F(MathRev, assignment_complex) {
  std::complex<stan::math::var> rhs{4, 5};

  std::complex<stan::math::var> lhs;
  lhs = rhs;
  EXPECT_FLOAT_EQ(4, lhs.real().val());
  EXPECT_FLOAT_EQ(5, lhs.imag().val());
  EXPECT_EQ(4, stan::math::ChainableStack::instance().var_stack_.size());
  EXPECT_EQ(rhs.real().vi_, lhs.real().vi_);
  EXPECT_EQ(rhs.imag().vi_, lhs.imag().vi_);
}

TEST_F(MathRev, real_component) {
  std::complex<stan::math::var> c{1, 2};

  EXPECT_EQ(1, c.real().val());
}

TEST_F(MathRev, imag_component) {
  std::complex<stan::math::var> c{1, 2};

  EXPECT_EQ(2, c.imag().val());
}

TEST_F(MathRev, member_operators) {
  std::complex<stan::math::var> x{1, 2}, y{3, 4};

  // operator+=
  EXPECT_EQ(4, stan::math::ChainableStack::instance().var_stack_.size());
  x += y;
  EXPECT_FLOAT_EQ(4, x.real().val());
  EXPECT_FLOAT_EQ(6, x.imag().val());
  EXPECT_FLOAT_EQ(3, y.real().val());
  EXPECT_FLOAT_EQ(4, y.imag().val());
  EXPECT_EQ(6, stan::math::ChainableStack::instance().var_stack_.size());
  stan::math::recover_memory();

  // operator-=
  x = std::complex<stan::math::var>{1, 2};
  y = std::complex<stan::math::var>{3, 4};
  EXPECT_EQ(4, stan::math::ChainableStack::instance().var_stack_.size());
  x -= y;
  EXPECT_FLOAT_EQ(-2, x.real().val());
  EXPECT_FLOAT_EQ(-2, x.imag().val());
  EXPECT_FLOAT_EQ(3, y.real().val());
  EXPECT_FLOAT_EQ(4, y.imag().val());
  EXPECT_EQ(6, stan::math::ChainableStack::instance().var_stack_.size());
  stan::math::recover_memory();

  // operator/=
  x = std::complex<stan::math::var>{1, 2};
  y = std::complex<stan::math::var>{3, 4};
  EXPECT_EQ(4, stan::math::ChainableStack::instance().var_stack_.size());
  x /= y;
  EXPECT_FLOAT_EQ(1.0 / 3.0, x.real().val());
  EXPECT_FLOAT_EQ(0.5, x.imag().val());
  EXPECT_FLOAT_EQ(3, y.real().val());
  EXPECT_FLOAT_EQ(4, y.imag().val());
  EXPECT_EQ(6, stan::math::ChainableStack::instance().var_stack_.size());
  stan::math::recover_memory();

  // operator*=
  x = std::complex<stan::math::var>{1, 2};
  y = std::complex<stan::math::var>{3, 4};
  EXPECT_EQ(4, stan::math::ChainableStack::instance().var_stack_.size());
  x *= y;
  EXPECT_FLOAT_EQ(3, x.real().val());
  EXPECT_FLOAT_EQ(8, x.imag().val());
  EXPECT_FLOAT_EQ(3, y.real().val());
  EXPECT_FLOAT_EQ(4, y.imag().val());
  EXPECT_EQ(6, stan::math::ChainableStack::instance().var_stack_.size());
  stan::math::recover_memory();
}

TEST_F(MathRev, unary_operators) {
  std::complex<stan::math::var> x{1, 2};
  EXPECT_EQ(2, stan::math::ChainableStack::instance().var_stack_.size());

  // operator+
  std::complex<stan::math::var> z = +x;
  EXPECT_EQ(2, stan::math::ChainableStack::instance().var_stack_.size());
  EXPECT_FLOAT_EQ(1, z.real().val());
  EXPECT_FLOAT_EQ(2, z.imag().val());
  EXPECT_FLOAT_EQ(1, x.real().val());
  EXPECT_FLOAT_EQ(2, x.imag().val());
  stan::math::recover_memory();

  // operator-
  x = std::complex<stan::math::var>{1, 2};
  z = -x;

  EXPECT_EQ(4, stan::math::ChainableStack::instance().var_stack_.size());
  EXPECT_FLOAT_EQ(-1, z.real().val());
  EXPECT_FLOAT_EQ(-2, z.imag().val());
  EXPECT_FLOAT_EQ(1, x.real().val());
  EXPECT_FLOAT_EQ(2, x.imag().val());
  stan::math::recover_memory();
}

TEST_F(MathRev, arithmetic) {
  std::complex<stan::math::var> x{1, 2}, y{3, 4};
  EXPECT_EQ(4, stan::math::ChainableStack::instance().var_stack_.size());

  // operator+
  std::complex<stan::math::var> z = x - y;
  EXPECT_EQ(6, stan::math::ChainableStack::instance().var_stack_.size());
  EXPECT_FLOAT_EQ(-2, z.real().val());
  EXPECT_FLOAT_EQ(-2, z.imag().val());
  EXPECT_FLOAT_EQ(1, x.real().val());
  EXPECT_FLOAT_EQ(2, x.imag().val());
  EXPECT_FLOAT_EQ(3, y.real().val());
  EXPECT_FLOAT_EQ(4, y.imag().val());
  stan::math::recover_memory();

  // operator-
  x = std::complex<stan::math::var>{1, 2};
  y = std::complex<stan::math::var>{3, 4};
  z = x + y;

  EXPECT_EQ(6, stan::math::ChainableStack::instance().var_stack_.size());
  EXPECT_FLOAT_EQ(4, z.real().val());
  EXPECT_FLOAT_EQ(6, z.imag().val());
  EXPECT_FLOAT_EQ(1, x.real().val());
  EXPECT_FLOAT_EQ(2, x.imag().val());
  EXPECT_FLOAT_EQ(3, y.real().val());
  EXPECT_FLOAT_EQ(4, y.imag().val());
  stan::math::recover_memory();

  // operator*
  x = std::complex<stan::math::var>{1, 2};
  y = std::complex<stan::math::var>{3, 4};
  z = x * y;

  EXPECT_EQ(6, stan::math::ChainableStack::instance().var_stack_.size());
  EXPECT_FLOAT_EQ(3, z.real().val());
  EXPECT_FLOAT_EQ(8, z.imag().val());
  EXPECT_FLOAT_EQ(1, x.real().val());
  EXPECT_FLOAT_EQ(2, x.imag().val());
  EXPECT_FLOAT_EQ(3, y.real().val());
  EXPECT_FLOAT_EQ(4, y.imag().val());
  stan::math::recover_memory();

  // operator/
  x = std::complex<stan::math::var>{1, 2};
  y = std::complex<stan::math::var>{3, 4};
  z = x / y;

  EXPECT_EQ(6, stan::math::ChainableStack::instance().var_stack_.size());
  EXPECT_FLOAT_EQ(1.0 / 3.0, z.real().val());
  EXPECT_FLOAT_EQ(0.5, z.imag().val());
  EXPECT_FLOAT_EQ(1, x.real().val());
  EXPECT_FLOAT_EQ(2, x.imag().val());
  EXPECT_FLOAT_EQ(3, y.real().val());
  EXPECT_FLOAT_EQ(4, y.imag().val());
  stan::math::recover_memory();
}

TEST_F(MathRev, comparison) {
  // complex and scalar; assume var and double as scalar?
  // operator==
  std::complex<stan::math::var> x1{1, 2}, y1{3, 4};
  std::complex<stan::math::var> x2{1, 0}, y2{3, 0};
  double x_d = 1, y_d = 3;
  stan::math::var x_v = 1, y_v = 3;

  const int stack_size
      = stan::math::ChainableStack::instance().var_stack_.size();
  EXPECT_TRUE(x1 == x1);
  EXPECT_EQ(stack_size,
            stan::math::ChainableStack::instance().var_stack_.size());
  EXPECT_FALSE(x1 == x2);
  EXPECT_EQ(stack_size,
            stan::math::ChainableStack::instance().var_stack_.size());
  EXPECT_FALSE(x1 == y1);
  EXPECT_EQ(stack_size,
            stan::math::ChainableStack::instance().var_stack_.size());

  EXPECT_FALSE(x1 == x_d);
  EXPECT_EQ(stack_size,
            stan::math::ChainableStack::instance().var_stack_.size());
  EXPECT_FALSE(x_d == x1);
  EXPECT_EQ(stack_size,
            stan::math::ChainableStack::instance().var_stack_.size());
  EXPECT_TRUE(x2 == x_d);
  EXPECT_EQ(stack_size,
            stan::math::ChainableStack::instance().var_stack_.size());
  EXPECT_TRUE(x_d == x2);
  EXPECT_EQ(stack_size,
            stan::math::ChainableStack::instance().var_stack_.size());
  EXPECT_FALSE(x1 == y_d);
  EXPECT_EQ(stack_size,
            stan::math::ChainableStack::instance().var_stack_.size());

  EXPECT_FALSE(x1 == x_v);
  EXPECT_EQ(stack_size,
            stan::math::ChainableStack::instance().var_stack_.size());
  EXPECT_FALSE(x_v == x1);
  EXPECT_EQ(stack_size,
            stan::math::ChainableStack::instance().var_stack_.size());
  EXPECT_TRUE(x2 == x_v);
  EXPECT_EQ(stack_size,
            stan::math::ChainableStack::instance().var_stack_.size());
  EXPECT_TRUE(x_v == x2);
  EXPECT_EQ(stack_size,
            stan::math::ChainableStack::instance().var_stack_.size());
  EXPECT_FALSE(x1 == y_v);
  EXPECT_EQ(stack_size,
            stan::math::ChainableStack::instance().var_stack_.size());

  // operator!=
  EXPECT_FALSE(x1 != x1);
  EXPECT_EQ(stack_size,
            stan::math::ChainableStack::instance().var_stack_.size());
  EXPECT_TRUE(x1 != x2);
  EXPECT_EQ(stack_size,
            stan::math::ChainableStack::instance().var_stack_.size());
  EXPECT_TRUE(x1 != y1);
  EXPECT_EQ(stack_size,
            stan::math::ChainableStack::instance().var_stack_.size());

  EXPECT_TRUE(x1 != x_d);
  EXPECT_EQ(stack_size,
            stan::math::ChainableStack::instance().var_stack_.size());
  EXPECT_TRUE(x_d != x1);
  EXPECT_EQ(stack_size,
            stan::math::ChainableStack::instance().var_stack_.size());
  EXPECT_FALSE(x2 != x_d);
  EXPECT_EQ(stack_size,
            stan::math::ChainableStack::instance().var_stack_.size());
  EXPECT_FALSE(x_d != x2);
  EXPECT_EQ(stack_size,
            stan::math::ChainableStack::instance().var_stack_.size());
  EXPECT_TRUE(x1 != y_d);
  EXPECT_EQ(stack_size,
            stan::math::ChainableStack::instance().var_stack_.size());

  EXPECT_TRUE(x1 != x_v);
  EXPECT_EQ(stack_size,
            stan::math::ChainableStack::instance().var_stack_.size());
  EXPECT_TRUE(x_v != x1);
  EXPECT_EQ(stack_size,
            stan::math::ChainableStack::instance().var_stack_.size());
  EXPECT_FALSE(x2 != x_v);
  EXPECT_EQ(stack_size,
            stan::math::ChainableStack::instance().var_stack_.size());
  EXPECT_FALSE(x_v != x2);
  EXPECT_EQ(stack_size,
            stan::math::ChainableStack::instance().var_stack_.size());
  EXPECT_TRUE(x1 != y_v);
  EXPECT_EQ(stack_size,
            stan::math::ChainableStack::instance().var_stack_.size());
}

TEST_F(MathRev, serialize_deserialize) {
  std::stringstream msg;

  std::complex<stan::math::var> x{1, 2};
  int stack_size = stan::math::ChainableStack::instance().var_stack_.size();

  // operator<<
  msg << x;
  EXPECT_EQ("(1,2)", msg.str());
  EXPECT_EQ(stack_size,
            stan::math::ChainableStack::instance().var_stack_.size());

  // operator>>
  std::complex<stan::math::var> y;
  stack_size = stan::math::ChainableStack::instance().var_stack_.size();
  msg >> y;
  EXPECT_EQ(1, y.real().val());
  EXPECT_EQ(2, y.imag().val());
  EXPECT_EQ(stack_size + 2,
            stan::math::ChainableStack::instance().var_stack_.size());
}

TEST_F(MathRev, real) {
  std::complex<stan::math::var> x{1, 2};

  EXPECT_EQ(2, stan::math::ChainableStack::instance().var_stack_.size());

  EXPECT_EQ(1, x.real().val());
  EXPECT_EQ(2, stan::math::ChainableStack::instance().var_stack_.size());

  EXPECT_NO_THROW(x.real(2));
  EXPECT_EQ(2, x.real().val());
  EXPECT_EQ(3, stan::math::ChainableStack::instance().var_stack_.size());

  stan::math::var y = 3;
  EXPECT_NO_THROW(x.real(y));
  EXPECT_EQ(y, x.real().val());
  EXPECT_EQ(4, stan::math::ChainableStack::instance().var_stack_.size());
}

TEST_F(MathRev, imag) {
  std::complex<stan::math::var> x{1, 2};

  EXPECT_EQ(2, stan::math::ChainableStack::instance().var_stack_.size());

  EXPECT_EQ(2, x.imag().val());
  EXPECT_EQ(2, stan::math::ChainableStack::instance().var_stack_.size());

  EXPECT_NO_THROW(x.imag(-1));
  EXPECT_EQ(-1, x.imag().val());
  EXPECT_EQ(3, stan::math::ChainableStack::instance().var_stack_.size());

  stan::math::var y = 3;
  EXPECT_NO_THROW(x.imag(y));
  EXPECT_EQ(y, x.imag().val());
  EXPECT_EQ(4, stan::math::ChainableStack::instance().var_stack_.size());
}

TEST_F(MathRev, abs) {
  std::complex<stan::math::var> z{3, 4};
  stan::math::var f = abs(z);
  EXPECT_EQ(5, f.val());
  std::vector<stan::math::var> x{real(z)};
  std::vector<double> g;
  f.grad(x, g);
  EXPECT_FLOAT_EQ(0.6, g[0]);
}

TEST_F(MathRev, arg) {
  std::complex<stan::math::var> z{1, stan::math::pi()};
  stan::math::var f = arg(z);
  EXPECT_EQ(f.val(), stan::math::atan2(std::imag(z), std::real(z)));
}

TEST_F(MathRev, isinf) {
  std::complex<stan::math::var> a{0, 0};
  std::complex<double> b{0, 0};
  EXPECT_EQ(std::isinf(b), std::isinf(a));
  EXPECT_FALSE(std::isinf(a));

  a.real(std::numeric_limits<double>::infinity());
  a.imag(0);
  b.real(std::numeric_limits<double>::infinity());
  b.imag(0);
  EXPECT_EQ(std::isinf(b), std::isinf(a));
  EXPECT_TRUE(std::isinf(a));

  a.real(0);
  a.imag(std::numeric_limits<double>::infinity());
  b.real(0);
  b.imag(std::numeric_limits<double>::infinity());
  EXPECT_EQ(std::isinf(b), std::isinf(a));
  EXPECT_FALSE(std::isinf(a));
}

TEST_F(MathRev, isnan) {
  std::complex<stan::math::var> a{0, 0};
  std::complex<double> b{0, 0};
  EXPECT_EQ(std::isnan(b), std::isnan(a));
  EXPECT_FALSE(std::isnan(a));

  a.real(std::numeric_limits<double>::quiet_NaN());
  a.imag(0);
  b.real(std::numeric_limits<double>::quiet_NaN());
  b.imag(0);
  EXPECT_EQ(std::isnan(b), std::isnan(a));
  EXPECT_TRUE(std::isnan(a));

  a.real(0);
  a.imag(std::numeric_limits<double>::quiet_NaN());
  b.real(0);
  b.imag(std::numeric_limits<double>::quiet_NaN());
  EXPECT_EQ(std::isnan(b), std::isnan(a));
  EXPECT_FALSE(std::isnan(a));
}

TEST_F(MathRev, norm) {
  std::complex<stan::math::var> z{3, 4};
  stan::math::var f = norm(z);
  EXPECT_EQ(f.val(), 25);

  std::vector<stan::math::var> x{imag(z)};
  std::vector<double> g;
  f.grad(x, g);
  EXPECT_FLOAT_EQ(g[0], 8);
}

TEST_F(MathRev, conj) {
  std::complex<stan::math::var> z(1, 2);
  std::complex<stan::math::var> f = conj(z);
  EXPECT_EQ(real(f).val(), 1);
  EXPECT_EQ(imag(f).val(), -2);
}

TEST_F(MathRev, proj) {
  std::complex<stan::math::var> z1(1, 2);
  std::complex<stan::math::var> f = proj(z1);
  EXPECT_TRUE(f == z1) << f << std::endl;

  using stan::math::positive_infinity;
  std::complex<stan::math::var> z2(positive_infinity(), -1);
  f = proj(z2);
  EXPECT_TRUE(f == std::complex<stan::math::var>(positive_infinity(), -0))
      << f << std::endl;

  using stan::math::negative_infinity;
  std::complex<stan::math::var> z3(0, negative_infinity());
  f = proj(z3);
  EXPECT_TRUE(f == std::complex<stan::math::var>(positive_infinity(), -0))
      << f << std::endl;

  std::complex<stan::math::var> z4(std::numeric_limits<double>::quiet_NaN(),
                                   -2.0);
  f = proj(z4);
  EXPECT_TRUE(stan::math::is_nan(f.real()));
  EXPECT_FLOAT_EQ(-2.0, f.imag().val());
}

TEST_F(MathRev, polar) {
  std::complex<stan::math::var> z = std::polar(1, 0);
  stan::math::var f = arg(z);
  EXPECT_EQ(f.val(), 0.0);
}
