#include <stan/math.hpp>
#include <stan/math/prim.hpp>
#include <boost/type_traits/has_nothrow_copy.hpp>
#include <test/unit/pretty_print_types.hpp>
#include <test/unit/math/rev/fun/util.hpp>
#include <test/unit/math/rev/core/gradable.hpp>
#include <gtest/gtest.h>
#include <string>
#include <vector>

struct AgradRev : public testing::Test {
  void SetUp() {
    // make sure memory's clean before starting each test
    stan::math::recover_memory();
  }
};

template <typename T, typename S>
void ctor_overloads_impl() {
  using stan::math::var_value;
  using stan::math::vari_value;
  using stan::math::test::type_name;
  // standard constructor
  EXPECT_FLOAT_EQ(3.7, var_value<T>(3.7).val())
      << "Failed For T: " << type_name<T>() << "\n";
  // make sure copy ctor is used rather than casting vari* to unsigned int
  EXPECT_FLOAT_EQ(12.3, var_value<T>(new vari_value<T>(12.3)).val())
      << "Failed For T: " << type_name<T>() << std::endl;
  // make sure rvalue var_value can be accepted
  EXPECT_FLOAT_EQ(12.3, var_value<T>(var_value<T>(12.3)).val())
      << "Failed For T: " << type_name<T>() << std::endl;
  // S type is preserved
  EXPECT_FLOAT_EQ(static_cast<S>(3.7), var_value<T>(static_cast<S>(3.7)).val())
      << "Failed For T: " << type_name<T>() << " and S: " << type_name<S>()
      << "\n";
  // Make sure integral types don't hold a nullptr instead of zero.
  EXPECT_FLOAT_EQ(0, var_value<T>(static_cast<S>(0)).val())
      << "Failed For T: " << type_name<T>() << " and S:" << type_name<S>()
      << "\n";
}

template <typename T>
void ctor_overloads() {
  ctor_overloads_impl<T, double>();
  ctor_overloads_impl<T, long double>();
  ctor_overloads_impl<T, float>();
  ctor_overloads_impl<T, bool>();
  ctor_overloads_impl<T, char>();
  ctor_overloads_impl<T, int>();
  ctor_overloads_impl<T, int16_t>();
  ctor_overloads_impl<T, int32_t>();
  ctor_overloads_impl<T, unsigned char>();
  ctor_overloads_impl<T, unsigned int>();
  ctor_overloads_impl<T, uint32_t>();
  ctor_overloads_impl<T, size_t>();
  ctor_overloads_impl<T, ptrdiff_t>();
}

TEST_F(AgradRev, ctorOverloads) {
  ctor_overloads<float>();
  ctor_overloads<double>();
  ctor_overloads<long double>();
}

TEST_F(AgradRev, a_eq_x) {
  AVAR a = 5.0;
  EXPECT_FLOAT_EQ(5.0, a.val());
}

TEST_F(AgradRev, a_of_x) {
  AVAR a(6.0);
  EXPECT_FLOAT_EQ(6.0, a.val());
}

TEST_F(AgradRev, a__a_eq_x) {
  AVAR a;
  a = 7.0;
  EXPECT_FLOAT_EQ(7.0, a.val());
}

TEST_F(AgradRev, eq_a) {
  AVAR a = 5.0;
  AVAR f = a;
  AVEC x = createAVEC(a);
  VEC dx;
  f.grad(x, dx);
  EXPECT_FLOAT_EQ(1.0, dx[0]);
}

TEST_F(AgradRev, a_ostream) {
  AVAR a = 6.0;
  std::ostringstream os;

  os << a;
  EXPECT_EQ("6", os.str());

  os.str("");
  a = 10.5;
  os << a;
  EXPECT_EQ("10.5", os.str());
}

TEST_F(AgradRev, smart_ptrs) {
  AVAR a = 2.0;
  EXPECT_FLOAT_EQ(2.0, (*a).val_);
  EXPECT_FLOAT_EQ(2.0, a->val_);

  EXPECT_FLOAT_EQ(2.0, (*a.vi_).val_);
  EXPECT_FLOAT_EQ(2.0, a.vi_->val_);
}

TEST_F(AgradRev, stackAllocation) {
  using stan::math::var;
  using stan::math::vari;

  vari ai(1.0);
  vari bi(2.0);

  var a(&ai);
  var b(&bi);

  AVEC x = createAVEC(a, b);
  var f = a * b;

  VEC g;
  f.grad(x, g);

  EXPECT_EQ(2U, g.size());
  EXPECT_FLOAT_EQ(2.0, g[0]);
  EXPECT_FLOAT_EQ(1.0, g[1]);
}

TEST_F(AgradRev, print) {
  using stan::math::var;

  std::ostringstream output;
  std::string str;

  var initialized_var(0);
  output << initialized_var;
  str = output.str();
  EXPECT_STREQ("0", output.str().c_str());

  output.clear();
  output.str("");
  var uninitialized_var;
  output << uninitialized_var;
  str = output.str();
  EXPECT_STREQ("uninitialized", output.str().c_str());
}

// should really be doing this test with a mock object using ctor
// vari_(double, bool);  as in:
//
// struct nostack_test_vari : public stan::math::vari {
//   nostack_test_vari(double x)
//   : stan::math::vari(x, false) {
//   }
//   void chain() {
//     // no op on the chain
//   }
// };

// struct both_test_vari : public stan::math::vari {
//   both_test_vari(vari* vi, vari* bi) {

//   }
// };

// var foo(var y, var z) {
//   return y *
// }

TEST_F(AgradRev, basicGradient1) {
  using stan::math::recover_memory;

  for (int i = 0; i < 100; ++i) {
    gradable g = setup_simple();
    g.test();
    recover_memory();
  }
}

TEST_F(AgradRev, basicGradient2) {
  using stan::math::recover_memory;

  for (int i = 0; i < 100; ++i) {
    gradable g = setup_quad_form();
    g.test();
    recover_memory();
  }
}

TEST_F(AgradRev, nestedGradient1) {
  using stan::math::recover_memory;
  using stan::math::recover_memory_nested;
  using stan::math::start_nested;

  gradable g0 = setup_simple();

  start_nested();
  gradable g1 = setup_quad_form();
  g1.test();
  recover_memory_nested();

  start_nested();
  gradable g2 = setup_simple();
  g2.test();
  recover_memory_nested();

  g0.test();
  recover_memory();
}

TEST_F(AgradRev, nestedGradient2) {
  using stan::math::recover_memory;
  using stan::math::recover_memory_nested;
  using stan::math::start_nested;

  gradable g0 = setup_quad_form();

  start_nested();
  gradable g1 = setup_simple();
  g1.test();
  recover_memory_nested();

  start_nested();
  gradable g2 = setup_quad_form();
  g2.test();
  recover_memory_nested();

  g0.test();
  recover_memory();
}

TEST_F(AgradRev, nestedGradient3) {
  using stan::math::recover_memory;
  using stan::math::recover_memory_nested;
  using stan::math::start_nested;

  start_nested();
  gradable g1 = setup_simple();
  start_nested();
  gradable g2 = setup_quad_form();
  start_nested();
  gradable g3 = setup_quad_form();
  start_nested();
  gradable g4 = setup_simple();
  g4.test();
  recover_memory_nested();
  g3.test();
  recover_memory_nested();
  g2.test();
  recover_memory_nested();
  g1.test();
  recover_memory_nested();
  recover_memory();
}

TEST_F(AgradRev, grad) {
  AVAR a = 5.0;
  AVAR b = 10.0;
  AVAR f = a * b + a;

  EXPECT_NO_THROW(f.grad()) << "testing the grad function with no args";

  EXPECT_FLOAT_EQ(5.0, a.val());
  EXPECT_FLOAT_EQ(10.0, b.val());
  EXPECT_FLOAT_EQ(55.0, f.val());

  EXPECT_FLOAT_EQ(1.0, f.adj());
  EXPECT_FLOAT_EQ(11.0, a.adj());
  EXPECT_FLOAT_EQ(5.0, b.adj());
}
