#include <stan/math/rev/core.hpp>
#include <test/unit/math/rev/core/gradable.hpp>
#include <gtest/gtest.h>
#include <vector>
#include <memory>

struct AgradLocalScoped : public testing::Test {
  void SetUp() {
    // make sure memory's clean before starting each test
    stan::math::recover_memory();
  }
};

TEST_F(AgradLocalScoped, scoped_chainablestack_base) {
  using stan::math::nested_rev_autodiff;
  using stan::math::ScopedChainableStack;
  using stan::math::var;

  ScopedChainableStack scoped_stack;

  {
    nested_rev_autodiff nested;

    double cgrad_a = scoped_stack.execute([] {
      stan::math::start_nested();
      var a = 2.0;
      var b = 4.0;
      var c = a * a + b;
      c.grad();
      EXPECT_GT(stan::math::nested_size(), 0);
      return a.adj();
    });
    EXPECT_FLOAT_EQ(cgrad_a, 4.0);

    EXPECT_EQ(stan::math::nested_size(), 0);
  }

  // the nested autodiff stack went out of scope, but that did not
  // touch the scoped stack

  scoped_stack.execute([] { EXPECT_GT(stan::math::nested_size(), 0); });
}

TEST_F(AgradLocalScoped, scoped_chainablestack_functor) {
  using stan::math::nested_rev_autodiff;
  using stan::math::ScopedChainableStack;
  using stan::math::var;

  ScopedChainableStack scoped_stack;

  struct scoped_functor {
    double a_val_;
    double cgrad_a_;

    explicit scoped_functor(double a_val) : a_val_(a_val), cgrad_a_(0.0) {}

    void operator()() {
      stan::math::start_nested();
      var a = a_val_;
      var b = 4.0;
      var c = a * a + b;
      c.grad();
      EXPECT_GT(stan::math::nested_size(), 0);
      cgrad_a_ = a.adj();
    }
  } worker(2.0);

  {
    nested_rev_autodiff nested;

    scoped_stack.execute(worker);
    EXPECT_FLOAT_EQ(worker.cgrad_a_, 4.0);

    EXPECT_EQ(stan::math::nested_size(), 0);
  }

  // the nested autodiff stack went out of scope, but that did not
  // touch the scoped stack

  scoped_stack.execute([] { EXPECT_GT(stan::math::nested_size(), 0); });
}

TEST_F(AgradLocalScoped, scoped_chainablestack_simple) {
  using stan::math::nested_rev_autodiff;
  using stan::math::ScopedChainableStack;

  ScopedChainableStack scoped_stack;

  nested_rev_autodiff nested;

  scoped_stack.execute([] {
    gradable g_out = setup_simple();
    g_out.test();
  });

  EXPECT_EQ(stan::math::nested_size(), 0);
}

TEST_F(AgradLocalScoped, scoped_chainablestack_variadic) {
  using stan::math::nested_rev_autodiff;
  using stan::math::ScopedChainableStack;

  ScopedChainableStack scoped_stack;

  nested_rev_autodiff nested;

  scoped_stack.execute(
      [](double aval, double bval) {
        stan::math::var a = aval;
        stan::math::var b = bval;
        std::vector<stan::math::var> x{a, b};
        stan::math::var f = 2 * a * b;
        Eigen::Matrix<double, Eigen::Dynamic, 1> g_expected(2);
        g_expected << 2 * bval, 2 * aval;
        gradable g_out(x, f, g_expected);
        g_out.test();
      },
      5.0, 8.0);

  EXPECT_EQ(stan::math::nested_size(), 0);
}

TEST_F(AgradLocalScoped, scoped_chainablestack_holder) {
  using stan::math::nested_rev_autodiff;
  using stan::math::ScopedChainableStack;

  ScopedChainableStack scoped_stack;

  nested_rev_autodiff nested;

  std::unique_ptr<gradable> holder = scoped_stack.execute(
      [] { return std::make_unique<gradable>(gradable(setup_simple())); });

  EXPECT_EQ(stan::math::nested_size(), 0);

  scoped_stack.execute([&] { holder->test(); });

  EXPECT_FLOAT_EQ(holder->adj(), 1.0);

  nested.set_zero_all_adjoints();

  EXPECT_FLOAT_EQ(holder->adj(), 1.0);

  scoped_stack.execute([] { stan::math::set_zero_all_adjoints(); });

  EXPECT_FLOAT_EQ(holder->adj(), 0.0);
}

TEST_F(AgradLocalScoped, scoped_chainablestack_nesting) {
  using stan::math::nested_rev_autodiff;
  using stan::math::ScopedChainableStack;
  using stan::math::var;

  ScopedChainableStack scoped_stack;
  ScopedChainableStack nested_scoped_stack;

  EXPECT_NO_THROW(scoped_stack.execute([&] {
    double cgrad_a = nested_scoped_stack.execute([] {
      var a = 2.0;
      var b = 4.0;
      var c = a * a + b;
      c.grad();
      return a.adj();
    });
    return cgrad_a;
  }));

  EXPECT_THROW(scoped_stack.execute([&] {
    double cgrad_a = scoped_stack.execute([] {
      var a = 2.0;
      var b = 4.0;
      var c = a * a + b;
      c.grad();
      return a.adj();
    });
    return cgrad_a;
  }),
               std::logic_error);
}
