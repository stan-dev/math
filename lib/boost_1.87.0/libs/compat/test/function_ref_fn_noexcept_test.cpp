// Copyright 2024 Christian Mazakas
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/config/pragma_message.hpp>

#if !defined(__cpp_noexcept_function_type)

BOOST_PRAGMA_MESSAGE("Test skipped, __cpp_noexcept_function_type is not defined")
int main() {}

#else

#include <boost/compat/function_ref.hpp>
#include <boost/config.hpp>
#include <boost/core/lightweight_test.hpp>
#include <boost/core/lightweight_test_trait.hpp>

int f0() { return -1; }

int f1(int x1) noexcept { return x1; }
int g1(int x1) noexcept { return x1 * 2; }

int f2(int x1, int x2) { return 10 * x1 + x2; }

int f3(int x1, int x2, int x3) noexcept { return 100 * x1 + 10 * x2 + x3; }

int g(std::unique_ptr<int> p, std::unique_ptr<int> q) { return 10 * *p + *q; }

struct X {
  int v = 0;

  X() = default;
  X(int v_) noexcept : v{v_} {}
};

struct Y {
  int v = 0;

  Y() = default;
  explicit Y(int v_) : v{v_} {}
};

struct Z {
  int v = 0;

  Z() = default;
  Z(int v_) : v{v_} {}
};

namespace compat = boost::compat;

int main() {
  BOOST_TEST_TRAIT_TRUE((std::is_trivially_copyable<compat::function_ref<int() noexcept>>));
  BOOST_TEST_TRAIT_TRUE((std::is_trivially_copyable<compat::function_ref<int() const noexcept>>));
  BOOST_TEST_TRAIT_TRUE((std::is_trivially_copyable<compat::function_ref<int(int) noexcept>>));
  BOOST_TEST_TRAIT_TRUE((std::is_trivially_copyable<compat::function_ref<int(int) const noexcept>>));
  BOOST_TEST_TRAIT_TRUE((std::is_trivially_copyable<compat::function_ref<int(int, int) noexcept>>));
  BOOST_TEST_TRAIT_TRUE((std::is_trivially_copyable<compat::function_ref<int(int, int) const noexcept>>));
  BOOST_TEST_TRAIT_TRUE((std::is_trivially_copyable<compat::function_ref<int(int, int, int) noexcept>>));
  BOOST_TEST_TRAIT_TRUE((std::is_trivially_copyable<compat::function_ref<int(int, int, int) const noexcept>>));
  BOOST_TEST_TRAIT_TRUE((std::is_trivially_copyable<compat::function_ref<X(int, int, int) noexcept>>));
  BOOST_TEST_TRAIT_TRUE((std::is_trivially_copyable<compat::function_ref<X(int, int, int) const noexcept>>));
  BOOST_TEST_TRAIT_TRUE((std::is_trivially_copyable<compat::function_ref<void(int, int, int) noexcept>>));
  BOOST_TEST_TRAIT_TRUE((std::is_trivially_copyable<compat::function_ref<void(int, int, int) const noexcept>>));

  struct W {
    int w_;
  };

  BOOST_TEST_TRAIT_FALSE((std::is_assignable<compat::function_ref<int() noexcept>, W>));
  BOOST_TEST_TRAIT_FALSE((std::is_assignable<compat::function_ref<int() const noexcept>, W>));

  BOOST_TEST_TRAIT_FALSE((std::is_assignable<compat::function_ref<int() noexcept>, compat::function_ref<int()>>));
  BOOST_TEST_TRAIT_FALSE((std::is_assignable<compat::function_ref<int()>, compat::function_ref<int() noexcept>>));
  BOOST_TEST_TRAIT_FALSE(
      (std::is_assignable<compat::function_ref<int() const noexcept>, compat::function_ref<int() const>>));

  BOOST_TEST_TRAIT_FALSE(
      (std::is_assignable<compat::function_ref<int() noexcept>, compat::function_ref<int(int) noexcept>>));
  BOOST_TEST_TRAIT_FALSE(
      (std::is_assignable<compat::function_ref<int(int) noexcept>, compat::function_ref<int() noexcept>>));

  // f0
  {
    BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<int() noexcept>, decltype(f0)>));
    BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<int() const noexcept>, decltype(f0)>));
  }

  // f1
  {
    compat::function_ref<int(int) noexcept> fv1(f1);
    BOOST_TEST_EQ(fv1(1), 1);

    compat::function_ref<int(int) const noexcept> fv2(f1);
    BOOST_TEST_EQ(fv2(1), 1);
  }

  // f2
  {
    BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<int(int, int) noexcept>, decltype(f2)>));
    BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<int(int, int) const noexcept>, decltype(f2)>));
  }

  // f3
  {
    compat::function_ref<int(int, int, int) noexcept> fv1(f3);
    BOOST_TEST_EQ(fv1(1, 2, 3), 123);

    compat::function_ref<int(int, int, int) const noexcept> fv2(f3);
    BOOST_TEST_EQ(fv2(1, 2, 3), 123);
  }

  // g
  {
    BOOST_TEST_TRAIT_FALSE(
        (std::is_constructible<compat::function_ref<int(std::unique_ptr<int>, std::unique_ptr<int>) noexcept>,
                               decltype(g)>));
    BOOST_TEST_TRAIT_FALSE(
        (std::is_constructible<compat::function_ref<int(std::unique_ptr<int>, std::unique_ptr<int>) const noexcept>,
                               decltype(g)>));
  }

  // invoke_r
  {
    compat::function_ref<X(int, int, int) noexcept> fv1(f3);
    BOOST_TEST_EQ(fv1(1, 2, 3).v, 123);
    BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<Y(int, int, int) noexcept>, decltype(f3)>));
    BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<Z(int, int, int) noexcept>, decltype(f3)>));

    compat::function_ref<X(int, int, int) const noexcept> fv2(f3);
    BOOST_TEST_EQ(fv2(1, 2, 3).v, 123);
    BOOST_TEST_TRAIT_FALSE(
        (std::is_constructible<compat::function_ref<Y(int, int, int) const noexcept>, decltype(f3)>));
    BOOST_TEST_TRAIT_FALSE(
        (std::is_constructible<compat::function_ref<Z(int, int, int) const noexcept>, decltype(f3)>));

    compat::function_ref<void(int, int, int) noexcept> fv3(f3);
    fv3(1, 2, 3);

    compat::function_ref<void(int, int, int) const noexcept> fv4(f3);
    fv4(1, 2, 3);
  }

  // copy construct, copy assign
  {
    compat::function_ref<int(int) noexcept> fv(f1);
    compat::function_ref<int(int) noexcept> fv2(fv);
    BOOST_TEST_EQ(fv(42), fv2(42));

    fv2 = g1;
    BOOST_TEST_EQ(fv2(12), 24);

    compat::function_ref<int(int) const noexcept> cfv(f1);
    compat::function_ref<int(int) const noexcept> cfv2(cfv);
    BOOST_TEST_EQ(cfv(42), cfv2(42));

    cfv2 = g1;
    BOOST_TEST_EQ(cfv2(24), 48);
  }

  return boost::report_errors();
}

#endif
