// Copyright 2024 Christian Mazakas
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#if defined(__GNUC__) && __GNUC__ == 7
#pragma GCC diagnostic ignored "-Wnoexcept-type"
#endif

#include <boost/compat/function_ref.hpp>

#include <boost/config.hpp>
#include <boost/config/pragma_message.hpp>
#include <boost/config/workaround.hpp>
#include <boost/core/lightweight_test.hpp>
#include <boost/core/lightweight_test_trait.hpp>

int f0() { return -1; }
int g0() { return 7331; }

int f1(int x1) noexcept { return x1; }

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

#if !defined(BOOST_NO_CXX11_HDR_TYPE_TRAITS)
  BOOST_TEST_TRAIT_TRUE((std::is_trivially_copyable<compat::function_ref<int()>>));
  BOOST_TEST_TRAIT_TRUE((std::is_trivially_copyable<compat::function_ref<int() const>>));
  BOOST_TEST_TRAIT_TRUE((std::is_trivially_copyable<compat::function_ref<int(int)>>));
  BOOST_TEST_TRAIT_TRUE((std::is_trivially_copyable<compat::function_ref<int(int) const>>));
  BOOST_TEST_TRAIT_TRUE((std::is_trivially_copyable<compat::function_ref<int(int, int)>>));
  BOOST_TEST_TRAIT_TRUE((std::is_trivially_copyable<compat::function_ref<int(int, int) const>>));
  BOOST_TEST_TRAIT_TRUE((std::is_trivially_copyable<compat::function_ref<int(int, int, int)>>));
  BOOST_TEST_TRAIT_TRUE((std::is_trivially_copyable<compat::function_ref<int(int, int, int) const>>));
  BOOST_TEST_TRAIT_TRUE(
      (std::is_trivially_copyable<compat::function_ref<int(std::unique_ptr<int>, std::unique_ptr<int>)>>));
  BOOST_TEST_TRAIT_TRUE(
      (std::is_trivially_copyable<compat::function_ref<int(std::unique_ptr<int>, std::unique_ptr<int>) const>>));
  BOOST_TEST_TRAIT_TRUE((std::is_trivially_copyable<compat::function_ref<X(int, int, int)>>));
  BOOST_TEST_TRAIT_TRUE((std::is_trivially_copyable<compat::function_ref<X(int, int, int) const>>));
  BOOST_TEST_TRAIT_TRUE((std::is_trivially_copyable<compat::function_ref<void(int, int, int)>>));
  BOOST_TEST_TRAIT_TRUE((std::is_trivially_copyable<compat::function_ref<void(int, int, int) const>>));

#if !BOOST_WORKAROUND(BOOST_MSVC, < 1910)
  struct W {
    int w_;
  };

  BOOST_TEST_TRAIT_FALSE((std::is_assignable<compat::function_ref<int()>, W>));
  BOOST_TEST_TRAIT_FALSE((std::is_assignable<compat::function_ref<int() const>, W>));

  BOOST_TEST_TRAIT_FALSE((std::is_assignable<compat::function_ref<int()>, compat::function_ref<int() const>>));
  BOOST_TEST_TRAIT_FALSE((std::is_assignable<compat::function_ref<int() const>, compat::function_ref<int()>>));

  BOOST_TEST_TRAIT_FALSE((std::is_assignable<compat::function_ref<int()>, compat::function_ref<int(int)>>));
  BOOST_TEST_TRAIT_FALSE((std::is_assignable<compat::function_ref<int(int)>, compat::function_ref<int()>>));
#endif

#else
  BOOST_PRAGMA_MESSAGE("<type_traits> is incomplete, skipping is_trivially_copyable checks")
#endif

  // f0
  {
    compat::function_ref<int()> fv1(f0);
    BOOST_TEST_EQ(fv1(), -1);

    compat::function_ref<int() const> fv2(f0);
    BOOST_TEST_EQ(fv2(), -1);
  }

  // f1
  {
    compat::function_ref<int(int)> fv1(f1);
    BOOST_TEST_EQ(fv1(1), 1);

    compat::function_ref<int(int) const> fv2(f1);
    BOOST_TEST_EQ(fv2(1), 1);
  }

  // f2
  {
    compat::function_ref<int(int, int)> fv1(f2);
    BOOST_TEST_EQ(fv1(1, 2), 12);

    compat::function_ref<int(int, int) const> fv2(f2);
    BOOST_TEST_EQ(fv2(1, 2), 12);
  }

  // f3
  {
    compat::function_ref<int(int, int, int)> fv1(f3);
    BOOST_TEST_EQ(fv1(1, 2, 3), 123);

    compat::function_ref<int(int, int, int) const> fv2(f3);
    BOOST_TEST_EQ(fv2(1, 2, 3), 123);
  }

  // g
  {
    using S1 = int(std::unique_ptr<int>, std::unique_ptr<int>);
    using S2 = int(std::unique_ptr<int>, std::unique_ptr<int>) const;

    compat::function_ref<S1> fv1(g);
    {
      auto p = std::unique_ptr<int>(new int{1});
      auto q = std::unique_ptr<int>(new int{2});
      BOOST_TEST_EQ(fv1(std::move(p), std::move(q)), 12);
    }

    compat::function_ref<S2> fv2(g);
    {
      auto p = std::unique_ptr<int>(new int{130});
      auto q = std::unique_ptr<int>(new int{37});
      BOOST_TEST_EQ(fv1(std::move(p), std::move(q)), 1337);
    }
  }

  // invoke_r
  {
    compat::function_ref<X(int, int, int)> fv1(f3);
    BOOST_TEST_EQ(fv1(1, 2, 3).v, 123);
    BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<Y(int, int, int)>, decltype(f3)>));
    BOOST_TEST_TRAIT_TRUE((std::is_constructible<compat::function_ref<Z(int, int, int)>, decltype(f3)>));

    compat::function_ref<X(int, int, int) const> fv2(f3);
    BOOST_TEST_EQ(fv2(1, 2, 3).v, 123);
    BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<Y(int, int, int) const>, decltype(f3)>));
    BOOST_TEST_TRAIT_TRUE((std::is_constructible<compat::function_ref<Z(int, int, int) const>, decltype(f3)>));

    compat::function_ref<void(int, int, int)> fv3(f3);
    fv3(1, 2, 3);

    compat::function_ref<void(int, int, int) const> fv4(f3);
    fv4(1, 2, 3);
  }

  // copy construct, copy assign
  {
    compat::function_ref<int()> fv(f0);
    compat::function_ref<int()> fv2(fv);
    BOOST_TEST_EQ(fv(), fv2());

    fv2 = g0;
    BOOST_TEST_EQ(fv2(), 7331);

    compat::function_ref<int() const> cfv(f0);
    compat::function_ref<int() const> cfv2(cfv);
    BOOST_TEST_EQ(cfv(), cfv2());

    cfv2 = g0;
    BOOST_TEST_EQ(cfv2(), 7331);
  }

  return boost::report_errors();
}
