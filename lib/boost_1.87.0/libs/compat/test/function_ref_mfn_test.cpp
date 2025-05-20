// Copyright 2024 Christian Mazakas
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#if defined(__GNUC__) && __GNUC__ == 7
#pragma GCC diagnostic ignored "-Wnoexcept-type"
#endif

#include <boost/compat/function_ref.hpp>

#if !defined(BOOST_COMPAT_HAS_AUTO_NTTP)
#include <boost/config/pragma_message.hpp>
BOOST_PRAGMA_MESSAGE("no support for placeholder NTTPs detected, skipping this test")
int main() {}
#else

#include <boost/core/lightweight_test.hpp>
#include <boost/core/lightweight_test_trait.hpp>

#include <memory>
#include <type_traits>

struct F1 {
  std::unique_ptr<int> p_;
  int q_ = 1337;

  F1() : p_(new int(1)) {}

  int m1() { return *p_ + -1; }
  int m2(int x1) noexcept { return 10 * *p_ + x1; }
  int m3(int x1, int x2) const { return 100 * *p_ + 10 * x1 + x2; }
  int m4(int x1, int x2, int x3) const noexcept { return 1000 * *p_ + 100 * x1 + 10 * x2 + x3; }
};

namespace compat = boost::compat;

int main() {
  {
    F1 f1;

    compat::function_ref<int(F1&)> fn1(compat::nontype_t<&F1::m1>{});
    BOOST_TEST_EQ(fn1(f1), 0);

    compat::function_ref<int(F1&, int)> fn2(compat::nontype_t<&F1::m2>{});
    BOOST_TEST_EQ(fn2(f1, 2), 12);

    compat::function_ref<int(F1&, int, int)> fn3(compat::nontype_t<&F1::m3>{});
    BOOST_TEST_EQ(fn3(f1, 2, 3), 123);

    compat::function_ref<int(F1&, int, int, int)> fn4(compat::nontype_t<&F1::m4>{});
    BOOST_TEST_EQ(fn4(f1, 2, 3, 4), 1234);

    compat::function_ref<int(F1&)> a1(compat::nontype_t<&F1::q_>{});
    BOOST_TEST_EQ(a1(f1), 1337);
  }

  {
    F1 f1;

    compat::function_ref<int(F1&) const> fn1(compat::nontype_t<&F1::m1>{});
    BOOST_TEST_EQ(fn1(f1), 0);

    compat::function_ref<int(F1&, int) const> fn2(compat::nontype_t<&F1::m2>{});
    BOOST_TEST_EQ(fn2(f1, 2), 12);

    compat::function_ref<int(F1&, int, int) const> fn3(compat::nontype_t<&F1::m3>{});
    BOOST_TEST_EQ(fn3(f1, 2, 3), 123);

    compat::function_ref<int(F1&, int, int, int) const> fn4(compat::nontype_t<&F1::m4>{});
    BOOST_TEST_EQ(fn4(f1, 2, 3, 4), 1234);
  }

  {
    BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<int(F1 const&)>, compat::nontype_t<&F1::m1>>));
    BOOST_TEST_TRAIT_FALSE(
        (std::is_constructible<compat::function_ref<int(F1 const&, int)>, compat::nontype_t<&F1::m2>>));

    F1 f1;

    compat::function_ref<int(F1 const&, int, int) const> fn3(compat::nontype_t<&F1::m3>{});
    BOOST_TEST_EQ(fn3(f1, 2, 3), 123);

    compat::function_ref<int(F1 const&, int, int, int) const> fn4(compat::nontype_t<&F1::m4>{});
    BOOST_TEST_EQ(fn4(f1, 2, 3, 4), 1234);

    auto const& f2 = f1;

    compat::function_ref<int(F1 const&, int, int) const> fn23(compat::nontype_t<&F1::m3>{});
    BOOST_TEST_EQ(fn3(f2, 2, 3), 123);

    compat::function_ref<int(F1 const&, int, int, int) const> fn24(compat::nontype_t<&F1::m4>{});
    BOOST_TEST_EQ(fn4(f2, 2, 3, 4), 1234);
  }

  {
    BOOST_TEST_TRAIT_TRUE((std::is_constructible<compat::function_ref<int(int)>, compat::nontype_t<&F1::m2>, F1&>));
    BOOST_TEST_TRAIT_FALSE(
        (std::is_constructible<compat::function_ref<int(int)>, compat::nontype_t<&F1::m2>, F1 const&>));
    BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<int(int)>, compat::nontype_t<&F1::m2>, F1&&>));
    BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<int(int)>, compat::nontype_t<&F1::m2>, F1>));

    BOOST_TEST_TRAIT_TRUE(
        (std::is_constructible<compat::function_ref<int(int, int) const>, compat::nontype_t<&F1::m3>, F1&>));
    BOOST_TEST_TRAIT_TRUE(
        (std::is_constructible<compat::function_ref<int(int, int) const>, compat::nontype_t<&F1::m3>, F1 const&>));
    BOOST_TEST_TRAIT_FALSE(
        (std::is_constructible<compat::function_ref<int(int, int) const>, compat::nontype_t<&F1::m3>, F1&&>));
    BOOST_TEST_TRAIT_FALSE(
        (std::is_constructible<compat::function_ref<int(int, int) const>, compat::nontype_t<&F1::m3>, F1>));
  }

  {
    F1 f1;

    compat::function_ref<int()> fn1(compat::nontype_t<&F1::m1>{}, f1);
    BOOST_TEST_EQ(fn1(), 0);

    compat::function_ref<int(int)> fn2(compat::nontype_t<&F1::m2>{}, f1);
    BOOST_TEST_EQ(fn2(2), 12);

    compat::function_ref<int(int, int)> fn3(compat::nontype_t<&F1::m3>{}, f1);
    BOOST_TEST_EQ(fn3(2, 3), 123);

    compat::function_ref<int(int, int, int)> fn4(compat::nontype_t<&F1::m4>{}, f1);
    BOOST_TEST_EQ(fn4(2, 3, 4), 1234);
  }

  {
    F1 const f1;

    BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<int() const>, compat::nontype_t<&F1::m1>, F1&>));
    BOOST_TEST_TRAIT_FALSE(
        (std::is_constructible<compat::function_ref<int(int) const>, compat::nontype_t<&F1::m2>, F1&>));

    compat::function_ref<int(int, int) const> fn3(compat::nontype_t<&F1::m3>{}, f1);
    BOOST_TEST_EQ(fn3(2, 3), 123);

    compat::function_ref<int(int, int, int) const> fn4(compat::nontype_t<&F1::m4>{}, f1);
    BOOST_TEST_EQ(fn4(2, 3, 4), 1234);
  }

  {
    F1 f1;

    compat::function_ref<int()> fn1(compat::nontype_t<&F1::m1>{}, &f1);
    BOOST_TEST_EQ(fn1(), 0);

    compat::function_ref<int(int)> fn2(compat::nontype_t<&F1::m2>{}, &f1);
    BOOST_TEST_EQ(fn2(2), 12);

    compat::function_ref<int(int, int)> fn3(compat::nontype_t<&F1::m3>{}, &f1);
    BOOST_TEST_EQ(fn3(2, 3), 123);

    compat::function_ref<int(int, int, int)> fn4(compat::nontype_t<&F1::m4>{}, &f1);
    BOOST_TEST_EQ(fn4(2, 3, 4), 1234);
  }

  {
    F1 const f1;

    BOOST_TEST_TRAIT_FALSE(
        (std::is_constructible<compat::function_ref<int() const>, compat::nontype_t<&F1::m1>, F1 const*>));
    BOOST_TEST_TRAIT_FALSE(
        (std::is_constructible<compat::function_ref<int(int) const>, compat::nontype_t<&F1::m2>, F1 const*>));

    compat::function_ref<int(int, int) const> fn3(compat::nontype_t<&F1::m3>{}, &f1);
    BOOST_TEST_EQ(fn3(2, 3), 123);

    compat::function_ref<int(int, int, int) const> fn4(compat::nontype_t<&F1::m4>{}, &f1);
    BOOST_TEST_EQ(fn4(2, 3, 4), 1234);
  }

  return boost::report_errors();
}

#endif
