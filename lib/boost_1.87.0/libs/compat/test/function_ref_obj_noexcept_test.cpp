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

struct F1 {
  int operator()() { return -1; }

  int operator()(int x1) noexcept { return x1; }

  int operator()(int x1, int x2) const { return 10 * x1 + x2; }

  int operator()(int x1, int x2, int x3) const noexcept { return 100 * x1 + 10 * x2 + x3; }
};

struct F2 {
  int operator()(int x1, int x2) & { return 100 * x1 + 10 * x2 + 1; }
  int operator()(int x1, int x2) const& { return 100 * x1 + 10 * x2 + 2; }
  int operator()(int x1, int x2) && { return 100 * x1 + 10 * x2 + 3; }
  int operator()(int x1, int x2) const&& { return 100 * x1 + 10 * x2 + 4; }
};

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
  {
    F1 f1;
    {
      using S1 = int() noexcept;
      using S2 = int(int) noexcept;
      using S3 = int(int, int) noexcept;
      using S4 = int(int, int, int) noexcept;

      BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<S1>, F1>));
      BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<S1>, F1&>));
      BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<S1>, F1&&>));
      BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<S1>, F1 const>));
      BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<S1>, F1 const&>));
      BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<S1>, F1 const&&>));

      compat::function_ref<S2> fv2(f1);

      BOOST_TEST_EQ(fv2(1), 1);

      BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<S2>, F1 const>));
      BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<S2>, F1 const&>));
      BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<S2>, F1 const&&>));

      BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<S3>, F1>));
      BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<S3>, F1&>));
      BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<S3>, F1&&>));
      BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<S3>, F1 const>));
      BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<S3>, F1 const&>));
      BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<S3>, F1 const&&>));

      compat::function_ref<S4> fv4(f1);

      BOOST_TEST_EQ(fv4(1, 2, 3), 123);

      BOOST_TEST_TRAIT_TRUE((std::is_constructible<compat::function_ref<S4>, F1>));
      BOOST_TEST_TRAIT_TRUE((std::is_constructible<compat::function_ref<S4>, F1&>));
      BOOST_TEST_TRAIT_TRUE((std::is_constructible<compat::function_ref<S4>, F1&&>));
      BOOST_TEST_TRAIT_TRUE((std::is_constructible<compat::function_ref<S4>, F1 const>));
      BOOST_TEST_TRAIT_TRUE((std::is_constructible<compat::function_ref<S4>, F1 const&>));
      BOOST_TEST_TRAIT_TRUE((std::is_constructible<compat::function_ref<S4>, F1 const&&>));
    }
    {
      using S1 = int() const noexcept;
      using S2 = int(int) const noexcept;
      using S3 = int(int, int) const noexcept;
      using S4 = int(int, int, int) const noexcept;

      BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<S1>, F1>));
      BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<S1>, F1&>));
      BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<S1>, F1&&>));
      BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<S1>, F1 const>));
      BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<S1>, F1 const&>));
      BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<S1>, F1 const&&>));

      BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<S2>, F1>));
      BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<S2>, F1&>));
      BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<S2>, F1&&>));
      BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<S2>, F1 const>));
      BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<S2>, F1 const&>));
      BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<S2>, F1 const&&>));

      BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<S3>, F1>));
      BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<S3>, F1&>));
      BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<S3>, F1&&>));
      BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<S3>, F1 const>));
      BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<S3>, F1 const&>));
      BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<S3>, F1 const&&>));

      compat::function_ref<S4> fv4(f1);

      BOOST_TEST_EQ(fv4(1, 2, 3), 123);

      BOOST_TEST_TRAIT_TRUE((std::is_constructible<compat::function_ref<S4>, F1>));
      BOOST_TEST_TRAIT_TRUE((std::is_constructible<compat::function_ref<S4>, F1&>));
      BOOST_TEST_TRAIT_TRUE((std::is_constructible<compat::function_ref<S4>, F1&&>));
      BOOST_TEST_TRAIT_TRUE((std::is_constructible<compat::function_ref<S4>, F1 const>));
      BOOST_TEST_TRAIT_TRUE((std::is_constructible<compat::function_ref<S4>, F1 const&>));
      BOOST_TEST_TRAIT_TRUE((std::is_constructible<compat::function_ref<S4>, F1 const&&>));
    }

    {
      using S2 = int(int) noexcept;
      using S4 = int(int, int, int) noexcept;

      auto& fref = f1;

      compat::function_ref<S2> fv2(fref);
      BOOST_TEST_EQ(fv2(1), 1);

      auto const& cfref = f1;
      compat::function_ref<S4> fv4(cfref);
      BOOST_TEST_EQ(fv4(1, 2, 3), 123);
    }

    {
      using S4 = int(int, int, int) const noexcept;

      auto const& cfref = f1;
      compat::function_ref<S4> fv4(cfref);
      BOOST_TEST_EQ(fv4(1, 2, 3), 123);
    }
  }

  {
    using S1 = int(int, int) noexcept;
    using S2 = int(int, int) const noexcept;

    BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<S1>, F2>));
    BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<S1>, F2&>));
    BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<S1>, F2&&>));
    BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<S1>, F2 const>));
    BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<S1>, F2 const&>));
    BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<S1>, F2 const&&>));

    BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<S2>, F2>));
    BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<S2>, F2&>));
    BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<S2>, F2&&>));
    BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<S2>, F2 const>));
    BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<S2>, F2 const&>));
    BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<S2>, F2 const&&>));
  }

  // invoke_r
  {
    F1 f;

    compat::function_ref<X(int, int, int) noexcept> fv1(f);
    BOOST_TEST_EQ(fv1(1, 2, 3).v, 123);
    BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<Y(int, int, int) noexcept>, F1>));
    BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<Z(int, int, int) noexcept>, F1>));

    compat::function_ref<X(int, int, int) const noexcept> fv2(f);
    BOOST_TEST_EQ(fv2(1, 2, 3).v, 123);
    BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<Y(int, int, int) const noexcept>, F1>));
    BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<Z(int, int, int) const noexcept>, F1>));

    compat::function_ref<void(int, int, int) noexcept> fv3(f);
    fv3(1, 2, 3);

    compat::function_ref<void(int, int, int) const noexcept> fv4(f);
    fv4(1, 2, 3);
  }

  // copy construct, copy assign
  {
    F1 f;

    auto id = [](int x) noexcept { return x + 7331; };

    compat::function_ref<int(int) noexcept> fv(f);
    compat::function_ref<int(int) noexcept> fv2(fv);
    BOOST_TEST_EQ(fv(12), fv2(12));

    fv2 = compat::function_ref<int(int) noexcept>(id);
    BOOST_TEST_EQ(fv2(1), 7332);

    auto add = [](int x, int y, int z) noexcept { return x + y + z; };
    compat::function_ref<int(int, int, int) const noexcept> cfv(f);
    compat::function_ref<int(int, int, int) const noexcept> cfv2(cfv);
    BOOST_TEST_EQ(cfv(1, 2, 3), cfv2(1, 2, 3));

    cfv2 = compat::function_ref<int(int, int, int) const noexcept>(add);
    BOOST_TEST_EQ(cfv2(1, 2, 3), 6);
  }

  return boost::report_errors();
}

#endif
