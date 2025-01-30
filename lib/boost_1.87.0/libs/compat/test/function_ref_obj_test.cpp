// Copyright 2024 Christian Mazakas
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

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
    F1 f;
    {
      using S1 = int();
      using S2 = int(int);
      using S3 = int(int, int);
      using S4 = int(int, int, int);

      compat::function_ref<S1> fv1(f);

      BOOST_TEST_EQ(fv1(), -1);

      BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<S1>, F1 const>));
      BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<S1>, F1 const&>));
      BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<S1>, F1 const&&>));

      compat::function_ref<S2> fv2(f);

      BOOST_TEST_EQ(fv2(1), 1);

      BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<S2>, F1 const>));
      BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<S2>, F1 const&>));

      compat::function_ref<S3> fv3(f);

      BOOST_TEST_EQ(fv3(1, 2), 12);

      BOOST_TEST_TRAIT_TRUE((std::is_constructible<compat::function_ref<S3>, F1>));
      BOOST_TEST_TRAIT_TRUE((std::is_constructible<compat::function_ref<S3>, F1&>));
      BOOST_TEST_TRAIT_TRUE((std::is_constructible<compat::function_ref<S3>, F1&&>));
      BOOST_TEST_TRAIT_TRUE((std::is_constructible<compat::function_ref<S3>, F1 const>));
      BOOST_TEST_TRAIT_TRUE((std::is_constructible<compat::function_ref<S3>, F1 const&>));
      BOOST_TEST_TRAIT_TRUE((std::is_constructible<compat::function_ref<S3>, F1 const&&>));

      compat::function_ref<S4> fv4(f);

      BOOST_TEST_EQ(fv4(1, 2, 3), 123);

      BOOST_TEST_TRAIT_TRUE((std::is_constructible<compat::function_ref<S4>, F1>));
      BOOST_TEST_TRAIT_TRUE((std::is_constructible<compat::function_ref<S4>, F1&>));
      BOOST_TEST_TRAIT_TRUE((std::is_constructible<compat::function_ref<S4>, F1&&>));
      BOOST_TEST_TRAIT_TRUE((std::is_constructible<compat::function_ref<S4>, F1 const>));
      BOOST_TEST_TRAIT_TRUE((std::is_constructible<compat::function_ref<S4>, F1 const&>));
      BOOST_TEST_TRAIT_TRUE((std::is_constructible<compat::function_ref<S4>, F1 const&&>));
    }

    {
      using S1 = int() const;
      using S2 = int(int) const;
      using S3 = int(int, int) const;
      using S4 = int(int, int, int) const;

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

      compat::function_ref<S3> fv3(f);

      BOOST_TEST_EQ(fv3(1, 2), 12);

      BOOST_TEST_TRAIT_TRUE((std::is_constructible<compat::function_ref<S3>, F1>));
      BOOST_TEST_TRAIT_TRUE((std::is_constructible<compat::function_ref<S3>, F1&>));
      BOOST_TEST_TRAIT_TRUE((std::is_constructible<compat::function_ref<S3>, F1&&>));
      BOOST_TEST_TRAIT_TRUE((std::is_constructible<compat::function_ref<S3>, F1 const>));
      BOOST_TEST_TRAIT_TRUE((std::is_constructible<compat::function_ref<S3>, F1 const&>));
      BOOST_TEST_TRAIT_TRUE((std::is_constructible<compat::function_ref<S3>, F1 const&&>));

      compat::function_ref<S4> fv4(f);

      BOOST_TEST_EQ(fv4(1, 2, 3), 123);

      BOOST_TEST_TRAIT_TRUE((std::is_constructible<compat::function_ref<S4>, F1>));
      BOOST_TEST_TRAIT_TRUE((std::is_constructible<compat::function_ref<S4>, F1&>));
      BOOST_TEST_TRAIT_TRUE((std::is_constructible<compat::function_ref<S4>, F1&&>));
      BOOST_TEST_TRAIT_TRUE((std::is_constructible<compat::function_ref<S4>, F1 const>));
      BOOST_TEST_TRAIT_TRUE((std::is_constructible<compat::function_ref<S4>, F1 const&>));
      BOOST_TEST_TRAIT_TRUE((std::is_constructible<compat::function_ref<S4>, F1 const&&>));
    }

    {
      using S1 = int();
      using S2 = int(int);

      auto& fref = f;

      compat::function_ref<S1> fv1(fref);
      BOOST_TEST_EQ(fv1(), -1);

      compat::function_ref<S2> fv2(fref);
      BOOST_TEST_EQ(fv2(1), 1);
    }

    {
      using S1 = int();
      using S2 = int(int);

      compat::function_ref<S1> fv1(std::move(f));
      BOOST_TEST_EQ(fv1(), -1);

      compat::function_ref<S2> fv2(std::move(f));
      BOOST_TEST_EQ(fv2(1), 1);
    }

    {
      using S3 = int(int, int) const;
      using S4 = int(int, int, int) const;

      auto const& fref = f;

      compat::function_ref<S3> fv3(fref);
      BOOST_TEST_EQ(fv3(1, 2), 12);

      compat::function_ref<S4> fv4(fref);
      BOOST_TEST_EQ(fv4(1, 2, 3), 123);
    }

    {
      using S3 = int(int, int) const;
      using S4 = int(int, int, int) const;

      auto const&& fref = std::move(f);

      compat::function_ref<S3> fv3(fref);
      BOOST_TEST_EQ(fv3(1, 2), 12);

      compat::function_ref<S4> fv4(fref);
      BOOST_TEST_EQ(fv4(1, 2, 3), 123);
    }
  }

  {
    F2 g;
    {
      compat::function_ref<int(int, int)> fv1(g);
      BOOST_TEST_EQ(fv1(3, 2), 321);

      auto& gref = g;
      compat::function_ref<int(int, int)> rfv1(gref);
      BOOST_TEST_EQ(rfv1(3, 2), 321);

      compat::function_ref<int(int, int)> fv2(std::move(g));
      BOOST_TEST_EQ(fv2(3, 2), 321);

      compat::function_ref<int(int, int) const> fv3(g);
      BOOST_TEST_EQ(fv3(3, 2), 322);

      auto const& gcref = g;
      compat::function_ref<int(int, int) const> crfv3(gcref);
      BOOST_TEST_EQ(fv3(3, 2), 322);

      compat::function_ref<int(int, int) const> fv4(std::move(g));
      BOOST_TEST_EQ(fv4(3, 2), 322);
    }
  }

  // invoke_r
  {
    F1 f;

    compat::function_ref<X(int, int, int)> fv1(f);
    BOOST_TEST_EQ(fv1(1, 2, 3).v, 123);
    BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<Y(int, int, int)>, F1>));
    BOOST_TEST_TRAIT_TRUE((std::is_constructible<compat::function_ref<Z(int, int, int)>, F1>));

    compat::function_ref<X(int, int, int) const> fv2(f);
    BOOST_TEST_EQ(fv2(1, 2, 3).v, 123);
    BOOST_TEST_TRAIT_FALSE((std::is_constructible<compat::function_ref<Y(int, int, int) const>, F1>));
    BOOST_TEST_TRAIT_TRUE((std::is_constructible<compat::function_ref<Z(int, int, int) const>, F1>));

    compat::function_ref<void(int, int, int)> fv3(f);
    fv3(1, 2, 3);

    compat::function_ref<void(int, int, int) const> fv4(f);
    fv4(1, 2, 3);
  }

  // copy construct, copy assign
  {
    F1 f;

    auto id = [] { return 7331; };

    compat::function_ref<int()> fv(f);
    compat::function_ref<int()> fv2(fv);
    BOOST_TEST_EQ(fv(), fv2());

    fv2 = compat::function_ref<int()>(id);
    BOOST_TEST_EQ(fv2(), 7331);

    auto add = [](int x, int y) { return x + y; };
    compat::function_ref<int(int, int) const> cfv(f);
    compat::function_ref<int(int, int) const> cfv2(cfv);
    BOOST_TEST_EQ(cfv(1, 2), cfv2(1, 2));

    cfv2 = compat::function_ref<int(int, int) const>(add);
    BOOST_TEST_EQ(cfv2(1, 2), 3);
  }

  return boost::report_errors();
}
