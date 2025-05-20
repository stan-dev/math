
// Copyright (C) 2023 Joaquin M Lopez Munoz
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "../helpers/unordered.hpp"

#include "../helpers/random_values.hpp"
#include "../objects/test.hpp"

#include <boost/config.hpp>
#include <boost/config/workaround.hpp>

#if BOOST_WORKAROUND(BOOST_MSVC, < 1910)
#pragma warning(disable:4714) /* marked as __forceinline not inlined */
#endif

namespace {

  test::seed_t initialize_seed(1002310);

  static std::size_t counted_pointer_count = 0;

  template <typename T>
  class counted_pointer {
   public:
    counted_pointer(T* p_ = nullptr) : p{p_} {
      ++counted_pointer_count;
    }
    counted_pointer(counted_pointer const& x) : p{x.p} {
      ++counted_pointer_count;
    }
    ~counted_pointer() { 
      --counted_pointer_count;
    }

    counted_pointer& operator=(counted_pointer const&) = default;

    counted_pointer& operator=(T* p_) {
      p = p_;
      return *this;
    }

    operator T*() const noexcept { return p; }

    template <typename Q = T>
    Q& operator*() const noexcept {
      return *p;
    }

    T* operator->() const noexcept { return p; }
    counted_pointer& operator++() noexcept {
      ++p;
      return *this;
    }
    counted_pointer operator++(int) noexcept {
      auto x = *this;
      ++p;
      return x;
    }
    counted_pointer& operator+=(std::ptrdiff_t n) noexcept {
      p += n;
      return *this;
    }
    counted_pointer& operator-=(std::ptrdiff_t n) noexcept {
      p -= n;
      return *this;
    }
    friend bool operator==(const counted_pointer& x, const counted_pointer& y)
    {
      return x.p == y.p;
    }
    friend bool operator!=(const counted_pointer& x, const counted_pointer& y)
    {
      return x.p != y.p;
    }
    friend bool operator<(const counted_pointer& x, const counted_pointer& y)
    {
      return x.p < y.p;
    }
    friend bool operator<=(const counted_pointer& x, const counted_pointer& y)
    {
      return x.p <= y.p;
    }
    friend bool operator>(const counted_pointer& x, const counted_pointer& y)
    {
      return x.p > y.p;
    }
    friend bool operator>=(const counted_pointer& x, const counted_pointer& y)
    {
      return x.p >= y.p;
    }

    template <typename Q = T>
    static counted_pointer pointer_to(Q& x) noexcept {
      return std::addressof(x);
    }

   private:
    T* p;
  };

  template <class T>
  struct counted_pointer_allocator {
    using value_type = T;
    using pointer = counted_pointer<T>;

    counted_pointer_allocator() = default;
    template <class U>
    counted_pointer_allocator(const counted_pointer_allocator<U>&) noexcept {}

    template <class U>
    bool operator==(const counted_pointer_allocator<U>&) const noexcept {
      return true;
    }

    template <class U>
    bool operator!=(const counted_pointer_allocator<U>&) const noexcept {
      return false;
    }

    pointer allocate(std::size_t n) const {
      return std::allocator<T>().allocate(n);
    }

    void deallocate(pointer p, std::size_t n) const noexcept {
      std::allocator<T>().deallocate(p, n);
    }
  };

  template <class T>
  void fancy_pointer_noleak_test(T*, test::random_generator const& generator)
  {
    // https://github.com/boostorg/unordered/issues/201

    auto const pointer_count = counted_pointer_count;
    {
      test::random_values<T> v(1000, generator);
      T x(v.begin(), v.end());
      (void)x.begin();
    }
    BOOST_TEST_EQ(pointer_count, counted_pointer_count);
  }

  using test::default_generator;
  using test::generate_collisions;
  using test::limited_range;

#ifdef BOOST_UNORDERED_FOA_TESTS
  boost::unordered_flat_set<test::object, test::hash, test::equal_to,
    counted_pointer_allocator<test::object> >* test_flat_set;
  boost::unordered_node_set<test::object, test::hash, test::equal_to,
    counted_pointer_allocator<test::object> >* test_node_set;
  boost::unordered_flat_map<test::object, test::object, test::hash,
    test::equal_to, counted_pointer_allocator<
      std::pair<test::object const,test::object> > >* test_flat_map;
  boost::unordered_node_map<test::object, test::object, test::hash,
    test::equal_to, counted_pointer_allocator<
      std::pair<test::object const,test::object> > >* test_node_map;

  UNORDERED_TEST(fancy_pointer_noleak_test,
    ((test_flat_set)(test_node_set)(test_flat_map)(test_node_map))
     ((default_generator)))
#else
  boost::unordered_set<test::object, test::hash, test::equal_to,
    counted_pointer_allocator<test::object> >* test_set;
  boost::unordered_multiset<test::object, test::hash, test::equal_to,
    counted_pointer_allocator<test::object> >* test_multiset;
  boost::unordered_map<test::object, test::object, test::hash,
    test::equal_to, counted_pointer_allocator<
      std::pair<test::object const,test::object> > >* test_map;
  boost::unordered_multimap<test::object, test::object, test::hash,
    test::equal_to, counted_pointer_allocator<
      std::pair<test::object const,test::object> > >* test_multimap;

  UNORDERED_TEST(fancy_pointer_noleak_test,
    ((test_set)(test_multiset)(test_map)(test_multimap))
     ((default_generator)))
#endif

} // namespace

RUN_TESTS_QUIET()
