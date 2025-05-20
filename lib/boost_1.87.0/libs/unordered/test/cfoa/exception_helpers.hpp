// Copyright (C) 2023 Christian Mazakas
// Copyright (C) 2023 Joaquin M Lopez Munoz
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_UNORDERED_TEST_CFOA_EXCEPTION_HELPERS_HPP
#define BOOST_UNORDERED_TEST_CFOA_EXCEPTION_HELPERS_HPP

#include "../helpers/generators.hpp"
#include "../helpers/test.hpp"
#include "common_helpers.hpp"

#include <boost/compat/latch.hpp>
#include <boost/container_hash/hash.hpp>
#include <boost/core/span.hpp>
#include <boost/unordered/unordered_flat_map.hpp>
#include <boost/unordered/unordered_flat_set.hpp>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <condition_variable>
#include <cstddef>
#include <iostream>
#include <mutex>
#include <random>
#include <thread>
#include <type_traits>
#include <vector>

static std::size_t const num_threads =
  std::max(2u, std::thread::hardware_concurrency());

std::atomic_bool should_throw{false};

constexpr std::uint32_t throw_threshold = 2300;
constexpr std::uint32_t alloc_throw_threshold = 10;

void enable_exceptions() { should_throw = true; }
void disable_exceptions() { should_throw = false; }

struct exception_tag
{
};

struct stateful_hash
{
  int x_ = -1;

  static std::atomic<std::uint32_t> c;

  void throw_helper() const
  {
    auto n = ++c;
    if (should_throw && ((n + 1) % throw_threshold == 0)) {
      throw exception_tag{};
    }
  }

  stateful_hash() {}
  stateful_hash(stateful_hash const& rhs) : x_(rhs.x_) {}
  stateful_hash(stateful_hash&& rhs) noexcept
  {
    auto tmp = x_;
    x_ = rhs.x_;
    rhs.x_ = tmp;
  }

  stateful_hash(int const x) : x_{x} {}

  template <class T> std::size_t operator()(T const& t) const
  {
    throw_helper();
    std::size_t h = static_cast<std::size_t>(x_);
    boost::hash_combine(h, t);
    return h;
  }

  bool operator==(stateful_hash const& rhs) const { return x_ == rhs.x_; }

  friend std::ostream& operator<<(std::ostream& os, stateful_hash const& rhs)
  {
    os << "{ x_: " << rhs.x_ << " }";
    return os;
  }

  friend void swap(stateful_hash& lhs, stateful_hash& rhs) noexcept
  {
    if (&lhs != &rhs) {
      std::swap(lhs.x_, rhs.x_);
    }
  }
};

std::atomic<std::uint32_t> stateful_hash::c{0};

struct stateful_key_equal
{
  int x_ = -1;
  static std::atomic<std::uint32_t> c;

  void throw_helper() const
  {
    auto n = ++c;
    if (should_throw && ((n + 1) % throw_threshold == 0)) {
      throw exception_tag{};
    }
  }

  stateful_key_equal() = default;
  stateful_key_equal(stateful_key_equal const&) = default;
  stateful_key_equal(stateful_key_equal&& rhs) noexcept
  {
    auto tmp = x_;
    x_ = rhs.x_;
    rhs.x_ = tmp;
  }

  stateful_key_equal(int const x) : x_{x} {}

  template <class T, class U> bool operator()(T const& t, U const& u) const
  {
    throw_helper();
    return t == u;
  }

  bool operator==(stateful_key_equal const& rhs) const { return x_ == rhs.x_; }

  friend std::ostream& operator<<(
    std::ostream& os, stateful_key_equal const& rhs)
  {
    os << "{ x_: " << rhs.x_ << " }";
    return os;
  }

  friend void swap(stateful_key_equal& lhs, stateful_key_equal& rhs) noexcept
  {
    if (&lhs != &rhs) {
      std::swap(lhs.x_, rhs.x_);
    }
  }
};
std::atomic<std::uint32_t> stateful_key_equal::c{0};

static std::atomic<std::uint32_t> allocator_c = {};

template <class T> struct stateful_allocator
{
  int x_ = -1;

  void throw_helper() const
  {
    auto n = ++allocator_c;
    if (should_throw && ((n + 1) % alloc_throw_threshold == 0)) {
      throw exception_tag{};
    }
  }

  using value_type = T;

  stateful_allocator() = default;
  stateful_allocator(stateful_allocator const&) = default;
  stateful_allocator(stateful_allocator&&) = default;

  stateful_allocator(int const x) : x_{x} {}

  template <class U>
  stateful_allocator(stateful_allocator<U> const& rhs) : x_{rhs.x_}
  {
  }

  T* allocate(std::size_t n)
  {
    throw_helper();
    return static_cast<T*>(::operator new(n * sizeof(T)));
  }

  void deallocate(T* p, std::size_t) { ::operator delete(p); }

  bool operator==(stateful_allocator const& rhs) const { return x_ == rhs.x_; }
  bool operator!=(stateful_allocator const& rhs) const { return x_ != rhs.x_; }
};

struct raii
{
  static std::atomic<std::uint32_t> default_constructor;
  static std::atomic<std::uint32_t> copy_constructor;
  static std::atomic<std::uint32_t> move_constructor;
  static std::atomic<std::uint32_t> destructor;

  static std::atomic<std::uint32_t> copy_assignment;
  static std::atomic<std::uint32_t> move_assignment;

  static std::atomic<std::uint32_t> c;
  void throw_helper() const
  {
    auto n = ++c;
    if (should_throw && ((n + 1) % throw_threshold == 0)) {
      throw exception_tag{};
    }
  }

  int x_ = -1;

  raii()
  {
    throw_helper();
    ++default_constructor;
  }

  raii(int const x) : x_{x}
  {
    throw_helper();
    ++default_constructor;
  }

  raii(raii const& rhs) : x_{rhs.x_}
  {
    throw_helper();
    ++copy_constructor;
  }
  raii(raii&& rhs) noexcept : x_{rhs.x_}
  {
    rhs.x_ = -1;
    ++move_constructor;
  }
  ~raii() { ++destructor; }

  raii& operator=(raii const& rhs)
  {
    throw_helper();
    ++copy_assignment;
    if (this != &rhs) {
      x_ = rhs.x_;
    }
    return *this;
  }

  raii& operator=(raii&& rhs) noexcept
  {
    ++move_assignment;
    if (this != &rhs) {
      x_ = rhs.x_;
      rhs.x_ = -1;
    }
    return *this;
  }

  friend bool operator==(raii const& lhs, raii const& rhs)
  {
    return lhs.x_ == rhs.x_;
  }

  friend bool operator!=(raii const& lhs, raii const& rhs)
  {
    return !(lhs == rhs);
  }

  friend bool operator==(raii const& lhs, int const x) { return lhs.x_ == x; }
  friend bool operator!=(raii const& lhs, int const x)
  {
    return !(lhs.x_ == x);
  }

  friend bool operator==(int const x, raii const& rhs) { return rhs.x_ == x; }

  friend bool operator!=(int const x, raii const& rhs)
  {
    return !(rhs.x_ == x);
  }

  friend std::ostream& operator<<(std::ostream& os, raii const& rhs)
  {
    os << "{ x_: " << rhs.x_ << " }";
    return os;
  }

  friend std::ostream& operator<<(
    std::ostream& os, std::pair<raii const, raii> const& rhs)
  {
    os << "pair<" << rhs.first << ", " << rhs.second << ">";
    return os;
  }

  static void reset_counts()
  {
    default_constructor = 0;
    copy_constructor = 0;
    move_constructor = 0;
    destructor = 0;
    copy_assignment = 0;
    move_assignment = 0;
    c = 0;

    stateful_hash::c = 0;
    stateful_key_equal::c = 0;
    allocator_c = 0;
  }

  friend void swap(raii& lhs, raii& rhs) { std::swap(lhs.x_, rhs.x_); }
};

std::atomic<std::uint32_t> raii::default_constructor{0};
std::atomic<std::uint32_t> raii::copy_constructor{0};
std::atomic<std::uint32_t> raii::move_constructor{0};
std::atomic<std::uint32_t> raii::destructor{0};
std::atomic<std::uint32_t> raii::copy_assignment{0};
std::atomic<std::uint32_t> raii::move_assignment{0};
std::atomic<std::uint32_t> raii::c{0};

std::size_t hash_value(raii const& r) noexcept
{
  boost::hash<int> hasher;
  return hasher(r.x_);
}

template <typename K>
struct exception_value_generator
{
  using value_type = raii;

  value_type operator()(test::random_generator rg)
  {
    int* p = nullptr;
    int a = generate(p, rg);
    return value_type(a);
  }
};

template <typename K, typename V>
struct exception_value_generator<std::pair<K, V> >
{
  static constexpr bool const_key = std::is_const<K>::value;
  static constexpr bool const_mapped = std::is_const<V>::value;
  using value_type = std::pair<
    typename std::conditional<const_key, raii const, raii>::type,
    typename std::conditional<const_mapped, raii const, raii>::type>;

  value_type operator()(test::random_generator rg)
  {
    int* p = nullptr;
    int a = generate(p, rg);
    int b = generate(p, rg);
    return std::make_pair(raii{a}, raii{b});
  }
};

struct exception_value_type_generator_factory_type
{
  template <typename Container>
  exception_value_generator<typename Container::value_type> get()
  {
    return {}; 
  }
} exception_value_type_generator_factory;

struct exception_init_type_generator_factory_type
{
  template <typename Container>
  exception_value_generator<typename Container::init_type> get() 
  {
    return {}; 
  }
} exception_init_type_generator_factory;

struct exception_init_type_generator_type
{
  std::pair<raii, raii> operator()(test::random_generator rg)
  {
    int* p = nullptr;
    int a = generate(p, rg);
    int b = generate(p, rg);
    return std::make_pair(raii{a}, raii{b});
  }
} exception_init_type_generator;

template <class T>
std::vector<boost::span<T> > split(
  boost::span<T> s, std::size_t const nt /* num threads*/)
{
  std::vector<boost::span<T> > subslices;
  subslices.reserve(nt);

  auto a = s.size() / nt;
  auto b = a;
  if (s.size() % nt != 0) {
    ++b;
  }

  auto num_a = nt;
  auto num_b = std::size_t{0};

  if (nt * b > s.size()) {
    num_a = nt * b - s.size();
    num_b = nt - num_a;
  }

  auto sub_b = s.subspan(0, num_b * b);
  auto sub_a = s.subspan(num_b * b);

  for (std::size_t i = 0; i < num_b; ++i) {
    subslices.push_back(sub_b.subspan(i * b, b));
  }

  for (std::size_t i = 0; i < num_a; ++i) {
    auto const is_last = i == (num_a - 1);
    subslices.push_back(
      sub_a.subspan(i * a, is_last ? boost::dynamic_extent : a));
  }

  return subslices;
}

template <class T, class F> void thread_runner(std::vector<T>& values, F f)
{
  boost::compat::latch latch(static_cast<std::ptrdiff_t>(num_threads));

  std::vector<std::thread> threads;
  auto subslices = split<T>(values, num_threads);

  for (std::size_t i = 0; i < num_threads; ++i) {
    threads.emplace_back([&f, &subslices, i, &latch] {
      latch.arrive_and_wait();

      auto s = subslices[i];
      f(s);
    });
  }

  for (auto& t : threads) {
    t.join();
  }
}

template <class T> using span_value_type = typename T::value_type;

void check_raii_counts()
{
  BOOST_TEST_GT(raii::destructor, 0u);

  BOOST_TEST_EQ(
    raii::default_constructor + raii::copy_constructor + raii::move_constructor,
    raii::destructor);
}

template <class T> void shuffle_values(std::vector<T>& v)
{
  std::random_device rd;
  std::mt19937 g(rd());

  std::shuffle(v.begin(), v.end(), g);
}

template <class F>
auto make_random_values(std::size_t count, F f) -> std::vector<decltype(f())>
{
  using vector_type = std::vector<decltype(f())>;

  vector_type v;
  v.reserve(count);
  for (std::size_t i = 0; i < count; ++i) {
    v.emplace_back(f());
  }
  return v;
}

#endif // BOOST_UNORDERED_TEST_CFOA_EXCEPTION_HELPERS_HPP
