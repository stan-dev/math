// Copyright Antony Polukhin, 2023-2024.
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <boost/stacktrace.hpp>

#include <iostream>
#include <thread>

#include <boost/core/lightweight_test.hpp>

namespace boost { namespace stacktrace { namespace impl {
  void assert_no_pending_traces() noexcept;
}}}

using boost::stacktrace::stacktrace;

struct test_no_pending_on_finish {
  ~test_no_pending_on_finish() {
    boost::stacktrace::impl::assert_no_pending_traces();
  }
};


BOOST_NOINLINE BOOST_SYMBOL_VISIBLE void in_test_throw_1(const char* msg) {
  std::string new_msg{msg};
  throw std::runtime_error(new_msg);
}

BOOST_NOINLINE BOOST_SYMBOL_VISIBLE void in_test_throw_2(const char* msg) {
  std::string new_msg{msg};
  throw std::logic_error(new_msg);
}

BOOST_NOINLINE BOOST_SYMBOL_VISIBLE void in_test_rethrow_1(const char* msg) {
  try {
    in_test_throw_1(msg);
  } catch (const std::exception&) {
    throw;
  }
}

BOOST_NOINLINE BOOST_SYMBOL_VISIBLE void in_test_rethrow_2(const char* msg) {
  try {
    in_test_throw_2(msg);
  } catch (const std::exception&) {
    try {
      in_test_throw_1(msg);
    } catch (const std::exception&) {}

    throw;
  }
}

BOOST_NOINLINE BOOST_SYMBOL_VISIBLE void test_no_exception() {
  auto trace = stacktrace::from_current_exception();
  BOOST_TEST(!trace);
}

BOOST_NOINLINE BOOST_SYMBOL_VISIBLE void test_trace_from_exception() {
  // const test_no_pending_on_finish guard{}; // something strange
  try {
    in_test_throw_1("testing basic");
  } catch (const std::exception&) {
    auto trace = stacktrace::from_current_exception();
    BOOST_TEST(trace);
    std::cout << "Tarce in test_trace_from_exception(): " << trace << '\n';
    BOOST_TEST(to_string(trace).find("in_test_throw_1") != std::string::npos);
  }
}

BOOST_NOINLINE BOOST_SYMBOL_VISIBLE void test_after_other_exception() {
  try {
    in_test_throw_1("test_other_exception_active");
  } catch (const std::exception&) {
    try {
      in_test_throw_2("test_other_exception_active 2");
    } catch (const std::exception&) {}

    auto trace = stacktrace::from_current_exception();
    BOOST_TEST(trace);
    std::cout << "Tarce in test_after_other_exception(): " << trace;
    BOOST_TEST(to_string(trace).find("in_test_throw_1") != std::string::npos);
    BOOST_TEST(to_string(trace).find("in_test_throw_2") == std::string::npos);
  }
}

BOOST_NOINLINE BOOST_SYMBOL_VISIBLE void test_rethrow() {
  try {
    in_test_rethrow_1("test rethrow");
  } catch (const std::exception&) {
    auto trace = stacktrace::from_current_exception();
    BOOST_TEST(trace);
    std::cout << "Tarce in test_rethrow(): " << trace << '\n';
    BOOST_TEST(to_string(trace).find("in_test_throw_1")   != std::string::npos);
    BOOST_TEST(to_string(trace).find("in_test_rethrow_1") != std::string::npos);
  }
}

BOOST_NOINLINE BOOST_SYMBOL_VISIBLE void test_rethrow_after_other_exception() {
  try {
    in_test_rethrow_2("test_rethrow_after_other_exception");
  } catch (const std::exception&) {
    auto trace = stacktrace::from_current_exception();
    BOOST_TEST(trace);
    std::cout << "Tarce in test_rethrow_after_other_exception(): " << trace << '\n';
    BOOST_TEST(to_string(trace).find("in_test_throw_1")   == std::string::npos);
    BOOST_TEST(to_string(trace).find("in_test_throw_2")   != std::string::npos);
    BOOST_TEST(to_string(trace).find("in_test_rethrow_2") != std::string::npos);
  }
}

BOOST_NOINLINE BOOST_SYMBOL_VISIBLE void test_nested() {
  try {
    in_test_throw_1("test_other_exception_active");
  } catch (const std::exception&) {
    try {
      in_test_throw_2("test_other_exception_active 2");
    } catch (const std::exception&) {
      auto trace = stacktrace::from_current_exception();
      BOOST_TEST(trace);
      std::cout << "Tarce in test_nested(): " << trace << '\n';
      BOOST_TEST(to_string(trace).find("in_test_throw_1") == std::string::npos);
      BOOST_TEST(to_string(trace).find("in_test_throw_2") != std::string::npos);
    }
  }
}

BOOST_NOINLINE BOOST_SYMBOL_VISIBLE void test_rethrow_nested() {
  std::exception_ptr ptr;

  try {
    in_test_throw_1("test_other_exception_active");
  } catch (const std::exception&) {
    try {
      in_test_throw_2("test_other_exception_active 2");
    } catch (const std::exception&) {
      ptr = std::current_exception();
    }
  }

  try {
    std::rethrow_exception(ptr);
  } catch (...) {
    auto trace = stacktrace::from_current_exception();
    BOOST_TEST(trace);
    std::cout << "Tarce in test_rethrow_nested(): " << trace << '\n';
    BOOST_TEST(to_string(trace).find("in_test_throw_1") == std::string::npos);
#if defined(BOOST_MSVC)
    BOOST_TEST(to_string(trace).find("in_test_throw_2") == std::string::npos);
#else
    BOOST_TEST(to_string(trace).find("in_test_throw_2") != std::string::npos);
#endif
  }
}

BOOST_NOINLINE BOOST_SYMBOL_VISIBLE void test_from_other_thread() {

// MinGW error: 'thread' is not a member of 'std'
#ifndef __MINGW32__
  std::exception_ptr ptr;

  std::thread t([&ptr]{
    try {
      in_test_throw_1("test_other_exception_active");
    } catch (const std::exception&) {
      try {
        in_test_throw_2("test_other_exception_active 2");
      } catch (const std::exception&) {
        ptr = std::current_exception();
      }
    }
  });
  t.join();

  try {
    std::rethrow_exception(ptr);
  } catch (...) {
    auto trace = stacktrace::from_current_exception();
    BOOST_TEST(trace);
    std::cout << "Tarce in test_rethrow_nested(): " << trace << '\n';
    BOOST_TEST(to_string(trace).find("in_test_throw_1") == std::string::npos);
#if defined(BOOST_MSVC)
    BOOST_TEST(to_string(trace).find("in_test_throw_2") == std::string::npos);
#else
    BOOST_TEST(to_string(trace).find("in_test_throw_2") != std::string::npos);
#endif
  }
#endif
}

int main() {
  const test_no_pending_on_finish guard{};

  BOOST_TEST(boost::stacktrace::this_thread::get_capture_stacktraces_at_throw());

  test_no_exception();
  test_trace_from_exception();
  test_after_other_exception();
  test_rethrow();
  test_rethrow_after_other_exception();
  test_nested();
  test_rethrow_nested();
  test_from_other_thread();

  return boost::report_errors();
}
