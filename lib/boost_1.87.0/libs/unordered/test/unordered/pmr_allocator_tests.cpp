//
// Copyright 2024 Braden Ganetsky.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "../helpers/helpers.hpp"
#include "../helpers/pmr.hpp"
#include "../helpers/test.hpp"
#include "../helpers/unordered.hpp"
#include <boost/config/pragma_message.hpp>
#include <string>

#ifdef BOOST_NO_CXX17_HDR_MEMORY_RESOURCE

BOOST_PRAGMA_MESSAGE(
  "Test skipped because C++17 header <memory_resource> is not available.")

#elif defined(__MAC_OS_X_VERSION_MIN_REQUIRED) &&                              \
  __MAC_OS_X_VERSION_MIN_REQUIRED < 140000

BOOST_PRAGMA_MESSAGE(
  "Test skipped because __MAC_OS_X_VERSION_MIN_REQUIRED < 140000");

#else

namespace pmr_allocator_tests {

  using pmr_string = std::basic_string<char, std::char_traits<char>,
    std::pmr::polymorphic_allocator<char> >;

#if defined(BOOST_UNORDERED_CFOA_TESTS)
  static boost::unordered::pmr::concurrent_flat_map<std::string, std::string>*
    test_string_flat_map;
  static boost::unordered::pmr::concurrent_flat_map<pmr_string, pmr_string>*
    test_pmr_string_flat_map;
  static boost::unordered::pmr::concurrent_node_map<std::string, std::string>*
    test_string_node_map;
  static boost::unordered::pmr::concurrent_node_map<pmr_string, pmr_string>*
    test_pmr_string_node_map;
  static boost::unordered::pmr::concurrent_flat_set<std::string>*
    test_string_flat_set;
  static boost::unordered::pmr::concurrent_flat_set<pmr_string>*
    test_pmr_string_flat_set;
  static boost::unordered::pmr::concurrent_node_set<std::string>*
    test_string_node_set;
  static boost::unordered::pmr::concurrent_node_set<pmr_string>*
    test_pmr_string_node_set;
#define PMR_ALLOCATOR_TESTS_ARGS                                               \
  ((test_string_flat_map)(test_pmr_string_flat_map)                            \
   (test_string_node_map)(test_pmr_string_node_map)                            \
   (test_string_flat_set)(test_pmr_string_flat_set)                            \
   (test_string_node_set)(test_pmr_string_node_set))
#elif defined(BOOST_UNORDERED_FOA_TESTS)
  static boost::unordered::pmr::unordered_flat_map<std::string, std::string>*
    test_string_flat_map;
  static boost::unordered::pmr::unordered_flat_map<pmr_string, pmr_string>*
    test_pmr_string_flat_map;
  static boost::unordered::pmr::unordered_node_map<std::string, std::string>*
    test_string_node_map;
  static boost::unordered::pmr::unordered_node_map<pmr_string, pmr_string>*
    test_pmr_string_node_map;
  static boost::unordered::pmr::unordered_flat_set<std::string>*
    test_string_flat_set;
  static boost::unordered::pmr::unordered_flat_set<pmr_string>*
    test_pmr_string_flat_set;
  static boost::unordered::pmr::unordered_node_set<std::string>*
    test_string_node_set;
  static boost::unordered::pmr::unordered_node_set<pmr_string>*
    test_pmr_string_node_set;
#define PMR_ALLOCATOR_TESTS_ARGS                                               \
  ((test_string_flat_map)(test_pmr_string_flat_map)(test_string_node_map)(test_pmr_string_node_map)(test_string_flat_set)(test_pmr_string_flat_set)(test_string_node_set)(test_pmr_string_node_set))
#else
  static boost::unordered::pmr::unordered_map<std::string, std::string>*
    test_string_map;
  static boost::unordered::pmr::unordered_map<pmr_string, pmr_string>*
    test_pmr_string_map;
  static boost::unordered::pmr::unordered_multimap<std::string, std::string>*
    test_string_multimap;
  static boost::unordered::pmr::unordered_multimap<pmr_string, pmr_string>*
    test_pmr_string_multimap;
  static boost::unordered::pmr::unordered_set<std::string>* test_string_set;
  static boost::unordered::pmr::unordered_set<pmr_string>* test_pmr_string_set;
  static boost::unordered::pmr::unordered_multiset<std::string>*
    test_string_multiset;
  static boost::unordered::pmr::unordered_multiset<pmr_string>*
    test_pmr_string_multiset;
#define PMR_ALLOCATOR_TESTS_ARGS                                               \
  ((test_string_map)(test_pmr_string_map)(test_string_multimap)(test_pmr_string_multimap)(test_string_set)(test_pmr_string_set)(test_string_multiset)(test_pmr_string_multiset))
#endif

  template <class X>
  typename std::enable_if<!test::is_map<X>::value, std::size_t>::type
  emplace_strings(X& x)
  {
    std::string_view sv =
      "this is a string that's longer than the SBO threshold";
    x.emplace(sv);
    // Return how many chars were allocated using a pmr allocator
    return std::is_same<typename X::key_type, pmr_string>::value ? sv.size() + 1
                                                                 : 0;
  }

  template <class X>
  typename std::enable_if<test::is_map<X>::value, std::size_t>::type
  emplace_strings(X& x)
  {
    std::string_view key =
      "this is a string that's longer than the SBO threshold";
    std::string_view value =
      "this is another long string that's longer than the SBO threshold";
    x.emplace(key, value);
    // Return how many chars were allocated using a pmr allocator
    return std::is_same<typename X::key_type, pmr_string>::value
             ? key.size() + value.size() + 2
             : 0;
  }

  void validate_resource(
    pmr_string const& str, test::counted_new_delete_resource const& resource)
  {
    BOOST_TEST_EQ(str.get_allocator().resource(), &resource);
  }

  void validate_resource(
    std::string const&, test::counted_new_delete_resource const&)
  {
    // Pass through
  }

  template <class X>
  typename std::enable_if<!test::is_map<X>::value>::type validate_resource(
    X& x, test::counted_new_delete_resource const& resource)
  {
#if defined(BOOST_UNORDERED_CFOA_TESTS)
    x.cvisit_all(
      [&resource](auto& element) { validate_resource(element, resource); });
#else
    for (auto& element : x) {
      validate_resource(element, resource);
    }
#endif
  }

  template <class X>
  typename std::enable_if<test::is_map<X>::value>::type validate_resource(
    X& x, test::counted_new_delete_resource const& resource)
  {
#if defined(BOOST_UNORDERED_CFOA_TESTS)
    x.cvisit_all([&resource](auto& element) {
      validate_resource(element.first, resource);
      validate_resource(element.second, resource);
    });
#else
    for (auto& element : x) {
      validate_resource(element.first, resource);
      validate_resource(element.second, resource);
    }
#endif
  }

  template <class X> static void pmr_emplace_erase(X*)
  {
    using container = X;
    using allocator_type = typename container::allocator_type;

    test::counted_new_delete_resource resource;

    {
      allocator_type alloc(&resource);
      container x(alloc);

      std::size_t num_chars = emplace_strings(x);
      BOOST_TEST_EQ(x.size(), 1);
      std::size_t count_after_emplace = resource.count();
      BOOST_TEST_GT(count_after_emplace,
        num_chars + sizeof(typename container::value_type));

      validate_resource(x, resource);

      x.clear();
      BOOST_TEST_LE(resource.count(), count_after_emplace);
    }

    BOOST_TEST_EQ(resource.count(), 0);
  }

  // clang-format off

  UNORDERED_TEST(
    pmr_emplace_erase,
    PMR_ALLOCATOR_TESTS_ARGS
  )

  // clang-format on

  enum operation
  {
    copy_op,
    move_op,
    swap_op,
  };

  template <class X> void do_operation(X& x1, X& x2, operation op)
  {
    switch (op) {
    case copy_op:
      x2 = x1;
      return;
    case move_op:
      x2 = std::move(x1);
      return;
    case swap_op:
      x1.swap(x1); // Swapping with non-equal non-pocs allocators is UB
      return;
    default:
      BOOST_TEST(false);
    }
  }

  template <class X> static void pmr_no_propagate_on_operation(X*, operation op)
  {
    using container = X;
    using allocator_type = typename container::allocator_type;
    BOOST_STATIC_ASSERT(
      !std::allocator_traits<
        allocator_type>::propagate_on_container_copy_assignment::value);
    BOOST_STATIC_ASSERT(
      !std::allocator_traits<
        allocator_type>::propagate_on_container_move_assignment::value);
    BOOST_STATIC_ASSERT(!std::allocator_traits<
                        allocator_type>::propagate_on_container_swap::value);

    test::counted_new_delete_resource resource1;
    test::counted_new_delete_resource resource2;

    allocator_type alloc1(&resource1);
    allocator_type alloc2(&resource2);

    container x1(alloc1);
    container x2(alloc2);
    bool allocators_not_equal = x1.get_allocator() != x2.get_allocator();
    BOOST_TEST(allocators_not_equal);

    emplace_strings(x1);
    do_operation(x1, x2, op);
    allocators_not_equal = x1.get_allocator() != x2.get_allocator();
    BOOST_TEST(allocators_not_equal);
  }

  // clang-format off

  UNORDERED_TEST(
    pmr_no_propagate_on_operation,
    PMR_ALLOCATOR_TESTS_ARGS
    ((copy_op)(move_op)(swap_op))
  )

  // clang-format on

} // namespace pmr_allocator_tests

#endif

RUN_TESTS()
