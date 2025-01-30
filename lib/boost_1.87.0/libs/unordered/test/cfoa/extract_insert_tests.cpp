// Copyright (C) 2024 Joaquin M Lopez Munoz
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "helpers.hpp"

#include <boost/assert.hpp>
#include <boost/unordered/concurrent_node_map.hpp>
#include <boost/unordered/concurrent_node_set.hpp>

using hasher = stateful_hash;
using key_equal = stateful_key_equal;

using node_map_type = boost::unordered::concurrent_node_map<raii, raii, hasher,
  key_equal, stateful_allocator<std::pair<raii const, raii> > >;

using node_set_type = boost::unordered::concurrent_node_set<raii, hasher,
  key_equal, stateful_allocator<raii> >;

node_map_type* test_node_map;
node_set_type* test_node_set;

namespace {

  template <class X, class GF>
  void extract_insert_tests(X*, GF gen_factory)
  {
    using value_type = typename X::value_type;
    using allocator_type = typename X::allocator_type;

    // set visit is always const access
    using arg_visit_type = typename std::conditional<
      std::is_same<typename X::key_type, typename X::value_type>::value,
      typename X::value_type const,
      typename X::value_type
    >::type;

    test::random_generator rg = test::sequential;
    auto gen = gen_factory.template get<X>();
    auto values = make_random_values(1024 * 16, [&] { return gen(rg); });


    X in(0,hasher(1),key_equal(2), allocator_type(3));
    std::vector<X> out(2,in);
    for(std::size_t i = 0; i < values.size(); ++i) {
      in.insert(values[i]);
      out[i % 3 == 0? 0 : 1].insert(values[i]);
    }

    raii::reset_counts();

    thread_runner(values, [&](boost::span<value_type> s) {
      std::size_t br1 = 0, br2 = 0, br3 = 0;

      for(auto const& v: s) {
        typename X::node_type nh;

        while (nh.empty()) {
          switch (br1++ % 3) {
          case 0:
            nh = in.extract(test::get_key<X>(v));
            BOOST_ASSERT(!nh.empty());
            break;
          case 1:
            nh = in.extract_if(
              test::get_key<X>(v), [&](arg_visit_type& v2) {
                BOOST_ASSERT(test::get_key<X>(v) == test::get_key<X>(v2));
                (void)v2;
                return false;
              });
            BOOST_ASSERT(nh.empty());
            break;
          case 2: default:
            nh = in.extract_if(
              test::get_key<X>(v), [&](arg_visit_type& v2) {
                BOOST_ASSERT(test::get_key<X>(v) == test::get_key<X>(v2));
                (void)v2;
                return true;
              });
            BOOST_ASSERT(!nh.empty());
            break;
          }
        }
        BOOST_ASSERT(nh.get_allocator() == in.get_allocator());

        while (!nh.empty()) {
          auto& o = out[br2++ % out.size()];
          typename X::insert_return_type r;
          switch (br3++ % 5) {
          case 0:
            r = o.insert(std::move(nh));
            break;
          case 1: 
            r = o.insert_or_visit(
              std::move(nh), [&](arg_visit_type& v2) {
                BOOST_ASSERT(test::get_key<X>(v) == test::get_key<X>(v2));
                (void)v2;
              });
            break;
          case 2:
            r = o.insert_or_cvisit(
              std::move(nh), [&](arg_visit_type const& v2) {
                BOOST_ASSERT(test::get_key<X>(v) == test::get_key<X>(v2));
                (void)v2;
              });
            break;
          case 3:
            r = o.insert_and_visit(
              std::move(nh), 
              [&](arg_visit_type& v2) {
                BOOST_ASSERT(test::get_key<X>(v) == test::get_key<X>(v2));
                (void)v2;
              },
              [&](arg_visit_type& v2) {
                BOOST_ASSERT(test::get_key<X>(v) == test::get_key<X>(v2));
                (void)v2;
              });
            break;
          case 4: default:
            r = o.insert_and_cvisit(
              std::move(nh), 
              [&](arg_visit_type& v2) {
                BOOST_ASSERT(test::get_key<X>(v) == test::get_key<X>(v2));
                (void)v2;
              },
              [&](arg_visit_type const& v2) {
                BOOST_ASSERT(test::get_key<X>(v) == test::get_key<X>(v2));
                (void)v2;
              });
            break;
          }
          BOOST_ASSERT(r.inserted || !r.node.empty());
          nh = std::move(r.node);
        }
      }
    });

    BOOST_TEST_EQ(in.size(), 0u);
    BOOST_TEST_EQ(out[0].size() + out[1].size(), 2 * values.size());
    BOOST_TEST_EQ(raii::default_constructor, 0u);
    BOOST_TEST_EQ(raii::copy_constructor, 0u);
    BOOST_TEST_EQ(raii::move_constructor, 0u);
    BOOST_TEST_EQ(raii::destructor, 0u);
  }

  template <class X>
  void insert_empty_node_tests(X*)
  {
    using value_type = typename X::value_type;
    using node_type = typename X::node_type ;

    X x;
    {
      node_type nh;
      auto r = x.insert(std::move(nh));
      BOOST_TEST(!r.inserted);
      BOOST_TEST(r.node.empty());
    }
    {
      node_type nh;
      auto r = x.insert_or_visit(std::move(nh), [](value_type const&) {});
      BOOST_TEST(!r.inserted);
      BOOST_TEST(r.node.empty());
    }
    {
      node_type nh;
      auto r = x.insert_or_cvisit(std::move(nh), [](value_type const&) {});
      BOOST_TEST(!r.inserted);
      BOOST_TEST(r.node.empty());
    }
    {
      node_type nh;
      auto r = x.insert_and_visit(
        std::move(nh), [](value_type const&) {}, [](value_type const&) {});
      BOOST_TEST(!r.inserted);
      BOOST_TEST(r.node.empty());
    }
    {
      node_type nh;
      auto r = x.insert_and_cvisit(
        std::move(nh), [](value_type const&) {}, [](value_type const&) {});
      BOOST_TEST(!r.inserted);
      BOOST_TEST(r.node.empty());
    }
  }

} // namespace

// clang-format off
UNORDERED_TEST(
  extract_insert_tests,
  ((test_node_map)(test_node_set))
  ((value_type_generator_factory)))

UNORDERED_TEST(
  insert_empty_node_tests,
  ((test_node_map)(test_node_set)))
// clang-format on

RUN_TESTS()
