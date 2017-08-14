
// Copyright 2016 Daniel James.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "../helpers/postfix.hpp"
#include "../helpers/prefix.hpp"
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

#include "../helpers/count.hpp"
#include "../helpers/helpers.hpp"
#include "../helpers/invariants.hpp"
#include "../helpers/random_values.hpp"
#include "../helpers/test.hpp"
#include "../helpers/tracker.hpp"
#include "../objects/test.hpp"
#include <boost/next_prior.hpp>

namespace merge_tests {

UNORDERED_AUTO_TEST(merge_set)
{
    boost::unordered_set<int> x;
    boost::unordered_set<int> y;

    x.merge(y);
    BOOST_TEST(x.empty());
    BOOST_TEST(y.empty());

    x.insert(10);
    x.merge(y);
    BOOST_TEST(x.size() == 1);
    BOOST_TEST(x.count(10) == 1);
    BOOST_TEST(y.empty());

    y.merge(x);
    BOOST_TEST(x.empty());
    BOOST_TEST(y.size() == 1);
    BOOST_TEST(y.count(10) == 1);

    x.insert(10);
    x.insert(50);
    y.insert(70);
    y.insert(80);
    x.merge(y);
    BOOST_TEST_EQ(x.size(), 4u);
    BOOST_TEST_EQ(y.size(), 1u);
    BOOST_TEST_EQ(x.count(10), 1u);
    BOOST_TEST_EQ(x.count(50), 1u);
    BOOST_TEST_EQ(x.count(70), 1u);
    BOOST_TEST_EQ(x.count(80), 1u);
    BOOST_TEST_EQ(y.count(10), 1u);
    BOOST_TEST_EQ(y.count(50), 0u);
    BOOST_TEST_EQ(y.count(70), 0u);
    BOOST_TEST_EQ(y.count(80), 0u);

    test::check_equivalent_keys(x);
    test::check_equivalent_keys(y);
}

UNORDERED_AUTO_TEST(merge_multiset)
{
    boost::unordered_multiset<int> x;
    boost::unordered_multiset<int> y;

    x.merge(y);
    BOOST_TEST(x.empty());
    BOOST_TEST(y.empty());

    x.insert(10);
    x.merge(y);
    BOOST_TEST(x.size() == 1);
    BOOST_TEST(x.count(10) == 1);
    BOOST_TEST(y.empty());

    y.merge(x);
    BOOST_TEST(x.empty());
    BOOST_TEST(y.size() == 1);
    BOOST_TEST(y.count(10) == 1);

    x.insert(10);
    x.insert(50);
    y.insert(70);
    y.insert(80);
    x.merge(y);
    BOOST_TEST_EQ(x.size(), 5u);
    BOOST_TEST_EQ(y.size(), 0u);
    BOOST_TEST_EQ(x.count(10), 2u);
    BOOST_TEST_EQ(x.count(50), 1u);
    BOOST_TEST_EQ(x.count(70), 1u);
    BOOST_TEST_EQ(x.count(80), 1u);
    BOOST_TEST_EQ(y.count(10), 0u);
    BOOST_TEST_EQ(y.count(50), 0u);
    BOOST_TEST_EQ(y.count(70), 0u);
    BOOST_TEST_EQ(y.count(80), 0u);

    test::check_equivalent_keys(x);
    test::check_equivalent_keys(y);
}

#if BOOST_UNORDERED_INTEROPERABLE_NODES
UNORDERED_AUTO_TEST(merge_set_and_multiset)
{
    boost::unordered_set<int> x;
    boost::unordered_multiset<int> y;

    x.merge(y);
    BOOST_TEST(x.empty());
    BOOST_TEST(y.empty());

    x.insert(10);
    x.merge(y);
    BOOST_TEST(x.size() == 1);
    BOOST_TEST(x.count(10) == 1);
    BOOST_TEST(y.empty());

    y.merge(x);
    BOOST_TEST(x.empty());
    BOOST_TEST(y.size() == 1);
    BOOST_TEST(y.count(10) == 1);

    x.insert(10);
    x.insert(50);
    y.insert(70);
    y.insert(80);
    x.merge(y);
    BOOST_TEST_EQ(x.size(), 4u);
    BOOST_TEST_EQ(y.size(), 1u);
    BOOST_TEST_EQ(x.count(10), 1u);
    BOOST_TEST_EQ(x.count(50), 1u);
    BOOST_TEST_EQ(x.count(70), 1u);
    BOOST_TEST_EQ(x.count(80), 1u);
    BOOST_TEST_EQ(y.count(10), 1u);
    BOOST_TEST_EQ(y.count(50), 0u);
    BOOST_TEST_EQ(y.count(70), 0u);
    BOOST_TEST_EQ(y.count(80), 0u);

    test::check_equivalent_keys(x);
    test::check_equivalent_keys(y);
}
#endif

template <class X> void merge_empty_test(X*, test::random_generator generator)
{
    test::check_instances check_;

    test::random_values<X> v(1000, generator);
    X x1(v.begin(), v.end()), x2;
    x1.merge(x2);
    test::check_container(x1, v);
    BOOST_TEST(x2.empty());
    test::check_equivalent_keys(x1);
    test::check_equivalent_keys(x2);
}

template <class X>
void merge_into_empty_test(X*, test::random_generator generator)
{
    test::check_instances check_;

    test::random_values<X> v(1000, generator);
    X x1, x2(v.begin(), v.end());
    x1.merge(x2);
    test::check_container(x1, v);
    BOOST_TEST(x2.empty());
    test::check_equivalent_keys(x1);
    test::check_equivalent_keys(x2);
}

template <class X> void unique_merge_test(X*, test::random_generator generator)
{
    test::check_instances check_;

    test::random_values<X> v1(1000, generator);
    test::random_values<X> v2(1000, generator);
    v1.insert(v2.begin(), boost::next(v2.begin(), 100));
    v2.insert(v1.begin(), boost::next(v1.begin(), 100));

    X x1(v1.begin(), v1.end()), x2(v2.begin(), v2.end());
    x1.merge(x2);

    test::ordered<X> tracker1 = test::create_ordered(x1);
    test::ordered<X> tracker2 = test::create_ordered(x2);
    test::ordered<X> tracker_tmp = test::create_ordered(x2);
    tracker1.insert(v1.begin(), v1.end());
    tracker_tmp.insert(v2.begin(), v2.end());
    for (BOOST_DEDUCED_TYPENAME test::ordered<X>::iterator it =
             tracker_tmp.begin();
         it != tracker_tmp.end(); ++it) {
        if (!tracker1.insert(*it).second) {
            tracker2.insert(*it);
        }
    }

    tracker1.compare(x1);
    tracker2.compare(x2);
    test::check_equivalent_keys(x1);
    test::check_equivalent_keys(x2);
}

template <class X> void equiv_merge_test(X*, test::random_generator generator)
{
    test::check_instances check_;

    test::random_values<X> v1(1000, generator);
    test::random_values<X> v2(1000, generator);
    v1.insert(v2.begin(), boost::next(v2.begin(), 100));
    v2.insert(v1.begin(), boost::next(v1.begin(), 100));

    X x1(v1.begin(), v1.end()), x2(v2.begin(), v2.end());
    x1.merge(x2);

    test::ordered<X> tracker1 = test::create_ordered(x1);
    tracker1.insert(v1.begin(), v1.end());
    tracker1.insert(v2.begin(), v2.end());

    tracker1.compare(x1);
    BOOST_TEST(x2.empty());
    test::check_equivalent_keys(x1);
    test::check_equivalent_keys(x2);
}

boost::unordered_set<test::movable, test::hash, test::equal_to,
    std::allocator<test::movable> >* test_set_std_alloc;
boost::unordered_multimap<test::object, test::object, test::hash,
    test::equal_to, std::allocator<test::object> >* test_multimap_std_alloc;

boost::unordered_set<test::object, test::hash, test::equal_to,
    test::allocator1<test::object> >* test_set;
boost::unordered_multiset<test::movable, test::hash, test::equal_to,
    test::allocator2<test::movable> >* test_multiset;
boost::unordered_map<test::movable, test::movable, test::hash, test::equal_to,
    test::allocator2<test::movable> >* test_map;
boost::unordered_multimap<test::object, test::object, test::hash,
    test::equal_to, test::allocator1<test::object> >* test_multimap;

using test::default_generator;
using test::generate_collisions;

UNORDERED_TEST(merge_empty_test,
    ((test_set_std_alloc)(test_multimap_std_alloc)(test_set)(test_multiset)(
        test_map)(test_multimap))((default_generator)(generate_collisions)))

UNORDERED_TEST(merge_into_empty_test,
    ((test_set_std_alloc)(test_multimap_std_alloc)(test_set)(test_multiset)(
        test_map)(test_multimap))((default_generator)(generate_collisions)))

UNORDERED_TEST(unique_merge_test,
    ((test_set_std_alloc)(test_set)(test_map))((default_generator)))

UNORDERED_TEST(equiv_merge_test, ((test_multimap_std_alloc)(test_multiset)(
                                     test_multimap))((default_generator)))
}

RUN_TESTS()
