#if !defined(BOOST_UNORDERED_FOA_TESTS)
#error "max_load_tests is currently only supported by open-addressed containers"
#else

#include "../helpers/unordered.hpp"

#include "../helpers/helpers.hpp"
#include "../helpers/random_values.hpp"
#include "../helpers/test.hpp"
#include "../helpers/tracker.hpp"
#include "../objects/test.hpp"

template <class X> void max_load_tests(X*, test::random_generator generator)
{
  typedef typename X::size_type size_type;

  test::reset_sequence();

  X x;
  size_type max_load = x.max_load();

  BOOST_TEST_EQ(max_load, 0u);

  x.reserve(1000);
  max_load = x.max_load();

  size_type bucket_count = x.bucket_count();
  BOOST_TEST_GE(bucket_count, 1000u);

  test::ordered<X> tracker;
  {
    test::random_values<X> v(max_load, generator);

    x.insert(v.begin(), v.end());
    tracker.insert_range(v.begin(), v.end());

    BOOST_TEST_EQ(x.bucket_count(), bucket_count);
    BOOST_TEST_EQ(x.max_load(), max_load);
    BOOST_TEST_EQ(x.size(), max_load);
  }

  {
    test::random_values<X> v(100, generator);

    x.insert(v.begin(), v.end());
    tracker.insert_range(v.begin(), v.end());

    BOOST_TEST_GT(x.bucket_count(), bucket_count);
    BOOST_TEST_GT(x.max_load(), max_load);
    BOOST_TEST_GT(x.size(), max_load);
  }

  tracker.compare(x);
}

using test::default_generator;
using test::generate_collisions;
using test::limited_range;
using test::sequential;

boost::unordered_flat_set<int>* int_set_ptr;
boost::unordered_flat_map<test::movable, test::movable, test::hash,
  test::equal_to, test::allocator2<test::movable> >* test_map_ptr;

boost::unordered_flat_set<test::object, test::hash, test::equal_to,
  test::allocator1<test::object> >* test_set_tracking;
boost::unordered_flat_map<test::object, test::object, test::hash,
  test::equal_to,
  test::allocator1<std::pair<test::object const, test::object> > >*
  test_map_tracking;

UNORDERED_TEST(max_load_tests,
  ((int_set_ptr)(test_map_ptr)(test_set_tracking)(test_map_tracking))(
    (sequential)))
#endif

RUN_TESTS()
