
// Copyright 2008-2009 Daniel James.
// Copyright 2022-2023 Christian Mazakas.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "../helpers/unordered.hpp"

#include "../helpers/test.hpp"
#include <boost/preprocessor/seq.hpp>
#include <list>

// TODO: this test needs to be someday cleaned up to not be so heavily
// macro-generated
//

namespace equality_tests {
  struct mod_compare
  {
    bool alt_hash_;

    explicit mod_compare(bool alt_hash = false) : alt_hash_(alt_hash) {}

    bool operator()(int x, int y) const { return x % 1000 == y % 1000; }

    std::size_t operator()(int x) const
    {
      return alt_hash_ ? static_cast<std::size_t>(x % 250)
                       : static_cast<std::size_t>((x + 5) % 250);
    }
  };

#ifdef BOOST_UNORDERED_FOA_TESTS
  using boost_unordered_set =
    boost::unordered_flat_set<int, mod_compare, mod_compare>;

  using boost_unordered_map =
    boost::unordered_flat_map<int, int, mod_compare, mod_compare>;

  using boost_unordered_node_set =
    boost::unordered_node_set<int, mod_compare, mod_compare>;

  using boost_unordered_node_map =
    boost::unordered_node_map<int, int, mod_compare, mod_compare>;

#define UNORDERED_EQUALITY_MULTISET_TEST(seq1, op, seq2)                       \
  {                                                                            \
  }

#define UNORDERED_EQUALITY_MULTIMAP_TEST(seq1, op, seq2)                       \
  {                                                                            \
  }
#else
  typedef boost::unordered_set<int, mod_compare, mod_compare>
    boost_unordered_set;

  typedef boost::unordered_map<int, int, mod_compare, mod_compare>
    boost_unordered_map;

#define UNORDERED_EQUALITY_MULTISET_TEST(seq1, op, seq2)                       \
  {                                                                            \
    boost::unordered_multiset<int, mod_compare, mod_compare> set1, set2;       \
    BOOST_PP_SEQ_FOR_EACH(UNORDERED_SET_INSERT, set1, seq1)                    \
    BOOST_PP_SEQ_FOR_EACH(UNORDERED_SET_INSERT, set2, seq2)                    \
    BOOST_TEST(set1 op set2);                                                  \
  }

#define UNORDERED_EQUALITY_MULTIMAP_TEST(seq1, op, seq2)                       \
  {                                                                            \
    boost::unordered_multimap<int, int, mod_compare, mod_compare> map1, map2;  \
    BOOST_PP_SEQ_FOR_EACH(UNORDERED_MAP_INSERT, map1, seq1)                    \
    BOOST_PP_SEQ_FOR_EACH(UNORDERED_MAP_INSERT, map2, seq2)                    \
    BOOST_TEST(map1 op map2);                                                  \
  }
#endif

#ifdef BOOST_UNORDERED_FOA_TESTS
#define UNORDERED_EQUALITY_SET_TEST(seq1, op, seq2)                            \
  {                                                                            \
    boost_unordered_set set1, set2;                                            \
    BOOST_PP_SEQ_FOR_EACH(UNORDERED_SET_INSERT, set1, seq1)                    \
    BOOST_PP_SEQ_FOR_EACH(UNORDERED_SET_INSERT, set2, seq2)                    \
    BOOST_TEST(set1 op set2);                                                  \
  }                                                                            \
  {                                                                            \
    boost_unordered_node_set set1, set2;                                       \
    BOOST_PP_SEQ_FOR_EACH(UNORDERED_SET_INSERT, set1, seq1)                    \
    BOOST_PP_SEQ_FOR_EACH(UNORDERED_SET_INSERT, set2, seq2)                    \
    BOOST_TEST(set1 op set2);                                                  \
  }
#else
#define UNORDERED_EQUALITY_SET_TEST(seq1, op, seq2)                            \
  {                                                                            \
    boost_unordered_set set1, set2;                                            \
    BOOST_PP_SEQ_FOR_EACH(UNORDERED_SET_INSERT, set1, seq1)                    \
    BOOST_PP_SEQ_FOR_EACH(UNORDERED_SET_INSERT, set2, seq2)                    \
    BOOST_TEST(set1 op set2);                                                  \
  }
#endif

#ifdef BOOST_UNORDERED_FOA_TESTS
#define UNORDERED_EQUALITY_MAP_TEST(seq1, op, seq2)                            \
  {                                                                            \
    boost_unordered_map map1, map2;                                            \
    BOOST_PP_SEQ_FOR_EACH(UNORDERED_MAP_INSERT, map1, seq1)                    \
    BOOST_PP_SEQ_FOR_EACH(UNORDERED_MAP_INSERT, map2, seq2)                    \
    BOOST_TEST(map1 op map2);                                                  \
  }                                                                            \
                                                                               \
  {                                                                            \
    boost_unordered_node_map map1, map2;                                       \
    BOOST_PP_SEQ_FOR_EACH(UNORDERED_MAP_INSERT, map1, seq1)                    \
    BOOST_PP_SEQ_FOR_EACH(UNORDERED_MAP_INSERT, map2, seq2)                    \
    BOOST_TEST(map1 op map2);                                                  \
  }
#else
#define UNORDERED_EQUALITY_MAP_TEST(seq1, op, seq2)                            \
  {                                                                            \
    boost_unordered_map map1, map2;                                            \
    BOOST_PP_SEQ_FOR_EACH(UNORDERED_MAP_INSERT, map1, seq1)                    \
    BOOST_PP_SEQ_FOR_EACH(UNORDERED_MAP_INSERT, map2, seq2)                    \
    BOOST_TEST(map1 op map2);                                                  \
  }
#endif

#define UNORDERED_SET_INSERT(r, set, item) set.insert(item);
#define UNORDERED_MAP_INSERT(r, map, item)                                     \
  map.insert(std::pair<int const, int> BOOST_PP_SEQ_TO_TUPLE(item));

  UNORDERED_AUTO_TEST (equality_size_tests) {
#ifdef BOOST_UNORDERED_FOA_TESTS
    {
      boost::unordered_flat_set<int> x1, x2;
      BOOST_TEST(x1 == x2);
      BOOST_TEST(!(x1 != x2));

      x1.insert(1);
      BOOST_TEST(x1 != x2);
      BOOST_TEST(!(x1 == x2));
      BOOST_TEST(x2 != x1);
      BOOST_TEST(!(x2 == x1));

      x2.insert(1);
      BOOST_TEST(x1 == x2);
      BOOST_TEST(!(x1 != x2));

      x2.insert(2);
      BOOST_TEST(x1 != x2);
      BOOST_TEST(!(x1 == x2));
      BOOST_TEST(x2 != x1);
      BOOST_TEST(!(x2 == x1));
    }

    {
      boost::unordered_node_set<int> x1, x2;
      BOOST_TEST(x1 == x2);
      BOOST_TEST(!(x1 != x2));

      x1.insert(1);
      BOOST_TEST(x1 != x2);
      BOOST_TEST(!(x1 == x2));
      BOOST_TEST(x2 != x1);
      BOOST_TEST(!(x2 == x1));

      x2.insert(1);
      BOOST_TEST(x1 == x2);
      BOOST_TEST(!(x1 != x2));

      x2.insert(2);
      BOOST_TEST(x1 != x2);
      BOOST_TEST(!(x1 == x2));
      BOOST_TEST(x2 != x1);
      BOOST_TEST(!(x2 == x1));
    }
#else
    boost::unordered_set<int> x1, x2;
    BOOST_TEST(x1 == x2);
    BOOST_TEST(!(x1 != x2));

    x1.insert(1);
    BOOST_TEST(x1 != x2);
    BOOST_TEST(!(x1 == x2));
    BOOST_TEST(x2 != x1);
    BOOST_TEST(!(x2 == x1));

    x2.insert(1);
    BOOST_TEST(x1 == x2);
    BOOST_TEST(!(x1 != x2));

    x2.insert(2);
    BOOST_TEST(x1 != x2);
    BOOST_TEST(!(x1 == x2));
    BOOST_TEST(x2 != x1);
    BOOST_TEST(!(x2 == x1));
#endif
  }

  UNORDERED_AUTO_TEST (equality_key_value_tests) {
    UNORDERED_EQUALITY_MULTISET_TEST((1), !=, (2))
    UNORDERED_EQUALITY_SET_TEST((2), ==, (2))
    UNORDERED_EQUALITY_MAP_TEST(((1)(1))((2)(1)), !=, ((1)(1))((3)(1)))
  }

  UNORDERED_AUTO_TEST (equality_collision_test) {
    UNORDERED_EQUALITY_MULTISET_TEST((1), !=, (501))
    UNORDERED_EQUALITY_MULTISET_TEST((1)(251), !=, (1)(501))
    UNORDERED_EQUALITY_MULTIMAP_TEST(((251)(1))((1)(1)), !=, ((501)(1))((1)(1)))
    UNORDERED_EQUALITY_MULTISET_TEST((1)(501), ==, (1)(501))
    UNORDERED_EQUALITY_SET_TEST((1)(501), ==, (501)(1))
  }

  UNORDERED_AUTO_TEST (equality_group_size_test) {
    UNORDERED_EQUALITY_MULTISET_TEST((10)(20)(20), !=, (10)(10)(20))
    UNORDERED_EQUALITY_MULTIMAP_TEST(
      ((10)(1))((20)(1))((20)(1)), !=, ((10)(1))((20)(1))((10)(1)))
    UNORDERED_EQUALITY_MULTIMAP_TEST(
      ((20)(1))((10)(1))((10)(1)), ==, ((10)(1))((20)(1))((10)(1)))
  }

  UNORDERED_AUTO_TEST (equality_map_value_test) {
    UNORDERED_EQUALITY_MAP_TEST(((1)(1)), !=, ((1)(2)))
    UNORDERED_EQUALITY_MAP_TEST(((1)(1)), ==, ((1)(1)))
    UNORDERED_EQUALITY_MULTIMAP_TEST(((1)(1)), !=, ((1)(2)))
    UNORDERED_EQUALITY_MULTIMAP_TEST(((1)(1))((1)(1)), !=, ((1)(1))((1)(2)))
    UNORDERED_EQUALITY_MULTIMAP_TEST(((1)(2))((1)(1)), ==, ((1)(1))((1)(2)))
    UNORDERED_EQUALITY_MULTIMAP_TEST(((1)(2))((1)(1)), !=, ((1)(1))((1)(3)))
  }

  UNORDERED_AUTO_TEST (equality_predicate_test) {
    UNORDERED_EQUALITY_SET_TEST((1), !=, (1001))
    UNORDERED_EQUALITY_MAP_TEST(((1)(2))((1001)(1)), !=, ((1001)(2))((1)(1)))
  }

  UNORDERED_AUTO_TEST (equality_multiple_group_test) {
    UNORDERED_EQUALITY_MULTISET_TEST(
      (1)(1)(1)(1001)(2001)(2001)(2)(1002)(3)(1003)(2003), ==,
      (3)(1003)(2003)(1002)(2)(2001)(2001)(1)(1001)(1)(1))
  }

  // Test that equality still works when the two containers have
  // different hash functions but the same equality predicate.

  UNORDERED_AUTO_TEST (equality_different_hash_test) {
#ifdef BOOST_UNORDERED_FOA_TESTS
    {
      typedef boost_unordered_set set;
      set set1(0, mod_compare(false), mod_compare(false));
      set set2(0, mod_compare(true), mod_compare(true));
      BOOST_TEST(set1 == set2);
      set1.insert(1);
      set2.insert(2);
      BOOST_TEST(set1 != set2);
      set1.insert(2);
      set2.insert(1);
      BOOST_TEST(set1 == set2);
      set1.insert(10);
      set2.insert(20);
      BOOST_TEST(set1 != set2);
      set1.insert(20);
      set2.insert(10);
      BOOST_TEST(set1 == set2);
    }

    {
      typedef boost_unordered_node_set set;
      set set1(0, mod_compare(false), mod_compare(false));
      set set2(0, mod_compare(true), mod_compare(true));
      BOOST_TEST(set1 == set2);
      set1.insert(1);
      set2.insert(2);
      BOOST_TEST(set1 != set2);
      set1.insert(2);
      set2.insert(1);
      BOOST_TEST(set1 == set2);
      set1.insert(10);
      set2.insert(20);
      BOOST_TEST(set1 != set2);
      set1.insert(20);
      set2.insert(10);
      BOOST_TEST(set1 == set2);
    }
#else
    typedef boost_unordered_set set;
    set set1(0, mod_compare(false), mod_compare(false));
    set set2(0, mod_compare(true), mod_compare(true));
    BOOST_TEST(set1 == set2);
    set1.insert(1);
    set2.insert(2);
    BOOST_TEST(set1 != set2);
    set1.insert(2);
    set2.insert(1);
    BOOST_TEST(set1 == set2);
    set1.insert(10);
    set2.insert(20);
    BOOST_TEST(set1 != set2);
    set1.insert(20);
    set2.insert(10);
    BOOST_TEST(set1 == set2);
#endif
  }
} // namespace equality_tests

RUN_TESTS()
