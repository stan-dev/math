//
// Copyright 2023-2024 Braden Ganetsky.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "../helpers/unordered.hpp"

#include "../helpers/count.hpp"
#include "../helpers/test.hpp"

namespace emplace_smf_tests {
  using test::smf_count;
  using test::smf_counted_object;

  using converting_key = smf_counted_object<struct cvt_key_tag_>;
  using converting_value = smf_counted_object<struct cvt_value_tag_>;

  class counted_key : public smf_counted_object<struct key_tag_>
  {
  public:
    using smf_counted_object::smf_counted_object;
    counted_key() = default;
    counted_key(const converting_key& k) : counted_key(k.index_) {}
  };
  class counted_value : public smf_counted_object<struct value_tag_>
  {
  public:
    using smf_counted_object::smf_counted_object;
    counted_value() = default;
    counted_value(const converting_value& v) : counted_value(v.index_) {}
  };

  class immovable_key : public smf_counted_object<struct imm_key_tag_>
  {
  public:
    using smf_counted_object::smf_counted_object;
    immovable_key(immovable_key&&) = delete;
    immovable_key& operator=(immovable_key&&) = delete;
  };
  class immovable_value : public smf_counted_object<struct imm_value_tag_>
  {
  public:
    using smf_counted_object::smf_counted_object;
    immovable_value(immovable_value&&) = delete;
    immovable_value& operator=(immovable_value&&) = delete;
  };

  void reset_counts()
  {
    counted_key::reset_count();
    counted_value::reset_count();
    converting_key::reset_count();
    converting_value::reset_count();
    immovable_key::reset_count();
    immovable_value::reset_count();
  }

#ifdef BOOST_UNORDERED_FOA_TESTS
  static boost::unordered_flat_map<counted_key, counted_value>* test_smf_map;
  static boost::unordered_node_map<counted_key, counted_value>*
    test_smf_node_map;
  static boost::unordered_flat_set<counted_value>* test_smf_set;
  static boost::unordered_node_set<counted_value>* test_smf_node_set;
#define EMPLACE_SMF_TESTS_MAP_ARGS ((test_smf_map)(test_smf_node_map))
#define EMPLACE_SMF_TESTS_SET_ARGS ((test_smf_set)(test_smf_node_set))
#else
  static boost::unordered_map<counted_key, counted_value>* test_smf_map;
  static boost::unordered_set<counted_value>* test_smf_set;
#define EMPLACE_SMF_TESTS_MAP_ARGS ((test_smf_map))
#define EMPLACE_SMF_TESTS_SET_ARGS ((test_smf_set))
#endif

  template <class X> static void emplace_smf_value_type_map(X*)
  {
    using container = X;
    using value_type = typename container::value_type;

    container x;

    {
      value_type val{counted_key{}, counted_value{}};
      x.clear();
      reset_counts();
      x.emplace(val);
      BOOST_TEST_EQ(counted_key::count, (smf_count{0, 1, 0, 0, 0, 0}));
      BOOST_TEST_EQ(counted_value::count, (smf_count{0, 1, 0, 0, 0, 0}));
    }

    {
      value_type val{counted_key{}, counted_value{}};
      x.clear();
      reset_counts();
      x.emplace(std::move(val));
      BOOST_TEST_EQ(counted_key::count, (smf_count{0, 1, 0, 0, 0, 0}));
      BOOST_TEST_EQ(counted_value::count, (smf_count{0, 0, 1, 0, 0, 0}));
    }

    {
      x.clear();
      reset_counts();
      x.emplace(value_type{counted_key{}, counted_value{}});
      BOOST_TEST_EQ(counted_key::count, (smf_count{1, 1, 1, 0, 0, 2}));
      BOOST_TEST_EQ(counted_value::count, (smf_count{1, 0, 2, 0, 0, 2}));
    }

    {
      counted_key key{};
      counted_value value{};
      x.clear();
      reset_counts();
      x.emplace(value_type{std::move(key), std::move(value)});
      BOOST_TEST_EQ(counted_key::count, (smf_count{0, 1, 1, 0, 0, 1}));
      BOOST_TEST_EQ(counted_value::count, (smf_count{0, 0, 2, 0, 0, 1}));
    }
  }

  UNORDERED_TEST(emplace_smf_value_type_map, EMPLACE_SMF_TESTS_MAP_ARGS)

  template <class X> static void emplace_smf_init_type_map(X*)
  {
    using container = X;
#ifdef BOOST_UNORDERED_FOA_TESTS
    using init_type = typename container::init_type;
#else
    using raw_key =
      typename std::remove_const<typename container::key_type>::type;
    using init_type = std::pair<raw_key, typename container::mapped_type>;
#endif

    container x;

    {
      init_type val{counted_key{}, counted_value{}};
      x.clear();
      reset_counts();
      x.emplace(val);
      BOOST_TEST_EQ(counted_key::count, (smf_count{0, 1, 0, 0, 0, 0}));
      BOOST_TEST_EQ(counted_value::count, (smf_count{0, 1, 0, 0, 0, 0}));
    }

    {
      init_type val{counted_key{}, counted_value{}};
      x.clear();
      reset_counts();
      x.emplace(std::move(val));
      BOOST_TEST_EQ(counted_key::count, (smf_count{0, 0, 1, 0, 0, 0}));
      BOOST_TEST_EQ(counted_value::count, (smf_count{0, 0, 1, 0, 0, 0}));
    }

    {
      x.clear();
      reset_counts();
      x.emplace(init_type{counted_key{}, counted_value{}});
      BOOST_TEST_EQ(counted_key::count, (smf_count{1, 0, 2, 0, 0, 2}));
      BOOST_TEST_EQ(counted_value::count, (smf_count{1, 0, 2, 0, 0, 2}));
    }

    {
      counted_key key{};
      counted_value value{};
      x.clear();
      reset_counts();
      x.emplace(init_type{std::move(key), std::move(value)});
      BOOST_TEST_EQ(counted_key::count, (smf_count{0, 0, 2, 0, 0, 1}));
      BOOST_TEST_EQ(counted_value::count, (smf_count{0, 0, 2, 0, 0, 1}));
    }
  }

  UNORDERED_TEST(emplace_smf_init_type_map, EMPLACE_SMF_TESTS_MAP_ARGS)

  template <class X> static void emplace_smf_value_type_set(X*)
  {
    using container = X;
    using value_type = typename container::value_type;
#ifdef BOOST_UNORDERED_FOA_TESTS
    BOOST_STATIC_ASSERT(
      std::is_same<value_type, typename container::init_type>::value);
#endif
    BOOST_STATIC_ASSERT(std::is_same<value_type, counted_value>::value);

    container x;

    {
      counted_value val{};
      x.clear();
      reset_counts();
      x.emplace(val);
      BOOST_TEST_EQ(counted_value::count, (smf_count{0, 1, 0, 0, 0, 0}));
    }

    {
      counted_value val{};
      x.clear();
      reset_counts();
      x.emplace(std::move(val));
      BOOST_TEST_EQ(counted_value::count, (smf_count{0, 0, 1, 0, 0, 0}));
    }

    {
      x.clear();
      reset_counts();
      x.emplace(counted_value{});
      BOOST_TEST_EQ(counted_value::count, (smf_count{1, 0, 1, 0, 0, 1}));
    }
  }

  UNORDERED_TEST(emplace_smf_value_type_set, EMPLACE_SMF_TESTS_SET_ARGS)

  enum emplace_kind
  {
    copy,
    move
  };

  enum emplace_status
  {
    fail,
    success
  };

  struct counted_key_checker_type
  {
    using key_type = counted_key;
    void operator()(emplace_kind kind, emplace_status status)
    {
      int copies = (kind == copy && status == success) ? 1 : 0;
      int moves = (kind == move && status == success) ? 1 : 0;
      BOOST_TEST_EQ(counted_key::count, (smf_count{0, copies, moves, 0, 0, 0}));
    }
  } counted_key_checker;

  struct converting_key_checker_type
  {
    using key_type = converting_key;
    void operator()(emplace_kind, emplace_status status)
    {
      int moves = (status == success) ? 1 : 0;
      BOOST_TEST_EQ(counted_key::count, (smf_count{1, 0, moves, 0, 0, 1}));
    }
  } converting_key_checker;

  struct counted_value_checker_type
  {
    using mapped_type = counted_value;
    void operator()(emplace_kind kind, emplace_status status)
    {
      int copies = (kind == copy && status == success) ? 1 : 0;
      int moves = (kind == move && status == success) ? 1 : 0;
      BOOST_TEST_EQ(
        counted_value::count, (smf_count{0, copies, moves, 0, 0, 0}));
    }
  } counted_value_checker;

  struct converting_value_checker_type
  {
    using mapped_type = converting_value;
    void operator()(emplace_kind, emplace_status status)
    {
      int ctors = (status == success) ? 1 : 0;
      BOOST_TEST_EQ(counted_value::count, (smf_count{ctors, 0, 0, 0, 0, 0}));
    }
  } converting_value_checker;

  template <class X, class KC, class VC>
  static void emplace_smf_key_value_map(
    X*, emplace_kind kind, KC key_checker, VC value_checker)
  {
    using container = X;
    using key_type = typename KC::key_type;
    using mapped_type = typename VC::mapped_type;

    container x;
    key_type key{};
    key_type key2 = key;
    mapped_type value{};
    mapped_type value2 = value;

    {
      reset_counts();
      auto ret = (kind == copy) ? x.emplace(key, value)
                                : x.emplace(std::move(key), std::move(value));
      BOOST_TEST_EQ(ret.second, true);
      key_checker(kind, success);
      value_checker(kind, success);
      BOOST_TEST_EQ(converting_key::count, (smf_count{0, 0, 0, 0, 0, 0}));
      BOOST_TEST_EQ(converting_value::count, (smf_count{0, 0, 0, 0, 0, 0}));
    }

    {
      reset_counts();
      auto ret = x.emplace(key2, value2);
      BOOST_TEST_EQ(ret.second, false);
      key_checker(kind, fail);
      value_checker(kind, fail);
      BOOST_TEST_EQ(converting_key::count, (smf_count{0, 0, 0, 0, 0, 0}));
      BOOST_TEST_EQ(converting_value::count, (smf_count{0, 0, 0, 0, 0, 0}));
    }

    {
      reset_counts();
      auto ret = x.emplace(std::move(key2), std::move(value2));
      BOOST_TEST_EQ(ret.second, false);
      key_checker(kind, fail);
      value_checker(kind, fail);
      BOOST_TEST_EQ(converting_key::count, (smf_count{0, 0, 0, 0, 0, 0}));
      BOOST_TEST_EQ(converting_value::count, (smf_count{0, 0, 0, 0, 0, 0}));
    }
  }

  // clang-format off

  UNORDERED_TEST(
    emplace_smf_key_value_map,
    EMPLACE_SMF_TESTS_MAP_ARGS
    ((copy)(move))
    ((counted_key_checker)(converting_key_checker))
    ((counted_value_checker)(converting_value_checker))
  )

  // clang-format on

  template <class X> static void emplace_smf_key_value_map_immovable_key(X*)
  {
#ifndef BOOST_UNORDERED_NO_INIT_TYPE_TESTS
    using container = X;
    BOOST_STATIC_ASSERT(
      std::is_same<immovable_key, typename X::key_type>::value);
    using mapped_type = typename X::mapped_type;

    container x;

    {
      reset_counts();
      auto ret = x.emplace(0, 0);
      BOOST_TEST_EQ(ret.second, true);
      BOOST_TEST_EQ(immovable_key::count, (smf_count{1, 0, 0, 0, 0, 0}));
      BOOST_TEST_EQ(mapped_type::count, (smf_count{1, 0, 0, 0, 0, 0}));
    }

    {
      reset_counts();
      auto ret = x.emplace(0, 1);
      BOOST_TEST_EQ(ret.second, false);
      BOOST_TEST_EQ(immovable_key::count, (smf_count{1, 0, 0, 0, 0, 1}));
      BOOST_TEST_EQ(mapped_type::count, (smf_count{1, 0, 0, 0, 0, 1}));
    }
#endif
  }

#ifdef BOOST_UNORDERED_FOA_TESTS
  static boost::unordered_node_map<immovable_key, counted_value>*
    test_smf_node_map_immovable_key_counted_value;
  static boost::unordered_node_map<immovable_key, immovable_value>*
    test_smf_node_map_immovable_key_immovable_value;
#define EMPLACE_SMF_TESTS_MAP_IMMOVABLE_ARGS                                   \
  ((test_smf_node_map_immovable_key_counted_value)(test_smf_node_map_immovable_key_immovable_value))
#else
  static boost::unordered_map<immovable_key, counted_value>*
    test_smf_map_immovable_key_counted_value;
  static boost::unordered_map<immovable_key, immovable_value>*
    test_smf_map_immovable_key_immovable_value;
#define EMPLACE_SMF_TESTS_MAP_IMMOVABLE_ARGS                                   \
  ((test_smf_map_immovable_key_counted_value)(test_smf_map_immovable_key_immovable_value))
#endif

  // clang-format off

  UNORDERED_TEST(
    emplace_smf_key_value_map_immovable_key,
    EMPLACE_SMF_TESTS_MAP_IMMOVABLE_ARGS
  )

  // clang-format on

} // namespace emplace_smf_tests

RUN_TESTS()
