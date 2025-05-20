// Copyright 2024 Braden Ganetsky
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// This is a file for testing of visualizations,
// such as Visual Studio Natvis or GDB pretty printers.
// Run this test and break at the label called `break_here`.
// Inspect the variables to test correctness.

#if 0 // Change to `#if 1` to test turning off SIMD optimizations
#define BOOST_UNORDERED_DISABLE_SSE2
#define BOOST_UNORDERED_DISABLE_NEON
#endif

#define BOOST_UNORDERED_ENABLE_STATS

#include <boost/core/lightweight_test.hpp>

#include <boost/interprocess/allocators/allocator.hpp>
#include <boost/interprocess/managed_shared_memory.hpp>
#include <boost/unordered/concurrent_flat_map.hpp>
#include <boost/unordered/concurrent_flat_set.hpp>
#include <boost/unordered/concurrent_node_map.hpp>
#include <boost/unordered/concurrent_node_set.hpp>
#include <boost/unordered/unordered_flat_map.hpp>
#include <boost/unordered/unordered_flat_set.hpp>
#include <boost/unordered/unordered_map.hpp>
#include <boost/unordered/unordered_node_map.hpp>
#include <boost/unordered/unordered_node_set.hpp>
#include <boost/unordered/unordered_set.hpp>
#include <boost/uuid/random_generator.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <memory>
#include <string>
#include <typeinfo>

// Prevent any "unused" errors
template <class... Args> void use(Args&&...) {}

using map_value_type = std::pair<const std::string, int>;
using set_value_type = std::string;

template <class Tester> void visualization_test(Tester& tester)
{
  // clang-format off
    auto fca_map_ptr      = tester.template construct_map<boost::unordered_map>();
    auto fca_multimap_ptr = tester.template construct_map<boost::unordered_multimap>();
    auto fca_set_ptr      = tester.template construct_set<boost::unordered_set>();
    auto fca_multiset_ptr = tester.template construct_set<boost::unordered_multiset>();
    auto& fca_map      = *fca_map_ptr;
    auto& fca_multimap = *fca_multimap_ptr;
    auto& fca_set      = *fca_set_ptr;
    auto& fca_multiset = *fca_multiset_ptr;

    auto foa_flat_map_ptr = tester.template construct_map<boost::unordered_flat_map>();
    auto foa_flat_set_ptr = tester.template construct_set<boost::unordered_flat_set>();
    auto foa_node_map_ptr = tester.template construct_map<boost::unordered_node_map>();
    auto foa_node_set_ptr = tester.template construct_set<boost::unordered_node_set>();
    auto& foa_flat_map = *foa_flat_map_ptr;
    auto& foa_flat_set = *foa_flat_set_ptr;
    auto& foa_node_map = *foa_node_map_ptr;
    auto& foa_node_set = *foa_node_set_ptr;

    auto cfoa_flat_map_ptr = tester.template construct_map<boost::concurrent_flat_map>();
    auto cfoa_flat_set_ptr = tester.template construct_set<boost::concurrent_flat_set>();
    auto cfoa_node_map_ptr = tester.template construct_map<boost::concurrent_node_map>();
    auto cfoa_node_set_ptr = tester.template construct_set<boost::concurrent_node_set>();
    auto& cfoa_flat_map = *cfoa_flat_map_ptr;
    auto& cfoa_flat_set = *cfoa_flat_set_ptr;
    auto& cfoa_node_map = *cfoa_node_map_ptr;
    auto& cfoa_node_set = *cfoa_node_set_ptr;
  // clang-format on

  for (int i = 0; i < 5; ++i) {
    const auto str = std::to_string(i * 2);
    const auto num = i * 11;

    fca_map.emplace(str, num);
    fca_multimap.emplace(str, num);
    fca_multimap.emplace(str, num + 1);
    foa_flat_map.emplace(str, num);
    foa_node_map.emplace(str, num);
    cfoa_flat_map.emplace(str, num);
    cfoa_node_map.emplace(str, num);

    fca_set.emplace(str);
    fca_multiset.emplace(str);
    fca_multiset.emplace(str);
    foa_flat_set.emplace(str);
    foa_node_set.emplace(str);
    cfoa_flat_set.emplace(str);
    cfoa_node_set.emplace(str);
  }

  auto fca_map_begin = fca_map.begin();
  auto fca_map_end = fca_map.end();
  auto fca_multimap_begin = fca_multimap.begin();
  auto fca_multimap_end = fca_multimap.end();
  auto fca_set_begin = fca_set.begin();
  auto fca_set_end = fca_set.end();
  auto fca_multiset_begin = fca_multiset.begin();
  auto fca_multiset_end = fca_multiset.end();

  auto foa_flat_map_begin = foa_flat_map.begin();
  auto foa_flat_map_end = foa_flat_map.end();
  auto foa_flat_set_begin = foa_flat_set.begin();
  auto foa_flat_set_end = foa_flat_set.end();
  auto foa_node_map_begin = foa_node_map.begin();
  auto foa_node_map_end = foa_node_map.end();
  auto foa_node_set_begin = foa_node_set.begin();
  auto foa_node_set_end = foa_node_set.end();

  use(cfoa_flat_map, cfoa_flat_set, cfoa_node_map, cfoa_node_set);
  use(fca_map_begin, fca_map_end, fca_multimap_begin, fca_multimap_end,
    fca_set_begin, fca_set_end, fca_multiset_begin, fca_multiset_end);
  use(foa_flat_map_begin, foa_flat_map_end, foa_flat_set_begin,
    foa_flat_set_end, foa_node_map_begin, foa_node_map_end, foa_node_set_begin,
    foa_node_set_end);

  goto break_here;
break_here:;
}

class offset_ptr_tester_
{
  static constexpr std::size_t SEGMENT_SIZE = 64 * 1024;
  std::string segment_name = to_string(boost::uuids::random_generator()());
  boost::interprocess::managed_shared_memory segment{
    boost::interprocess::create_only, segment_name.c_str(), SEGMENT_SIZE};

  using map_allocator = boost::interprocess::allocator<map_value_type,
    boost::interprocess::managed_shared_memory::segment_manager>;
  using set_allocator = boost::interprocess::allocator<set_value_type,
    boost::interprocess::managed_shared_memory::segment_manager>;

  template <template <class...> class MapTemplate>
  using map_type = MapTemplate<std::string, int, boost::hash<std::string>,
    std::equal_to<std::string>, map_allocator>;
  template <template <class...> class SetTemplate>
  using set_type = SetTemplate<std::string, boost::hash<std::string>,
    std::equal_to<std::string>, set_allocator>;

public:
  offset_ptr_tester_()
  {
    boost::interprocess::shared_memory_object::remove(segment_name.c_str());
  }
  ~offset_ptr_tester_()
  {
    boost::interprocess::shared_memory_object::remove(segment_name.c_str());
  }

  template <template <class...> class MapTemplate>
  map_type<MapTemplate>* construct_map()
  {
    return &*segment.construct<map_type<MapTemplate> >(
      typeid(map_type<MapTemplate>).name())(
      map_allocator(segment.get_segment_manager()));
  }
  template <template <class...> class SetTemplate>
  set_type<SetTemplate>* construct_set()
  {
    return &*segment.construct<set_type<SetTemplate> >(
      typeid(set_type<SetTemplate>).name())(
      set_allocator(segment.get_segment_manager()));
  }
} offset_ptr_tester;

class default_tester_
{
  template <template <class...> class MapTemplate>
  using map_type = MapTemplate<std::string, int, boost::hash<std::string>,
    std::equal_to<std::string>, std::allocator<map_value_type> >;
  template <template <class...> class SetTemplate>
  using set_type = SetTemplate<std::string, boost::hash<std::string>,
    std::equal_to<std::string>, std::allocator<set_value_type> >;

public:
  template <template <class...> class MapTemplate>
  std::unique_ptr<map_type<MapTemplate> > construct_map()
  {
    using type = map_type<MapTemplate>;
    return std::unique_ptr<type>(new type());
  }
  template <template <class...> class SetTemplate>
  std::unique_ptr<set_type<SetTemplate> > construct_set()
  {
    using type = set_type<SetTemplate>;
    return std::unique_ptr<type>(new type());
  }
} default_tester;

int main()
{
  visualization_test(default_tester);
  visualization_test(offset_ptr_tester);

  return boost::report_errors();
}
