// Copyright 2023 Christian Mazakas.
// Copyright 2023 Joaquin M Lopez Munoz.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/config.hpp>

#if defined(BOOST_CLANG_VERSION) && BOOST_CLANG_VERSION < 30900
#include <boost/config/pragma_message.hpp>
BOOST_PRAGMA_MESSAGE(
  "This version of clang is incompatible with Boost.Process");
int main() {}
#else
#include "../helpers/test.hpp"

#include <boost/unordered/unordered_flat_map.hpp>
#include <boost/unordered/unordered_map.hpp>
#include <boost/unordered/unordered_node_map.hpp>

#include <boost/unordered/unordered_flat_set.hpp>
#include <boost/unordered/unordered_node_set.hpp>
#include <boost/unordered/unordered_set.hpp>

#include <boost/unordered/concurrent_flat_map.hpp>
#include <boost/unordered/concurrent_flat_set.hpp>

#include <boost/interprocess/allocators/allocator.hpp>
#include <boost/interprocess/containers/string.hpp>
#include <boost/interprocess/managed_shared_memory.hpp>

#include <boost/process/child.hpp>
#include <boost/process/filesystem.hpp>

#include <boost/uuid/random_generator.hpp>
#include <boost/uuid/uuid_io.hpp>

#include <algorithm>
#include <iostream>
#include <type_traits>
#include <vector>

#ifndef BOOST_UNORDERED_FOA_MMAP_MAP_TYPE
#error "this requires a class template be passed as a macro"
#endif

using char_allocator = boost::interprocess::allocator<char,
  boost::interprocess::managed_shared_memory::segment_manager>;

using string_type = boost::interprocess::basic_string<char,
  std::char_traits<char>, char_allocator>;

using pair_type = std::pair<const string_type, string_type>;

using string_pair_type = std::pair<string_type, string_type>;

using string_pair_allocator = boost::interprocess::allocator<string_pair_type,
  boost::interprocess::managed_shared_memory::segment_manager>;

using pair_allocator = boost::interprocess::allocator<pair_type,
  boost::interprocess::managed_shared_memory::segment_manager>;

template <template <class, class, class, class, class> class Map,
  class MapType = Map<string_type, string_type, boost::hash<string_type>,
    std::equal_to<string_type>, pair_allocator> >
typename std::enable_if<
  !std::is_same<typename MapType::value_type, string_pair_type>::value,
  MapType>::type
get_container_type()
{
  return {};
}

template <template <class, class, class, class> class Set,
  class SetType = Set<string_pair_type, boost::hash<string_pair_type>,
    std::equal_to<string_pair_type>, string_pair_allocator> >
typename std::enable_if<
  std::is_same<typename SetType::value_type, string_pair_type>::value,
  SetType>::type
get_container_type()
{
  return {};
}

using concurrent_map = decltype(
  get_container_type<boost::concurrent_flat_map>());

using concurrent_set = decltype(
  get_container_type<boost::concurrent_flat_set>());

template <class C>
struct is_concurrent_container: std::false_type {};

template <typename... Args>
struct is_concurrent_container<boost::concurrent_flat_map<Args...> >:
  std::true_type {};

template <typename... Args>
struct is_concurrent_container<boost::concurrent_flat_set<Args...> >:
  std::true_type {};

static char const* shm_map_name = "shared_map";

template <class C>
typename std::enable_if<!is_concurrent_container<C>::value, void>::type
parent(std::string const& shm_name_, char const* exe_name, C*)
{
  struct shm_remove
  {
    char const* shm_name;

    shm_remove(char const* shm_name_) : shm_name(shm_name_)
    {
      boost::interprocess::shared_memory_object::remove(shm_name);
    }
    ~shm_remove()
    {
      boost::interprocess::shared_memory_object::remove(shm_name);
    }
  } remover{shm_name_.c_str()};

  using container_type = C;

  std::size_t const shm_size = 64 * 1024;

  auto shm_name = remover.shm_name;
  boost::interprocess::managed_shared_memory segment(
    boost::interprocess::create_only, shm_name, shm_size);

  auto segment_mngr = segment.get_segment_manager();
  char_allocator char_alloc(segment_mngr);
  pair_allocator pair_alloc(segment_mngr);

  using iterator_type = typename C::iterator;

  container_type* c =
    segment.construct<container_type>(shm_map_name)(pair_alloc);

  iterator_type* it = segment.construct<iterator_type>("shared_iterator")();

  auto const old_bc = c->bucket_count();

  BOOST_TEST(c->empty());

  boost::process::child child(exe_name, shm_name);
  child.wait();
  int ret = child.exit_code();

  BOOST_TEST(*it == c->begin());
  BOOST_TEST((**it) == *(c->begin()));

  c->erase(*it);

  auto inputs =
    std::vector<std::pair<std::string, std::string> >{{"hello", "world"}};

  for (const auto& kvp : inputs) {
    auto const& key = kvp.first;
    auto const& value = kvp.second;
    c->emplace(string_type(key.data(), char_alloc),
      string_type(value.data(), char_alloc));
  }

  BOOST_TEST_EQ(ret, 0);
  BOOST_TEST_GT(c->size(), inputs.size());
  BOOST_TEST_GT(c->bucket_count(), old_bc);

  segment.destroy<iterator_type>("shared_iterator");
  segment.destroy<container_type>(shm_map_name);
}

template <class C>
typename std::enable_if<!is_concurrent_container<C>::value, void>::type
child(std::string const& shm_name, C*)
{
  using container_type = C;
  using iterator = typename container_type::iterator;

  boost::interprocess::managed_shared_memory segment(
    boost::interprocess::open_only, shm_name.c_str());

  container_type* c = segment.find<container_type>(shm_map_name).first;

  iterator* it = segment.find<iterator>("shared_iterator").first;

  BOOST_TEST(c->empty());

  std::vector<std::pair<std::string, std::string> > inputs = {
    {"aaa", "AAA"}, {"bbb", "BBB"}, {"ccc", "CCCC"}};

  if (BOOST_TEST_NE(c, nullptr)) {
    auto a = segment.get_segment_manager();

    for (const auto& input : inputs) {
      c->emplace(string_type(input.first.c_str(), a),
        string_type(input.second.c_str(), a));
    }

    c->rehash(c->bucket_count() + 1);

    if (BOOST_TEST_NE(it, nullptr)) {
      *it = c->begin();
    }
  }
}

template <class C>
typename std::enable_if<is_concurrent_container<C>::value, void>::type
parent(std::string const& shm_name_, char const* exe_name, C*)
{
  struct shm_remove
  {
    char const* shm_name;

    shm_remove(char const* shm_name_) : shm_name(shm_name_)
    {
      boost::interprocess::shared_memory_object::remove(shm_name);
    }
    ~shm_remove()
    {
      boost::interprocess::shared_memory_object::remove(shm_name);
    }
  } remover{shm_name_.c_str()};

  using container_type = C;

  std::size_t const shm_size = 64 * 1024;

  auto shm_name = remover.shm_name;
  boost::interprocess::managed_shared_memory segment(
    boost::interprocess::create_only, shm_name, shm_size);

  auto segment_mngr = segment.get_segment_manager();
  char_allocator char_alloc(segment_mngr);
  pair_allocator pair_alloc(segment_mngr);

  container_type* c =
    segment.construct<container_type>(shm_map_name)(pair_alloc);

  auto const old_bc = c->bucket_count();

  BOOST_TEST(c->empty());

  boost::process::child child(exe_name, shm_name);
  child.wait();
  int ret = child.exit_code();

  auto inputs =
    std::vector<std::pair<std::string, std::string> >{{"hello", "world"}};
  for (const auto& kvp : inputs) {
    auto const& key = kvp.first;
    auto const& value = kvp.second;
    c->emplace(string_type(key.data(), char_alloc),
      string_type(value.data(), char_alloc));
  }

  BOOST_TEST_EQ(ret, 0);
  BOOST_TEST_GT(c->size(), inputs.size());
  BOOST_TEST_GT(c->bucket_count(), old_bc);

  segment.destroy<container_type>(shm_map_name);
}

template <class C>
typename std::enable_if<is_concurrent_container<C>::value, void>::type
child(std::string const& shm_name, C*)
{
  using container_type = C;

  boost::interprocess::managed_shared_memory segment(
    boost::interprocess::open_only, shm_name.c_str());

  container_type* c = segment.find<container_type>(shm_map_name).first;

  BOOST_TEST(c->empty());

  std::vector<std::pair<std::string, std::string> > inputs = {
    {"aaa", "AAA"}, {"bbb", "BBB"}, {"ccc", "CCCC"}};

  if (BOOST_TEST_NE(c, nullptr)) {
    auto a = segment.get_segment_manager();

    for (const auto& input : inputs) {
      c->emplace(string_type(input.first.c_str(), a),
        string_type(input.second.c_str(), a));
    }

    c->rehash(c->bucket_count() + 1);
  }
}

std::string shm_name_sanitize(std::string const& exe_name)
{
  std::string s(exe_name);
  auto pos = std::remove_if(s.begin(), s.end(), [](char const c) {
    switch (c) {
    case '/':
    case '.':
    case '\\':
    case '-':
    case '_':
    case '+':
      return true;

    default:
      return false;
    }
  });
  s.erase(pos, s.end());
  s = "/" + s;
  return s.substr(0, 255);
}

void mmap_test(int argc, char const** argv)
{
  using container_type =
    decltype(get_container_type<BOOST_UNORDERED_FOA_MMAP_MAP_TYPE>());

  if (argc == 1) {
    auto uuid = to_string(boost::uuids::random_generator()());
    auto exe_name = argv[0];
    auto shm_name = shm_name_sanitize(std::string(exe_name) + uuid);
    container_type* p = nullptr;
    parent(shm_name, exe_name, p);
  } else {
    auto shm_name = std::string(argv[1]);
    container_type* p = nullptr;
    child(shm_name, p);
  }
}

int main(int argc, char const** argv)
{
  mmap_test(argc, argv);
  return boost::report_errors();
}
#endif
