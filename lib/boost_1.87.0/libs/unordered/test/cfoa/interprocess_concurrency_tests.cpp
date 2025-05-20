// Copyright (C) 2024 Joaquin M Lopez Munoz
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/unordered/concurrent_flat_map.hpp>
#include <atomic>
#include <boost/core/lightweight_test.hpp>
#include <boost/interprocess/allocators/allocator.hpp>
#include <boost/interprocess/managed_shared_memory.hpp>
#include <boost/process/child.hpp>
#include <boost/uuid/random_generator.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <chrono>
#include <cstdlib>
#include <random>
#include <string>
#include <thread>
#include <vector>

namespace bip=boost::interprocess;

static char const* container_name = "shared_map";
static char const* start_name = "shared_start";
static constexpr int NUM_OPS_PER_CHILD = 10000;

using value_type = std::pair<const int, int>;
using allocator = bip::allocator<
  value_type, bip::managed_shared_memory::segment_manager>;
using container = boost::concurrent_flat_map<
  int, int, boost::hash<int>, std::equal_to<int>, allocator>;

int parent(const char* exe_name)
{
  static constexpr int NUM_CHILDS = 10;
  static constexpr std::size_t SEGMENT_SIZE = 64*1024;

  std::string segment_name_str = to_string(boost::uuids::random_generator()());
  auto* segment_name = segment_name_str.c_str();

  struct segment_remover
  {
    char const* name;

    segment_remover(char const* name_) : name(name_)
    {
      bip::shared_memory_object::remove(name);
    }
    ~segment_remover()
    {
      bip::shared_memory_object::remove(name);
    }
  } remover(segment_name);

  bip::managed_shared_memory segment(
    bip::create_only, segment_name, SEGMENT_SIZE);
  container& c = *segment.construct<container>(container_name)(
    allocator(segment.get_segment_manager()));
  std::atomic_int& start = *segment.construct<std::atomic_int>(start_name)(0);

  std::vector<boost::process::child> children;
  for (int i = 0; i < NUM_CHILDS; ++i) {
    children.emplace_back(exe_name, std::to_string(i), segment_name);
  }

  start.store(1);

  for (auto& child : children) {
    child.wait();
    BOOST_TEST_EQ(child.exit_code(), 0);
  }

  int num_ops = 0;
  c.cvisit_all([&](const value_type& x) {
    num_ops += x.second;
  });
  BOOST_TEST_EQ(num_ops, NUM_CHILDS * NUM_OPS_PER_CHILD);

  return boost::report_errors();
}

int child(int id,const char* segment_name)
{
  bip::managed_shared_memory segment(bip::open_only, segment_name);
  container& c = *segment.find<container>(container_name).first;
  std::atomic_int& start = *segment.find<std::atomic_int>(start_name).first;

  while(!start.load()){
    std::this_thread::sleep_for(std::chrono::milliseconds(1));
  }

  std::mt19937 rnd((unsigned int) id);
  std::geometric_distribution<int> d(0.15);

  for(int i = 0; i < NUM_OPS_PER_CHILD; ++i) {
    c.emplace_or_visit(d(rnd), 1, [&](value_type& x) {
      ++x.second;

      // artificially increase contention
      volatile unsigned int n = 10000;
      while(n--) ;
    });
  }
  return 0;
}

int main(int argc, char** argv)
{
  if (argc == 1) {
    return parent(argv[0]);
  }
  else { 
    return child(std::atoi(argv[1]),argv[2]);
  }
}
