////////////////////////////////////////
//
// (C) Copyright Ion Gaztanaga 2024. Distributed under the Boost
// Software License, Version 1.0. (See accompanying file
// LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// See http://www.boost.org/libs/container for documentation.
//
////////////////////////////////////////

#ifndef BOOST_CONTAINER_TEST_FLAT_MAP_TEST_HEADER
#define BOOST_CONTAINER_TEST_FLAT_MAP_TEST_HEADER

#include <boost/container/detail/config_begin.hpp>
#include "check_equal_containers.hpp"
#include "print_container.hpp"
#include <boost/move/iterator.hpp>
#include <boost/move/utility_core.hpp>
#include <boost/move/adl_move_swap.hpp>
#include <cstdlib>

#include <map>
#include "map_test.hpp"

namespace boost{
namespace container {
namespace test {

template< class RandomIt >
void random_shuffle( RandomIt first, RandomIt last )
{
   typedef typename boost::container::iterator_traits<RandomIt>::difference_type difference_type;
   difference_type n = last - first;
   for (difference_type i = n-1; i > 0; --i) {
      difference_type j = std::rand() % (i+1);
      if(j != i) {
         boost::adl_move_swap(first[i], first[j]);
      }
   }
}

template<class IntFlatMap, class IntFlatMultimap>
bool flat_tree_extract_adopt_test()
{
   using namespace boost::container;
   const std::size_t NumElements = 100;

   typedef typename IntFlatMap::sequence_type seq_t;

   //extract/adopt map
   {
      //Construction insertion
      IntFlatMap fmap, fmap_copy;

      for(std::size_t i = 0; i != NumElements; ++i){
         fmap.emplace(static_cast<int>(i), -static_cast<int>(i));
         fmap_copy.emplace(static_cast<int>(i), -static_cast<int>(i));
      }

      seq_t seq((BOOST_RV_REF(seq_t))fmap.extract_sequence());
      if(!fmap.empty())
         return false;
      if(!CheckEqualContainers(seq, fmap_copy))
         return false;

      #ifdef BOOST_CONTAINER_STD_PAIR_IS_MOVABLE
      for(std::size_t i = 0; i != NumElements; ++i){
         seq.emplace_back(static_cast<int>(i), -static_cast<int>(i));
      }

      boost::container::test::random_shuffle(seq.begin(), seq.end());
      #endif

      fmap.adopt_sequence(boost::move(seq));
      if(!CheckEqualContainers(fmap, fmap_copy))
         return false;
      if (!CheckEqualContainers(fmap.sequence(), fmap_copy.sequence()))
         return false;
   }

   //extract/adopt map, ordered_unique_range
   {
      //Construction insertion
      IntFlatMap fmap, fmap_copy;

      for(std::size_t i = 0; i != NumElements; ++i){
         fmap.emplace(static_cast<int>(i), -static_cast<int>(i));
         fmap_copy.emplace(static_cast<int>(i), -static_cast<int>(i));
      }

      seq_t seq((BOOST_RV_REF(seq_t))fmap.extract_sequence());

      if(!fmap.empty())
         return false;
      if(!CheckEqualContainers(seq, fmap_copy))
         return false;

      fmap.adopt_sequence(ordered_unique_range, boost::move(seq));
      if(!CheckEqualContainers(fmap, fmap_copy))
         return false;
      if (!CheckEqualContainers(fmap.sequence(), fmap_copy.sequence()))
         return false;
   }

   //extract/adopt multimap
   {
      //Construction insertion
      IntFlatMultimap fmmap, fmmap_copy;

      for(std::size_t i = 0; i != NumElements; ++i){
         fmmap.emplace(static_cast<int>(i), -static_cast<int>(i));
         fmmap_copy.emplace(static_cast<int>(i), -static_cast<int>(i));
         fmmap.emplace(static_cast<int>(i), -static_cast<int>(i));
         fmmap_copy.emplace(static_cast<int>(i), -static_cast<int>(i));
      }

      seq_t seq((BOOST_RV_REF(seq_t))fmmap.extract_sequence());
      if(!fmmap.empty())
         return false;
      if(!CheckEqualContainers(seq, fmmap_copy))
         return false;

      #ifdef BOOST_CONTAINER_STD_PAIR_IS_MOVABLE
      boost::container::test::random_shuffle(seq.begin(), seq.end());
      #endif

      fmmap.adopt_sequence(boost::move(seq));
      if(!CheckEqualContainers(fmmap, fmmap_copy))
         return false;
      if (!CheckEqualContainers(fmmap.sequence(), fmmap_copy.sequence()))
         return false;
   }

   //extract/adopt multimap, ordered_range
   {
      //Construction insertion
      IntFlatMultimap fmmap, fmmap_copy;

      for(std::size_t i = 0; i != NumElements; ++i){
         fmmap.emplace(static_cast<int>(i), -static_cast<int>(i));
         fmmap_copy.emplace(static_cast<int>(i), -static_cast<int>(i));
         fmmap.emplace(static_cast<int>(i), -static_cast<int>(i));
         fmmap_copy.emplace(static_cast<int>(i), -static_cast<int>(i));
      }

      seq_t seq((BOOST_RV_REF(seq_t))fmmap.extract_sequence());
      if(!fmmap.empty())
         return false;
      if(!CheckEqualContainers(seq, fmmap_copy))
         return false;

      fmmap.adopt_sequence(ordered_range, boost::move(seq));
      if(!CheckEqualContainers(fmmap, fmmap_copy))
         return false;
      if (!CheckEqualContainers(fmmap.sequence(), fmmap_copy.sequence()))
         return false;
   }
   return true;
}


template<class IntFlatMap, class IntFlatMultimap>
bool flat_tree_ordered_insertion_test()
{
   using namespace boost::container;
   const std::size_t NumElements = 100;

   //Ordered insertion multimap
   {
      std::multimap<int, int> int_mmap;
      for (std::size_t i = 0; i != NumElements; ++i) {
         int_mmap.insert(std::multimap<int, int>::value_type(static_cast<int>(i), static_cast<int>(i)));
      }
      //Construction insertion
      IntFlatMultimap fmmap(ordered_range, int_mmap.begin(), int_mmap.end());
      if (!CheckEqualContainers(int_mmap, fmmap))
         return false;
      //Insertion when empty
      fmmap.clear();
      fmmap.insert(ordered_range, int_mmap.begin(), int_mmap.end());
      if (!CheckEqualContainers(int_mmap, fmmap))
         return false;
      //Re-insertion
      fmmap.insert(ordered_range, int_mmap.begin(), int_mmap.end());
      std::multimap<int, int> int_mmap2(int_mmap);
      int_mmap2.insert(int_mmap.begin(), int_mmap.end());
      if (!CheckEqualContainers(int_mmap2, fmmap))
         return false;
      //Re-re-insertion
      fmmap.insert(ordered_range, int_mmap2.begin(), int_mmap2.end());
      std::multimap<int, int> int_mmap4(int_mmap2);
      int_mmap4.insert(int_mmap2.begin(), int_mmap2.end());
      if (!CheckEqualContainers(int_mmap4, fmmap))
         return false;
      //Re-re-insertion of even
      std::multimap<int, int> int_even_mmap;
      for (std::size_t i = 0; i < NumElements; i += 2) {
         int_mmap.insert(std::multimap<int, int>::value_type(static_cast<int>(i), static_cast<int>(i)));
      }
      fmmap.insert(ordered_range, int_even_mmap.begin(), int_even_mmap.end());
      int_mmap4.insert(int_even_mmap.begin(), int_even_mmap.end());
      if (!CheckEqualContainers(int_mmap4, fmmap))
         return false;
   }

   //Ordered insertion map
   {
      std::map<int, int> int_map;
      for (std::size_t i = 0; i != NumElements; ++i) {
         int_map.insert(std::map<int, int>::value_type(static_cast<int>(i), static_cast<int>(i)));
      }
      //Construction insertion
      IntFlatMap fmap(ordered_unique_range, int_map.begin(), int_map.end());
      if (!CheckEqualContainers(int_map, fmap))
         return false;
      //Insertion when empty
      fmap.clear();
      fmap.insert(ordered_unique_range, int_map.begin(), int_map.end());
      if (!CheckEqualContainers(int_map, fmap))
         return false;
      //Re-insertion
      fmap.insert(ordered_unique_range, int_map.begin(), int_map.end());
      std::map<int, int> int_map2(int_map);
      int_map2.insert(int_map.begin(), int_map.end());
      if (!CheckEqualContainers(int_map2, fmap))
         return false;
      //Re-re-insertion
      fmap.insert(ordered_unique_range, int_map2.begin(), int_map2.end());
      std::map<int, int> int_map4(int_map2);
      int_map4.insert(int_map2.begin(), int_map2.end());
      if (!CheckEqualContainers(int_map4, fmap))
         return false;
      //Re-re-insertion of even
      std::map<int, int> int_even_map;
      for (std::size_t i = 0; i < NumElements; i += 2) {
         int_map.insert(std::map<int, int>::value_type(static_cast<int>(i), static_cast<int>(i)));
      }
      fmap.insert(ordered_unique_range, int_even_map.begin(), int_even_map.end());
      int_map4.insert(int_even_map.begin(), int_even_map.end());
      if (!CheckEqualContainers(int_map4, fmap))
         return false;
   }

   return true;
}

template<class MyBoostFlatMap
        ,class MyStdMap
        ,class MyBoostFlatMultimap
        ,class MyStdMultiMap>
int flat_map_test()
{
   if(0 != map_test<MyBoostFlatMap, MyStdMap, MyBoostFlatMultimap, MyStdMultiMap>())
      return 1;

   if (!flat_tree_ordered_insertion_test<MyBoostFlatMap, MyBoostFlatMultimap>())
      return 1;

   if (!flat_tree_extract_adopt_test<MyBoostFlatMap, MyBoostFlatMultimap>())
      return 1;
   return 0;
}

}  //namespace test{
}  //namespace container {
}  //namespace boost{

#include <boost/container/detail/config_end.hpp>

#endif   //#ifndef BOOST_CONTAINER_TEST_FLAT_MAP_TEST_HEADER
