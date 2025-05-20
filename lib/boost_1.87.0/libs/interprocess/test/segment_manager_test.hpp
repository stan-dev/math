//////////////////////////////////////////////////////////////////////////////
//
// (C) Copyright Ion Gaztanaga 2004-2019. Distributed under the Boost
// Software License, Version 1.0. (See accompanying file
// LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// See http://www.boost.org/libs/interprocess for documentation.
//
//////////////////////////////////////////////////////////////////////////////

#include <boost/interprocess/indexes/flat_map_index.hpp>
#include <boost/interprocess/indexes/map_index.hpp>
#include <boost/interprocess/indexes/null_index.hpp>
#include <boost/interprocess/indexes/iset_index.hpp>
#include <boost/interprocess/indexes/iunordered_set_index.hpp>

#include <boost/interprocess/mem_algo/simple_seq_fit.hpp>
#include <boost/interprocess/mem_algo/rbtree_best_fit.hpp>
#include <boost/interprocess/mapped_region.hpp>
#include <boost/interprocess/segment_manager.hpp>
#include <boost/interprocess/shared_memory_object.hpp>
#include <boost/interprocess/sync/mutex_family.hpp>
#include <boost/interprocess/exceptions.hpp>
#include <boost/interprocess/detail/utilities.hpp>
#include "get_process_id_name.hpp"
#include <cstddef>
#include <new>
#include <cstring>
#include <typeinfo>
#include <iostream>

using namespace boost::interprocess;

template<std::size_t Align>
struct IntLike;

//Old GCC versions emit incorrect warnings like
//"requested alignment 256 is larger than 128"
#if defined(BOOST_GCC) && (BOOST_GCC <= 100000)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
#endif

#define BOOST_INTERPROCESS_ALIGNED_INTLIKE(A)\
template<>\
struct IntLike<A>\
{\
   IntLike(){}\
\
   IntLike(int i) : data(i) {}\
\
   BOOST_ALIGNMENT(A) int data;\
\
   operator int() const {  return data;   }\
};\
//

//Up to 4K alignment (typical page size)
BOOST_INTERPROCESS_ALIGNED_INTLIKE(16)
BOOST_INTERPROCESS_ALIGNED_INTLIKE(32)
BOOST_INTERPROCESS_ALIGNED_INTLIKE(64)
BOOST_INTERPROCESS_ALIGNED_INTLIKE(128)
BOOST_INTERPROCESS_ALIGNED_INTLIKE(256)

#undef BOOST_INTERPROCESS_ALIGNED_INTLIKE

#if defined(BOOST_GCC) && (BOOST_GCC <= 100000)
#pragma GCC diagnostic pop
#endif

template <class SegmentManager>
struct atomic_func_test
{
   SegmentManager &rsm;
   int *object;

   atomic_func_test(SegmentManager &sm)
      : rsm(sm), object()
   {}

   void operator()()
   {
      object = rsm.template find<int>("atomic_func_find_object").first;
   }
   private:
   atomic_func_test operator=(const atomic_func_test&);
   atomic_func_test(const atomic_func_test&);
};


template <class SegmentManager>
bool test_allocate_deallocate(SegmentManager* seg_mgr, mapped_region& mapping)
{
   typedef typename SegmentManager::size_type size_type;
   const std::size_t MappedRegionSize = mapping.get_size();

   for (std::size_t size = 1; size <= (MappedRegionSize / 2); size <<= 1 ) {
      const std::size_t free_mem_before = seg_mgr->get_free_memory();

      //Allocate memory
      void* mem = seg_mgr->allocate(size + 1);
      const size_type free_mem = seg_mgr->get_free_memory();
      if (free_mem >= (free_mem_before-size))
         return false;
      if (seg_mgr->all_memory_deallocated())
         return false;
      //Allocate half of the rest
      const size_type Size2 = free_mem / 2;
      void* mem2 = seg_mgr->allocate(size_type(Size2 + 1), std::nothrow);

      //Sanity checks
      if (seg_mgr->get_free_memory() >= Size2)
         return false;
      if (seg_mgr->size(mem) < (size + 1))
         return false;
      if (seg_mgr->size(mem2) < (Size2 + 1))
         return false;

      //Deallocate both
      seg_mgr->deallocate(mem);
      seg_mgr->deallocate(mem2);

      //Sanity checks again
      if (!seg_mgr->all_memory_deallocated())
         return false;
      if (seg_mgr->get_free_memory() != free_mem_before)
         return false;

      //Try an imposible size to test error is signalled
      bool operation_throws = false;
      BOOST_INTERPROCESS_TRY{ seg_mgr->allocate(MappedRegionSize * 2); }
      BOOST_INTERPROCESS_CATCH(interprocess_exception&) { operation_throws = true; }
      BOOST_INTERPROCESS_CATCH_END
      if (!operation_throws)
         return false;

      if (seg_mgr->get_free_memory() != free_mem_before)
         return false;

      if (seg_mgr->allocate(MappedRegionSize*2u, std::nothrow))
         return false;

      if (seg_mgr->get_free_memory() != free_mem_before)
         return false;
   }
   return true;
}

template <class SegmentManager>
bool test_allocate_aligned(SegmentManager* seg_mgr, mapped_region& mapping)
{
   const std::size_t MappedRegionSize = mapping.get_size();
   const std::size_t free_mem_before = seg_mgr->get_free_memory();
   const std::size_t InitialAlignment = SegmentManager::memory_algorithm::Alignment;
   const std::size_t RegionAlignment = mapped_region::get_page_size();

   for( std::size_t alignment = InitialAlignment
      ; (alignment <= MappedRegionSize/8 && alignment <= RegionAlignment/4)
      ; alignment <<= 1u) {

      //Allocate two buffers and test the alignment inside the mapped region
      void *mem = seg_mgr->allocate_aligned(MappedRegionSize/8, alignment);
      if(seg_mgr->all_memory_deallocated())
         return false;

      std::size_t offset = static_cast<std::size_t>
         (static_cast<const char *>(mem) -  static_cast<const char *>(mapping.get_address()));
      if(offset & (alignment -1))
         return false;
      void *mem2 = seg_mgr->allocate_aligned(MappedRegionSize/4, alignment, std::nothrow);
      std::size_t offset2 = static_cast<std::size_t>
         (static_cast<const char *>(mem2) -  static_cast<const char *>(mapping.get_address()));
      if(offset2 & (alignment -1))
         return false;

      //Deallocate them
      seg_mgr->deallocate(mem);
      seg_mgr->deallocate(mem2);
      if(!seg_mgr->all_memory_deallocated())
         return false;
      if(seg_mgr->get_free_memory() != free_mem_before)
         return false;

      //Try an imposible size to test error is signalled
      bool operation_throws = false;
      BOOST_INTERPROCESS_TRY{  seg_mgr->allocate_aligned(MappedRegionSize*2, alignment);  }
      BOOST_INTERPROCESS_CATCH(interprocess_exception&){ operation_throws = true; }
      BOOST_INTERPROCESS_CATCH_END
      if (!operation_throws)
         return false;

      if (seg_mgr->allocate_aligned(MappedRegionSize*2, alignment, std::nothrow))
         return false;
      if(seg_mgr->get_free_memory() != free_mem_before)
         return false;
      if(seg_mgr->allocate_aligned(MappedRegionSize*2, alignment, std::nothrow))
         return false;
      if(seg_mgr->get_free_memory() != free_mem_before)
         return false;
   }
   return true;
}

template <class SegmentManager>
bool test_shrink_to_fit(SegmentManager* seg_mgr, mapped_region &)
{
   typedef typename SegmentManager::size_type size_type;
   const std::size_t free_mem_before = seg_mgr->get_free_memory();
   std::size_t size_before = seg_mgr->get_size();
   seg_mgr->shrink_to_fit();
   if (!seg_mgr->all_memory_deallocated())
      return false;
   std::size_t empty_shrunk_size = seg_mgr->get_size();
   std::size_t empty_shrunk_free_mem = seg_mgr->get_free_memory();
   if (empty_shrunk_size >= size_before)
      return false;
   if (empty_shrunk_free_mem >= size_before)
      return false;
   seg_mgr->grow(size_type(size_before - empty_shrunk_size));
   if (seg_mgr->get_size() != size_before)
      return false;
   if (seg_mgr->get_free_memory() != free_mem_before)
      return false;
   if (!seg_mgr->all_memory_deallocated())
      return false;

   return true;
}

template <class SegmentManager>
bool test_zero_free_memory(SegmentManager* seg_mgr, mapped_region &mapping)
{
   typedef typename SegmentManager::size_type size_type;
   const std::size_t MappedRegionSize = mapping.get_size();
   const std::size_t free_mem_before = seg_mgr->get_free_memory();
   const size_type Size(MappedRegionSize / 2 + 1), Size2(MappedRegionSize / 8);
   void* mem = seg_mgr->allocate(Size);
   void* mem2 = seg_mgr->allocate(Size2);
   //Mark memory to non-zero
   std::memset(mem, 0xFF, Size);
   std::memset(mem2, 0xFF, Size2);
   //Deallocate and check still non-zero
   seg_mgr->deallocate(mem);
   seg_mgr->deallocate(mem2);
   {  //Use byte per byte comparison as "static unsigned char zerobuf[Size]"
      //seems to be problematic in some compilers
      unsigned char* const mem_uch_ptr = static_cast<unsigned char*>(mem);
      unsigned char* const mem2_uch_ptr = static_cast<unsigned char*>(mem2);
      size_type zeroes = 0;
      for (size_type i = 0; i != Size; ++i) {
         if (!mem_uch_ptr[i])
            ++zeroes;
      }
      if (zeroes == Size)
         return false;

      zeroes = 0;
      for (size_type i = 0; i != Size2; ++i) {
         if (!mem2_uch_ptr[i])
            ++zeroes;
      }
      if (zeroes == Size2)
         return false;
   }
   //zero_free_memory and check it's zeroed
   seg_mgr->zero_free_memory();
   //TODO: some parts are not zeroed because they are used
   //as internal metadata, find a way to test this
   //if(std::memcmp(mem,  zerobuf, Size))
      //return false;
   //if(std::memcmp(mem2, zerobuf, Size2))
      //return false;
   if (seg_mgr->get_free_memory() != free_mem_before)
      return false;
   if (!seg_mgr->all_memory_deallocated())
      return false;
   return true;
}


template <class IntLike, class SegmentManager>
bool test_anonymous_object_type(SegmentManager* seg_mgr, mapped_region& mapping)
{
   const std::size_t MappedRegionSize = mapping.get_size();
   const std::size_t free_mem_before = seg_mgr->get_free_memory();

   //Construct single object
   {
      IntLike* int_object = seg_mgr->template construct<IntLike>(anonymous_instance)();
      BOOST_ASSERT(is_ptr_aligned(int_object, boost::move_detail::alignment_of<IntLike>::value));
      if (!is_ptr_aligned(int_object, boost::move_detail::alignment_of<IntLike>::value))
         return false;
      if (1 != seg_mgr->get_instance_length(int_object))
         return false;
      if (anonymous_type != seg_mgr->get_instance_type(int_object))
         return false;
      if (seg_mgr->get_instance_name(int_object))
         return false;
      seg_mgr->destroy_ptr(int_object);
   }
   {
      //Construct array object
      IntLike* int_array = seg_mgr->template construct_it<IntLike>(anonymous_instance, std::nothrow)[5]();
      BOOST_ASSERT(is_ptr_aligned(int_array, boost::move_detail::alignment_of<IntLike>::value));
      if (!is_ptr_aligned(int_array, boost::move_detail::alignment_of<IntLike>::value))
         return false;
      if (5 != seg_mgr->get_instance_length(int_array))
         return false;
      if (anonymous_type != seg_mgr->get_instance_type(int_array))
         return false;
      if (seg_mgr->get_instance_name(int_array))
         return false;
      seg_mgr->destroy_ptr(int_array);
   }
   {
      //Construct array object from it
      const signed char int_array_values[10] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
      IntLike* int_array = seg_mgr->template construct_it<IntLike>(anonymous_instance, std::nothrow)[10](&int_array_values[0]);
      BOOST_ASSERT(is_ptr_aligned(int_array, boost::move_detail::alignment_of<IntLike>::value));
      if (!is_ptr_aligned(int_array, boost::move_detail::alignment_of<IntLike>::value))
         return false;
      if (10 != seg_mgr->get_instance_length(int_array))
         return false;
      if (anonymous_type != seg_mgr->get_instance_type(int_array))
         return false;
      if (seg_mgr->get_instance_name(int_array))
         return false;
      seg_mgr->destroy_ptr(int_array);
   }

   //Try an imposible size to test error is signalled
   {
      bool operation_throws = false;
      BOOST_INTERPROCESS_TRY{ seg_mgr->template construct<IntLike>(anonymous_instance)[MappedRegionSize](); }
      BOOST_INTERPROCESS_CATCH(interprocess_exception&) { operation_throws = true; }
      BOOST_INTERPROCESS_CATCH_END
      if (!operation_throws)
         return false;
      if (seg_mgr->get_free_memory() != free_mem_before)
         return false;
   }
   {
      if (seg_mgr->template construct<IntLike>(anonymous_instance, std::nothrow)[MappedRegionSize]())
      if (seg_mgr->get_free_memory() != free_mem_before)
         return false;
   }
   {
      const signed char int_array_values[10] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
      if (seg_mgr->template construct_it<IntLike>(anonymous_instance, std::nothrow)[MappedRegionSize](&int_array_values[0]))
         return false;
      if (seg_mgr->get_free_memory() != free_mem_before)
         return false;
   }

   if (!seg_mgr->all_memory_deallocated())
      return false;
   return true;
}

template <class SegmentManager>
bool test_anonymous_object(SegmentManager* seg_mgr, mapped_region& mapping)
{
   if (!test_anonymous_object_type<signed char>(seg_mgr, mapping))
      return false;

   if (!test_anonymous_object_type<short int>(seg_mgr, mapping))
      return false;

   if (!test_anonymous_object_type<int>(seg_mgr, mapping))
      return false;

   if (!test_anonymous_object_type<long int>(seg_mgr, mapping))
      return false;
   #if (BOOST_INTERPROCESS_SEGMENT_MANAGER_ABI >= 2)
   if (!test_anonymous_object_type<long long int>(seg_mgr, mapping))
      return false;

   if (!test_anonymous_object_type<IntLike<16> >(seg_mgr, mapping))
      return false;

   if (!test_anonymous_object_type<IntLike<32> >(seg_mgr, mapping))
      return false;

   if (!test_anonymous_object_type<IntLike<64> >(seg_mgr, mapping))
      return false;

   if (!test_anonymous_object_type<IntLike<128> >(seg_mgr, mapping))
      return false;
   #endif   //#if (BOOST_INTERPROCESS_SEGMENT_MANAGER_ABI >= 2)
   return true;
}

template <class IntLike, class SegmentManager>
bool test_named_object_type(SegmentManager* seg_mgr, mapped_region& mapping)
{
   const std::size_t MappedRegionSize = mapping.get_size();
   const std::size_t free_mem_before = seg_mgr->get_free_memory();

   const char* const object1_name = "object1";
   const char* const object2_name = "object2";
   const signed char int_array_values[10] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };

   for (std::size_t i = 0; i != 4; ++i) {
      if (seg_mgr->template find<IntLike>(object1_name).first)
         return false;
      //Single element construction
      IntLike* int_object = 0;
      switch (i) {
      case 0:
         int_object = seg_mgr->template construct<IntLike>(object1_name)();
         break;
      case 1:
         int_object = seg_mgr->template construct<IntLike>(object1_name, std::nothrow)();
         break;
      case 2:
         int_object = seg_mgr->template find_or_construct<IntLike>(object1_name)();
         break;
      case 3:
         int_object = seg_mgr->template find_or_construct<IntLike>(object1_name, std::nothrow)();
         break;
      }


      if (!is_ptr_aligned(int_object, boost::move_detail::alignment_of<IntLike>::value)){
         std::cout << "\ntype/alignment: " << typeid(IntLike).name() << "/" << boost::move_detail::alignment_of<IntLike>::value << "\n segment_manager: " << typeid(SegmentManager).name()
                   << "\nmem alignment: " << SegmentManager::memory_algorithm::Alignment <<std::endl;
         return false;
      }

      std::pair<IntLike*, std::size_t> find_ret = seg_mgr->template find<IntLike>(object1_name);
      if (int_object != find_ret.first)
         return false;
      if (1 != find_ret.second)
         return false;
      if (1 != seg_mgr->get_instance_length(int_object))
         return false;
      if (named_type != seg_mgr->get_instance_type(int_object))
         return false;
      if (std::strcmp(object1_name, seg_mgr->get_instance_name(int_object)))
         return false;

      //Array construction
      if (seg_mgr->template find<IntLike>(object2_name).first)
         return false;
      IntLike* int_array = 0;
      switch (i) {
      case 0:
         int_array = seg_mgr->template construct_it<IntLike>(object2_name)[10](&int_array_values[0]);
         break;
      case 1:
         int_array = seg_mgr->template construct_it<IntLike>(object2_name, std::nothrow)[10](&int_array_values[0]);
         break;
      case 2:
         int_array = seg_mgr->template find_or_construct_it<IntLike>(object2_name)[10](&int_array_values[0]);
         break;
      case 3:
         int_array = seg_mgr->template find_or_construct_it<IntLike>(object2_name, std::nothrow)[10](&int_array_values[0]);
         break;
      }

      BOOST_ASSERT(is_ptr_aligned(int_array, boost::move_detail::alignment_of<IntLike>::value));
      if (!is_ptr_aligned(int_array, boost::move_detail::alignment_of<IntLike>::value))
         return false;

      std::pair<IntLike*, std::size_t> find_ret2 = seg_mgr->template find<IntLike>(object2_name);
      if (int_array != find_ret2.first)
         return false;
      if (10 != find_ret2.second)
         return false;
      if (10 != seg_mgr->get_instance_length(int_array))
         return false;
      if (named_type != seg_mgr->get_instance_type(int_array))
         return false;
      if (std::strcmp(object2_name, seg_mgr->get_instance_name(int_array)))
         return false;
      if (seg_mgr->get_num_named_objects() != 2)
         return false;
      typename SegmentManager::const_named_iterator nb(seg_mgr->named_begin());
      typename SegmentManager::const_named_iterator ne(seg_mgr->named_end());
      for (std::size_t j = 0, imax = seg_mgr->get_num_named_objects(); j != imax; ++j) { ++nb; }
      if (nb != ne)
         return false;
      seg_mgr->destroy_ptr(int_object);
      seg_mgr->template destroy<IntLike>(object2_name);
   }

   //Try an imposible size to test error is signalled
   {
      bool operation_throws = false;
      BOOST_INTERPROCESS_TRY{ seg_mgr->template construct<IntLike>(object1_name)[MappedRegionSize](); }
      BOOST_INTERPROCESS_CATCH(interprocess_exception&) { operation_throws = true;}
      BOOST_INTERPROCESS_CATCH_END
      if (!operation_throws)
         return false;

      if (seg_mgr->template construct<IntLike>(object2_name, std::nothrow)[MappedRegionSize]())
         return false;
   }
   {
      bool operation_throws = false;
      BOOST_INTERPROCESS_TRY{ seg_mgr->template construct_it<IntLike>(object1_name)[MappedRegionSize](&int_array_values[0]); }
      BOOST_INTERPROCESS_CATCH(interprocess_exception&) { operation_throws = true; }
      BOOST_INTERPROCESS_CATCH_END
      if (!operation_throws)
         return false;

      if (seg_mgr->template construct_it<IntLike>(object2_name, std::nothrow)[MappedRegionSize](&int_array_values[0]))
         return false;
   }

   seg_mgr->shrink_to_fit_indexes();
   if (seg_mgr->get_free_memory() != free_mem_before)
      return false;
   if (!seg_mgr->all_memory_deallocated())
      return false;
   seg_mgr->reserve_named_objects(1);

   //In indexes with no capacity() memory won't be allocated so don't check anything was allocated.
   //if(seg_mgr->all_memory_deallocated())  return false;
   seg_mgr->shrink_to_fit_indexes();
   if (!seg_mgr->all_memory_deallocated())
      return false;
   return true;
}

template <class SegmentManager>
bool test_named_object(SegmentManager* seg_mgr, mapped_region& mapping)
{
   if (!test_named_object_type<signed char>(seg_mgr, mapping))
      return false;

   if (!test_named_object_type<short int>(seg_mgr, mapping))
      return false;

   if (!test_named_object_type<int>(seg_mgr, mapping))
      return false;

   if (!test_named_object_type<long int>(seg_mgr, mapping))
      return false;
   #if (BOOST_INTERPROCESS_SEGMENT_MANAGER_ABI >= 2)
   if (!test_named_object_type<long long int>(seg_mgr, mapping))
      return false;

   if (!test_named_object_type<IntLike<16> >(seg_mgr, mapping))
      return false;

   if (!test_named_object_type<IntLike<32> >(seg_mgr, mapping))
      return false;

   if (!test_named_object_type<IntLike<64> >(seg_mgr, mapping))
      return false;

   if (!test_named_object_type<IntLike<128> >(seg_mgr, mapping))
      return false;
   #endif   //#if (BOOST_INTERPROCESS_SEGMENT_MANAGER_ABI >= 2)
   return true;
}

template <class IntLike, class SegmentManager>
bool test_unique_object_type(SegmentManager* seg_mgr, mapped_region& mapping)
{
   const std::size_t MappedRegionSize = mapping.get_size();
   const std::size_t free_mem_before = seg_mgr->get_free_memory();

   const signed char int_array_values[10] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };

   for (std::size_t i = 0; i != 4; ++i) {
      if (seg_mgr->template find<IntLike>(unique_instance).first)
         return false;
      //Single element construction
      IntLike* int_object = 0;
      switch (i) {
      case 0:
         int_object = seg_mgr->template construct<IntLike>(unique_instance)();
         break;
      case 1:
         int_object = seg_mgr->template construct<IntLike>(unique_instance, std::nothrow)();
         break;
      case 2:
         int_object = seg_mgr->template find_or_construct<IntLike>(unique_instance)();
         break;
      case 3:
         int_object = seg_mgr->template find_or_construct<IntLike>(unique_instance, std::nothrow)();
         break;
      }

      BOOST_ASSERT(is_ptr_aligned(int_object, boost::move_detail::alignment_of<IntLike>::value));
      if (!is_ptr_aligned(int_object, boost::move_detail::alignment_of<IntLike>::value))
         return false;

      std::pair<IntLike*, std::size_t> find_ret = seg_mgr->template find<IntLike>(unique_instance);
      if (int_object != find_ret.first)
         return false;
      if (1 != find_ret.second)
         return false;
      if (1 != seg_mgr->get_instance_length(int_object))
         return false;
      if (unique_type != seg_mgr->get_instance_type(int_object))
         return false;
      if (std::strcmp(typeid(IntLike).name(), seg_mgr->get_instance_name(int_object)))
         return false;
      //Array construction
      if (!seg_mgr->template find<IntLike>(unique_instance).first)
         return false;

      seg_mgr->destroy_ptr(int_object);

      IntLike* int_array = 0;
      switch (i) {
      case 0:
         int_array = seg_mgr->template construct_it<IntLike>(unique_instance)[10](&int_array_values[0]);
         break;
      case 1:
         int_array = seg_mgr->template construct_it<IntLike>(unique_instance, std::nothrow)[10](&int_array_values[0]);
         break;
      case 2:
         int_array = seg_mgr->template find_or_construct_it<IntLike>(unique_instance)[10](&int_array_values[0]);
         break;
      case 3:
         int_array = seg_mgr->template find_or_construct_it<IntLike>(unique_instance, std::nothrow)[10](&int_array_values[0]);
         break;
      }

      BOOST_ASSERT(is_ptr_aligned(int_array, boost::move_detail::alignment_of<IntLike>::value));
      if (!is_ptr_aligned(int_array, boost::move_detail::alignment_of<IntLike>::value))
         return false;

      std::pair<IntLike*, std::size_t> find_ret2 = seg_mgr->template find<IntLike>(unique_instance);
      if (int_array != find_ret2.first)
         return false;
      if (10 != find_ret2.second)
         return false;
      if (10 != seg_mgr->get_instance_length(int_array))
         return false;
      if (unique_type != seg_mgr->get_instance_type(int_array))
         return false;
      if (std::strcmp(typeid(IntLike).name(), seg_mgr->get_instance_name(int_array)))
         return false;
      if (seg_mgr->get_num_unique_objects() != 1)
         return false;
      typename SegmentManager::const_unique_iterator nb(seg_mgr->unique_begin());
      typename SegmentManager::const_unique_iterator ne(seg_mgr->unique_end());
      for (std::size_t j = 0, imax = seg_mgr->get_num_unique_objects(); j != imax; ++j) { ++nb; }
      if (nb != ne)
         return false;
      seg_mgr->template destroy<IntLike>(unique_instance);
   }
   //Try an imposible size to test error is signalled
   {
      bool operation_throws = false;
      BOOST_INTERPROCESS_TRY{ seg_mgr->template construct<IntLike>(unique_instance)[MappedRegionSize](); }
      BOOST_INTERPROCESS_CATCH(interprocess_exception&) { operation_throws = true; }
      BOOST_INTERPROCESS_CATCH_END
      if (!operation_throws)
         return false;
      if (seg_mgr->template construct<IntLike>(unique_instance, std::nothrow)[MappedRegionSize]())
         return false;
   }
   {
      bool operation_throws = false;
      BOOST_INTERPROCESS_TRY{ seg_mgr->template construct_it<IntLike>(unique_instance)[MappedRegionSize](&int_array_values[0]); }
      BOOST_INTERPROCESS_CATCH(interprocess_exception&) { operation_throws = true; }
      BOOST_INTERPROCESS_CATCH_END
      if (!operation_throws)
         return false;
      if (seg_mgr->template construct_it<IntLike>(unique_instance, std::nothrow)[MappedRegionSize](&int_array_values[0]))
         return false;
   }

   seg_mgr->shrink_to_fit_indexes();

   if (seg_mgr->get_free_memory() != free_mem_before)
      return false;
   if (!seg_mgr->all_memory_deallocated())
      return false;

   seg_mgr->reserve_unique_objects(1);

   //In indexes with no capacity() memory won't be allocated so don't check anything was allocated.
   //if(seg_mgr->all_memory_deallocated())  return false;
   seg_mgr->shrink_to_fit_indexes();
   if (!seg_mgr->all_memory_deallocated())
      return false;
   return true;
}

template <class SegmentManager>
bool test_unique_object(SegmentManager* seg_mgr, mapped_region& mapping)
{
   if (!test_unique_object_type<signed char>(seg_mgr, mapping))
      return false;

   if (!test_unique_object_type<short int>(seg_mgr, mapping))
      return false;

   if (!test_unique_object_type<int>(seg_mgr, mapping))
      return false;

   if (!test_unique_object_type<long int>(seg_mgr, mapping))
      return false;
   #if (BOOST_INTERPROCESS_SEGMENT_MANAGER_ABI >= 2)
   if (!test_unique_object_type<long long int>(seg_mgr, mapping))
      return false;

   if (!test_unique_object_type<IntLike<16> >(seg_mgr, mapping))
      return false;

   if (!test_unique_object_type<IntLike<32> >(seg_mgr, mapping))
      return false;

   if (!test_unique_object_type<IntLike<64> >(seg_mgr, mapping))
      return false;

   if (!test_unique_object_type<IntLike<128> >(seg_mgr, mapping))
      return false;
   #endif   //#if (BOOST_INTERPROCESS_SEGMENT_MANAGER_ABI >= 2)
   return true;
}

template <class SegmentManager>
bool test_atomic_func(SegmentManager* seg_mgr, mapped_region& )
{
   if (!seg_mgr->all_memory_deallocated())
      return false;
   int* int_object = seg_mgr->template construct<int>("atomic_func_find_object")();
   atomic_func_test<SegmentManager> func(*seg_mgr);
   seg_mgr->atomic_func(func);
   if (int_object != func.object)
      return 1;
   seg_mgr->destroy_ptr(int_object);
   seg_mgr->shrink_to_fit_indexes();
   if (!seg_mgr->all_memory_deallocated())
      return false;
   return true;
}

template <class SegmentManager>
bool test_allocator_deleter(SegmentManager* seg_mgr, mapped_region&)
{//test allocator/deleter
   if (!seg_mgr->all_memory_deallocated())
      return false;
   typedef typename SegmentManager::template allocator<float>::type allocator_t;

   allocator_t alloc(seg_mgr->template get_allocator<float>());

   if (!seg_mgr->all_memory_deallocated())
      return false;
   offset_ptr<float> f = alloc.allocate(50);
   if (seg_mgr->all_memory_deallocated())
      return false;
   alloc.deallocate(f, 50);
   if (!seg_mgr->all_memory_deallocated())
      return false;
   typedef typename SegmentManager::template deleter<float>::type deleter_t;
   deleter_t delet(seg_mgr->template get_deleter<float>());
   delet(seg_mgr->template construct<float>(anonymous_instance)());
   if (!seg_mgr->all_memory_deallocated())
      return false;
   return true;
}

template <class SegmentManager>
bool test_get_memory_algorithm(SegmentManager* seg_mgr, mapped_region&)
{
   {
      typename SegmentManager::memory_algorithm& mem_algo =
         seg_mgr->get_memory_algorithm();
      const typename SegmentManager::memory_algorithm& const_mem_algo =
         const_cast<const SegmentManager*>(seg_mgr)->get_memory_algorithm();
      if (&mem_algo != &const_mem_algo)
         return false;
   }
   return true;
}


template <class SegmentManager>
bool test_segment_manager()
{
   const unsigned int MappedRegionSize = 1024*64u;
   std::string shmname(test::get_process_id_name());

   shared_memory_object::remove(shmname.c_str());
   shared_memory_object sh_mem( create_only, shmname.c_str(), read_write );
   sh_mem.truncate( MappedRegionSize );
   mapped_region mapping( sh_mem, read_write );

   //Remove shared memory to minimize risk of garbage on crash
   shared_memory_object::remove(shmname.c_str());

   SegmentManager* seg_mgr = new( mapping.get_address() ) SegmentManager( MappedRegionSize );
   std::size_t size_before = seg_mgr->get_size();

   if(size_before != MappedRegionSize)
      return false;
   if(!seg_mgr->all_memory_deallocated())
      return false;
   if(seg_mgr->get_min_size() >= MappedRegionSize)
      return false;

   if (!test_allocate_deallocate(seg_mgr, mapping))
      return false;

   if (!test_allocate_aligned(seg_mgr, mapping))
      return false;

   if (!test_shrink_to_fit(seg_mgr, mapping))
      return false;

   if (!test_zero_free_memory(seg_mgr, mapping))
      return false;

   if (!test_anonymous_object(seg_mgr, mapping))
      return false;

   if (!test_named_object(seg_mgr, mapping))
      return false;

   if (!test_unique_object(seg_mgr, mapping))
      return false;

   if (!test_allocator_deleter(seg_mgr, mapping))
      return false;

   if (!test_atomic_func(seg_mgr, mapping))
      return false;

   if (!test_allocator_deleter(seg_mgr, mapping))
      return false;
   
   if (!test_get_memory_algorithm(seg_mgr, mapping))
      return false;

   return true;
}

template<class MemoryAlgorithm>
bool test_each_algo()
{
   {
      typedef segment_manager< char, MemoryAlgorithm, flat_map_index > segment_manager_t;
      if(!test_segment_manager<segment_manager_t>())
         return false;
   }
   {
      typedef segment_manager< char, MemoryAlgorithm, map_index > segment_manager_t;
      if(!test_segment_manager<segment_manager_t>())
         return false;
   }
   {
      typedef segment_manager< char, MemoryAlgorithm, iset_index > segment_manager_t;
      if(!test_segment_manager<segment_manager_t>())
         return false;
   }
   {
      typedef segment_manager< char, MemoryAlgorithm, iunordered_set_index > segment_manager_t;
      if(!test_segment_manager<segment_manager_t>())
         return false;
   }
   return true;
}

int main()
{
   if(!test_each_algo< simple_seq_fit< null_mutex_family > >())
      return 1;
   if(!test_each_algo< rbtree_best_fit< null_mutex_family > >())
      return 1;

   return 0;
}
