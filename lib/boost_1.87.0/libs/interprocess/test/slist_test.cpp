//////////////////////////////////////////////////////////////////////////////
//
// (C) Copyright Ion Gaztanaga 2004-2012. Distributed under the Boost
// Software License, Version 1.0. (See accompanying file
// LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// See http://www.boost.org/libs/interprocess for documentation.
//
//////////////////////////////////////////////////////////////////////////////
#include <boost/interprocess/managed_shared_memory.hpp>
#include <boost/container/slist.hpp>
#include <boost/interprocess/allocators/allocator.hpp>
#include <boost/interprocess/offset_ptr.hpp>
#include "dummy_test_allocator.hpp"
#include "list_test.hpp"
#include "movable_int.hpp"
#include "emplace_test.hpp"

using namespace boost::interprocess;

typedef allocator<int, managed_shared_memory::segment_manager> ShmemAllocator;
typedef boost::container::slist<int, ShmemAllocator> MyList;



typedef allocator<test::movable_int, managed_shared_memory::segment_manager> ShmemMoveAllocator;
typedef boost::container::slist<test::movable_int, ShmemMoveAllocator> MyMoveList;

typedef allocator<test::movable_and_copyable_int, managed_shared_memory::segment_manager> ShmemCopyMoveAllocator;
typedef boost::container::slist<test::movable_and_copyable_int, ShmemCopyMoveAllocator> MyCopyMoveList;

typedef allocator<test::copyable_int, managed_shared_memory::segment_manager> ShmemCopyAllocator;
typedef boost::container::slist<test::copyable_int, ShmemCopyAllocator> MyCopyList;

int main ()
{
   if(test::list_test<managed_shared_memory, MyList, false>())
      return 1;
   
   if(test::list_test<managed_shared_memory, MyMoveList, false>())
      return 1;

   if(test::list_test<managed_shared_memory, MyCopyMoveList, false>())
      return 1;

   if(test::list_test<managed_shared_memory, MyCopyList, false>())
      return 1;

   const test::EmplaceOptions Options = (test::EmplaceOptions)
      (test::EMPLACE_FRONT | test::EMPLACE_AFTER | test::EMPLACE_BEFORE  | test::EMPLACE_AFTER);

   if(!boost::interprocess::test::test_emplace
      < boost::container::slist<test::EmplaceInt>, Options>())
      return 1;
   return 0;
}

