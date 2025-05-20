//////////////////////////////////////////////////////////////////////////////
//
// (C) Copyright Ion Gaztanaga 2006-2012. Distributed under the Boost
// Software License, Version 1.0. (See accompanying file
// LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// See http://www.boost.org/libs/interprocess for documentation.
//
//////////////////////////////////////////////////////////////////////////////

#include <boost/interprocess/detail/workaround.hpp>
//[doc_cont
#include <boost/container/vector.hpp>
#include <boost/interprocess/allocators/allocator.hpp>
#include <boost/interprocess/managed_shared_memory.hpp>
//<-
#include "../test/get_process_id_name.hpp"
//->

int main ()
{
   using namespace boost::interprocess;
   //Remove shared memory on construction and destruction
   struct shm_remove
   {
      shm_remove() { shared_memory_object::remove(test::get_process_id_name()); }
      ~shm_remove(){ shared_memory_object::remove(test::get_process_id_name()); }
   } remover;
   //<-
   (void)remover;
   //->

   //A managed shared memory where we can construct objects
   //associated with a c-string
   managed_shared_memory segment(create_only,test::get_process_id_name(), 65536);

   //Alias an STL-like allocator of ints that allocates ints from the segment
   typedef allocator<int, managed_shared_memory::segment_manager>
      ShmemAllocator;

   //Alias a vector that uses the previous STL-like allocator
   typedef boost::container::vector<int, ShmemAllocator> MyVector;

   int initVal[]        = {0, 1, 2, 3, 4, 5, 6 };
   const int *begVal    = initVal;
   const int *endVal    = initVal + sizeof(initVal)/sizeof(initVal[0]);

   //Initialize the STL-like allocator
   const ShmemAllocator alloc_inst (segment.get_segment_manager());

   //Construct the vector in the shared memory segment with the STL-like allocator
   //from a range of iterators
   MyVector *myvector =
      segment.construct<MyVector>
         ("MyVector")/*object name*/
         (begVal     /*first ctor parameter*/
         ,endVal     /*second ctor parameter*/
         ,alloc_inst /*third ctor parameter*/);

   //Use vector as your want
   std::sort(myvector->rbegin(), myvector->rend());
   // . . .
   //When done, destroy and delete vector from the segment
   segment.destroy<MyVector>("MyVector");
   return 0;
}
//]

