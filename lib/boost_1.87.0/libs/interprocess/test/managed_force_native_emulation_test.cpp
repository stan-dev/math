//////////////////////////////////////////////////////////////////////////////
//
// (C) Copyright Ion Gaztanaga 2024-2024. Distributed under the Boost
// Software License, Version 1.0. (See accompanying file
// LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// See http://www.boost.org/libs/interprocess for documentation.
//
//////////////////////////////////////////////////////////////////////////////

#define BOOST_INTERPROCESS_FORCE_NATIVE_EMULATION 1 // This is important, doesn't crash otherwise.

#include <boost/interprocess/managed_shared_memory.hpp>
#include "get_process_id_name.hpp"

using namespace boost::interprocess;

int main()
{
   const int MemSize = 65536;
   const char* const MemName = test::get_process_id_name();

   {
      managed_shared_memory mshm(open_or_create, MemName, MemSize);
      shared_memory_object::remove(MemName);

      for (std::size_t i = 0; i != 1000; ++i) {
         mshm.deallocate(mshm.allocate(i));
      }
   }  //Force early destruction of shared memory

   return 0;
}

