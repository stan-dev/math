//////////////////////////////////////////////////////////////////////////////
//
// (C) Copyright Ion Gaztanaga 2006-2012. Distributed under the Boost
// Software License, Version 1.0. (See accompanying file
// LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// See http://www.boost.org/libs/interprocess for documentation.
//
//////////////////////////////////////////////////////////////////////////////

//[doc_named_mutex
#include <boost/interprocess/sync/scoped_lock.hpp>
#include <boost/interprocess/sync/named_mutex.hpp>
#include <fstream>
#include <iostream>
#include <cstdio>

//<-
#include "../test/get_process_id_name.hpp"
//->

int main ()
{
   using namespace boost::interprocess;
   BOOST_INTERPROCESS_TRY{
      struct file_remove
      {
         file_remove() { std::remove(test::get_process_id_name()); }
         ~file_remove(){ std::remove(test::get_process_id_name()); }
      } file_remover;
      struct mutex_remove
      {
         mutex_remove() { named_mutex::remove(test::get_process_id_name()); }
         ~mutex_remove(){ named_mutex::remove(test::get_process_id_name()); }
      } remover;
      //<-
      (void)file_remover;
      (void)remover;
      //->

      //Open or create the named mutex
      named_mutex mutex(open_or_create, test::get_process_id_name());

      std::ofstream file(test::get_process_id_name());

      for(int i = 0; i < 10; ++i){

         //Do some operations...

         //Write to file atomically
         scoped_lock<named_mutex> lock(mutex);
         file << "Process name, ";
         file << "This is iteration #" << i;
         file << std::endl;
      }
   }
   BOOST_INTERPROCESS_CATCH(interprocess_exception &ex){
      std::cout << ex.what() << std::endl;
      return 1;
   } BOOST_INTERPROCESS_CATCH_END
   return 0;
}
//]

