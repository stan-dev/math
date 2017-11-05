#pragma once

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/shared_ptr.hpp>

/* define virtual class which gets send around by mpi. The run method
   contains the actual work to be executed.
 */

namespace stan {
  namespace math {
    
    struct mpi_command {
      friend class boost::serialization::access;
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version) {}
      virtual void run() const = 0;
    };

  }
}


BOOST_SERIALIZATION_ASSUME_ABSTRACT( stan::math::mpi_command )

