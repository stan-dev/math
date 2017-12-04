#ifndef STAN_MATH_PRIM_ARR_FUNCTOR_MPI_COMMAND_HPP
#define STAN_MATH_PRIM_ARR_FUNCTOR_MPI_COMMAND_HPP

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/export.hpp>

namespace stan {
  namespace math {

    /**
     * Base class for all MPI commands which defines the virtual run
     * method. All MPI commands must inherit from this class and
     * implement the run method which defines what will be executed on
     * each remote.
     *
     * <p>Note that a concrete command must also register itself with
     * the serialization library using BOOST_CLASS_EXPORT export
     * macro. For optimal performance it is also recommended to use
     * the BOOST_CLASS_TRACKING to disable version tracking which is
     * not relevant in this context
     */
    struct mpi_command {
      friend class boost::serialization::access;
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version) {}
      virtual void run() const = 0;
    };

  }
}


BOOST_SERIALIZATION_ASSUME_ABSTRACT( stan::math::mpi_command )

#define STAN_REGISTER_MPI_COMMAND(command) \
  BOOST_CLASS_IMPLEMENTATION(command, boost::serialization::object_serializable) \
  BOOST_CLASS_EXPORT(command) \
  BOOST_CLASS_TRACKING(command,track_never)


#endif
