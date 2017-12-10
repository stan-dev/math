#ifndef STAN_MATH_PRIM_ARR_FUNCTOR_MPI_DISTRIBUTED_APPLY_HPP
#define STAN_MATH_PRIM_ARR_FUNCTOR_MPI_DISTRIBUTED_APPLY_HPP

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>

#include <stan/math/prim/arr/functor/mpi_command.hpp>

namespace stan {
  namespace math {

    /**
     * MPI command template which calls the static method
     * distributed_apply of the given class F on remotes.
     *
     * @tparam F type of functor containing static distributed_apply method.
     */
    template<typename T>
    struct mpi_distributed_apply : public mpi_command {
      friend class boost::serialization::access;
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version) {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(mpi_command);
      }
      void run() const {
        T::distributed_apply();
      }
    };
    
  }
}

#define STAN_REGISTER_MPI_DISTRIBUTED_APPLY(APPLY_FUNCTOR) \
  STAN_REGISTER_MPI_COMMAND(stan::math::mpi_distributed_apply< APPLY_FUNCTOR >)

#endif
