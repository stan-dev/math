#ifdef STAN_MPI

#include <stan/math/prim/functor/mpi_cluster.hpp>

// register stop worker command (instantiates boost serialization
// templates)
STAN_REGISTER_MPI_COMMAND(stan::math::mpi_stop_worker)

#endif
