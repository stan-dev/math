#ifndef STAN_MATH_MPI_ENVIONMENT_HPP
#define STAN_MATH_MPI_ENVIONMENT_HPP

#include <stan/math/prim/err/check_greater_or_equal.hpp>
#ifdef STAN_LANG_MPI
#include <boost/mpi.hpp>
#endif

// default comm to world comm, in case stan needs to be
// called as library.
#define MPI_COMM_STAN MPI_COMM_WORLD

namespace stan {
namespace math {
  namespace mpi {

    /*
     * MPI Evionment that initializes and finalizes the MPI
     */
    struct Envionment {
        Envionment() {
          init();
        }
        ~Envionment() {
          finalize();
        }

        static void init() {
#ifdef STAN_LANG_MPI
          int flag;
          MPI_Initialized(&flag);
          if(!flag) {
            int provided;
            MPI_Init_thread(NULL, NULL, MPI_THREAD_SINGLE, &provided);
            // print provided when needed
          }
#endif
        }

        static void finalize() {
#ifdef STAN_LANG_MPI
          int flag;
          MPI_Finalized(&flag);
          if(!flag) MPI_Finalize();
#endif
        }
    };

    /*
     * MPI Communicators. With default constructor disabled,
     * a communicator can only be created through duplication.
     */
    struct Communicator {
    private:
      Communicator();

    public:
      MPI_Comm c;
      int s;
      int r;

      /*
       * communicator constructor using @c Envionment and @c MPI_Comm
       */
      explicit Communicator(MPI_Comm other) :
        c(MPI_COMM_NULL), s(0), r(-1) {
        if (other != MPI_COMM_NULL) {
          MPI_Comm_dup(other, &c);
          MPI_Comm_size(c, &s);
          MPI_Comm_rank(c, &r);        
        }
      }

      /*
       * copy constructor is deep
       */
      explicit Communicator(const Communicator& other) :
        Communicator(other.comm())
      {}

      /*
       * type-cast to MPI_Comm object
       */
      operator MPI_Comm() {
        return this -> c;
      }

      const MPI_Comm& comm() const {
        return c;
      }
      
      int size() const {return s;}

      int rank() const {return r;}

      void set_comm(const MPI_Comm& other) {
        if (other != MPI_COMM_NULL) {
          c = other;
          MPI_Comm_size(c, &s);
          MPI_Comm_rank(c, &r);
        }
      }

      /*
       * destructor needs to free MPI_Comm
       */
      ~Communicator() {
        if (c != MPI_COMM_NULL) {
          MPI_Comm_free(&c);
        }
      }
    };
  
    /* MPI communicator wrapper for RAII. Note that no
     * MPI's predfined comm such as @c MPI_COMM_WOLRD are allowed.*/
    struct Session {
      static Envionment env;
      static MPI_Comm MPI_COMM_INTER_CHAIN;
      static MPI_Comm MPI_COMM_INTRA_CHAIN;
      static Communicator inter_chain;
      static Communicator intra_chain;

      /*
       * Return inter-chain communicator. When called first
       * time, the communicator is populated with proper
       * members.
       */
      static const MPI_Comm& mpi_comm_inter_chain(int num_mpi_chains) {
        if(MPI_COMM_INTER_CHAIN == MPI_COMM_NULL) {
          int world_size;
          MPI_Comm_size(MPI_COMM_STAN, &world_size);
          stan::math::check_greater_or_equal("MPI inter-chain session",
                                             "number of procs", world_size,
                                             num_mpi_chains);

          MPI_Group stan_group, new_group;
          MPI_Comm_group(MPI_COMM_STAN, &stan_group);
          int num_chain_with_extra_proc = world_size % num_mpi_chains;
          int num_proc_per_chain = world_size / num_mpi_chains;
          std::vector<int> ranks(num_mpi_chains);
          if (num_chain_with_extra_proc == 0) {
            for (int i = 0, j = 0; i < world_size; i += num_proc_per_chain, ++j) {
              ranks[j] = i;
            }
          } else {
            num_proc_per_chain++;
            int i = 0;
            for (int j = 0; j < num_chain_with_extra_proc; ++j) {
              ranks[j] = i;
              i += num_proc_per_chain;
            }
            num_proc_per_chain--;
            for (int j = num_chain_with_extra_proc; j < num_mpi_chains; ++j) {
              ranks[j] = i;
              i += num_proc_per_chain;
            }
          }

          MPI_Group_incl(stan_group, num_mpi_chains, ranks.data(), &new_group);
          MPI_Comm_create_group(MPI_COMM_STAN, new_group, 99, &MPI_COMM_INTER_CHAIN);
          MPI_Group_free(&new_group);
          MPI_Group_free(&stan_group);
        }
        
        return MPI_COMM_INTER_CHAIN;                                                
      }

      /*
       * Return intra-chain communicator. When called first
       * time, the communicator is populated with proper
       * members.
       */
      static const MPI_Comm& mpi_comm_intra_chain(int num_mpi_chains) {
        // Envionment::env.init();

        if (MPI_COMM_INTRA_CHAIN == MPI_COMM_NULL) {
          int world_size, world_rank, color;
          MPI_Comm_size(MPI_COMM_STAN, &world_size);
          MPI_Comm_rank(MPI_COMM_STAN, &world_rank);

          int num_chain_with_extra_proc = world_size % num_mpi_chains;
          const int n_proc = world_size / num_mpi_chains;
          if (num_chain_with_extra_proc == 0) {
            color = world_rank / n_proc;
          } else {
            int i = 0;
            for (int j = 0; j < num_mpi_chains; ++j) {
              const int n = j < num_chain_with_extra_proc ? (n_proc + 1) : n_proc;
              if (world_rank >= i && world_rank < i + n) {
                color = i;
                break;
              }
              i += n;
            }
          }
          MPI_Comm_split(MPI_COMM_STAN, color, world_rank, &MPI_COMM_INTRA_CHAIN);
        }

        return MPI_COMM_INTRA_CHAIN;
      }

      static const Communicator& inter_chain_comm(int num_mpi_chains) {
        if (inter_chain.comm() == MPI_COMM_NULL && intra_chain.comm() == MPI_COMM_NULL) {
          mpi_comm_inter_chain(num_mpi_chains);
          mpi_comm_intra_chain(num_mpi_chains);
          if (MPI_COMM_INTER_CHAIN != MPI_COMM_NULL) {
            inter_chain.set_comm(MPI_COMM_INTER_CHAIN);
          }
        }
        return inter_chain;
      }
        
      static const Communicator& intra_chain_comm(int num_mpi_chains) {
        if (intra_chain.comm() == MPI_COMM_NULL) {
          mpi_comm_inter_chain(num_mpi_chains);
          mpi_comm_intra_chain(num_mpi_chains);
          intra_chain.set_comm(MPI_COMM_INTRA_CHAIN);
        }
        return intra_chain;
      }
        
      static bool is_in_inter_chain_comm(int num_mpi_chains) {
        return inter_chain_comm(num_mpi_chains).rank() >= 0;
      }

    };

    Envionment Session::env;
    MPI_Comm Session::MPI_COMM_INTER_CHAIN(MPI_COMM_NULL);
    MPI_Comm Session::MPI_COMM_INTRA_CHAIN(MPI_COMM_NULL);
    Communicator Session::inter_chain(MPI_COMM_NULL);
    Communicator Session::intra_chain(MPI_COMM_NULL);
  }
}
}

#endif
