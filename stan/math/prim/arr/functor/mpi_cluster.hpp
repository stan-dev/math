#ifndef STAN_MATH_PRIM_ARR_FUNCTOR_MPI_CLUSTER_HPP
#define STAN_MATH_PRIM_ARR_FUNCTOR_MPI_CLUSTER_HPP


#include <stan/math/prim/arr/functor/mpi_command.hpp>

#include <boost/mpi.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/shared_ptr.hpp>

#include <mutex>

namespace stan {
  namespace math {

    // MPI command which will shut down a child gracefully
    class mpi_stop_listen : public std::exception {
      virtual const char* what() const throw()
      {
        return "Stopping MPI listening mode.";
      }
    };
      
    struct mpi_stop_worker : public mpi_command {
      friend class boost::serialization::access;
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version) {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(mpi_command);
      }
      void run() const {
        boost::mpi::communicator world;
        std::cout << "Terminating worker " << world.rank() << std::endl;
        //MPI_Finalize();
        throw mpi_stop_listen();
        //std::exit(0);
      }
    };

    // map jobs of given chunk size to workers; used for deterministic
    // scheduling
    std::vector<int>
    mpi_map_chunks(std::size_t num_jobs, std::size_t chunk_size = 1) {
      boost::mpi::communicator world;
      const std::size_t world_size = world.size();
        
      std::vector<int> chunks(world_size, num_jobs / world_size);
      
      for(std::size_t r = 0; r != num_jobs % world_size; r++)
        ++chunks[r+1];

      for(std::size_t i = 0; i != world_size; i++)
        chunks[i] *= chunk_size;

      return(chunks);
    }

    template<typename T>
    void mpi_broadcast_command();

    struct mpi_cluster {
      boost::mpi::environment env;
      boost::mpi::communicator world_;
      std::size_t const rank_ = world_.rank();
      static std::mutex command_mutex;
      
      mpi_cluster() : cluster_listens_(false) {}

      ~mpi_cluster() {
        // the destructor will ensure that the childs are being
        // shutdown
        stop_listen();
      }

      void listen() {
        cluster_listens_ = true;
        if(rank_ == 0) {
          return;
        }
        std::cout << "Worker " << rank_ << " listening for commands..." << std::endl;
        
        try {
          while(1) {
            boost::shared_ptr<mpi_command> work;
            
            boost::mpi::broadcast(world_, work, 0);

            work->run();
          }
        } catch(const mpi_stop_listen& e) {
          std::cout << "Wrapping up MPI on worker " << rank_ << " ..." << std::endl;
        }
      }

      void stop_listen() {
        if(rank_ == 0 && cluster_listens_) {
          mpi_broadcast_command<mpi_stop_worker>();
        }
        cluster_listens_ = false;
     }

    private:
      bool cluster_listens_;
    };

    std::mutex mpi_cluster::command_mutex;

    template<typename T>
    void mpi_broadcast_command() {
      boost::mpi::communicator world;
      
      if(world.rank() != 0)
        throw std::runtime_error("only root may broadcast commands.");

      // used to lock the mpi cluster during execution of some
      // distributed task
      std::lock_guard<std::mutex> lock_cluster(mpi_cluster::command_mutex);
      
      boost::shared_ptr<mpi_command> command(new T);
      
      boost::mpi::broadcast(world, command, 0);
    }

  }
}

STAN_REGISTER_MPI_COMMAND(stan::math::mpi_stop_worker)

#endif
