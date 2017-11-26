#include <stan/math.hpp>
#include <gtest/gtest.h>

#include <test/unit/math/prim/mat/functor/hard_work.hpp>
#include <test/unit/math/prim/mat/functor/faulty_functor.hpp>

#include <iostream>

#ifdef STAN_HAS_MPI

typedef stan::math::map_rect_reduce<hard_work, stan::math::var, stan::math::var> hard_work_reducer_vv;
typedef stan::math::map_rect_combine<hard_work, stan::math::var, stan::math::var> hard_work_combiner_vv;
typedef stan::math::mpi_parallel_call<hard_work_reducer_vv,hard_work_combiner_vv> hard_work_parallel_call;
BOOST_CLASS_EXPORT(stan::math::mpi_distributed_apply<hard_work_parallel_call>);
BOOST_CLASS_TRACKING(stan::math::mpi_distributed_apply<hard_work_parallel_call>,track_never);

typedef stan::math::map_rect_reduce<faulty_functor, stan::math::var, stan::math::var> faulty_functor_reducer_vv;
typedef stan::math::map_rect_combine<faulty_functor, stan::math::var, stan::math::var> faulty_functor_combiner_vv;
typedef stan::math::mpi_parallel_call<faulty_functor_reducer_vv,faulty_functor_combiner_vv> faulty_functor_parallel_call;
BOOST_CLASS_EXPORT(stan::math::mpi_distributed_apply<faulty_functor_parallel_call>);
BOOST_CLASS_TRACKING(stan::math::mpi_distributed_apply<faulty_functor_parallel_call>,track_never);

#endif

struct MpiJob : public ::testing::Test {
  stan::math::vector_v shared_params_v;
  std::vector<stan::math::vector_v> job_params_v;
  const std::size_t N = 10;
  std::vector<std::vector<double> > x_r = std::vector<std::vector<double>>(N, std::vector<double>(1,1.0));
  std::vector<std::vector<int> > x_i = std::vector<std::vector<int>>(N, std::vector<int>(1,0));

   virtual void SetUp() {
     shared_params_v.resize(2);
     shared_params_v << 2, 0;
     
     for(std::size_t n = 0; n != N; ++n) {
       x_i[n][0] = n;
       stan::math::vector_v job_v(2);
       job_v << n+1.0, n * n;
       job_params_v.push_back(job_v);
     }
   }
  
};

TEST_F(MpiJob, hard_work_vv) {
  stan::math::vector_v result = stan::math::map_rect<hard_work>(shared_params_v, job_params_v, x_r, x_i, 0);

  EXPECT_EQ(result.rows(), 2*N );
} 

TEST_F(MpiJob, always_faulty_functor_vv) {

  stan::math::vector_v result;

  EXPECT_NO_THROW(result = stan::math::map_rect<faulty_functor>(shared_params_v, job_params_v, x_r, x_i, 1));

  // faulty functor throws on theta(0) being -1.0
  // throwing during the first evaluation is quite severe and will
  // lead to a respective runtime error
  job_params_v[0](0) = -1;

  // upon the second evaluation throwing is handled internally different
  EXPECT_ANY_THROW(result = stan::math::map_rect<faulty_functor>(shared_params_v, job_params_v, x_r, x_i, 1));

  // thorwing on the very first evaluation
  EXPECT_ANY_THROW(result = stan::math::map_rect<faulty_functor>(shared_params_v, job_params_v, x_r, x_i, 2));

}

// from :
// http://www.parresianz.com/mpi/c++/mpi-unit-testing-googletests-cmake/
// but does not work well with our MPI design
/*
class MPIEnvironment : public ::testing::Environment
{
public:
  virtual void SetUp() {
    char** argv;
    int argc = 0;
    int mpiError = MPI_Init(&argc, &argv);
    ASSERT_FALSE(mpiError);
  }
  virtual void TearDown() {
    int mpiError = MPI_Finalize();
    ASSERT_FALSE(mpiError);
  }
  virtual ~MPIEnvironment() {}
};
*/

int main(int argc, char* argv[]) {
    int result = 0;
    ::testing::InitGoogleTest(&argc, argv);
    //::testing::AddGlobalTestEnvironment(new MPIEnvironment);

    stan::math::mpi_cluster cluster;

    //boost::mpi::environment env;
    //boost::mpi::communicator world;

    const std::size_t rank = cluster.rank_;

    // disable output listeners on all workers
    ::testing::TestEventListeners& listeners =
        ::testing::UnitTest::GetInstance()->listeners();
    if (rank != 0) {
      delete listeners.Release(listeners.default_result_printer());
    }

    // send workers into listen mode
    cluster.listen();
    
    // run all tests on the root
    if (rank == 0) {
      result = RUN_ALL_TESTS();
    }

    return result;
}
