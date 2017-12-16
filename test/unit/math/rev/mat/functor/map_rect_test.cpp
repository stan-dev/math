#include <stan/math.hpp>
#include <gtest/gtest.h>

#include <test/unit/math/prim/mat/functor/hard_work.hpp>
#include <test/unit/math/prim/mat/functor/faulty_functor.hpp>

#include <iostream>

struct MapRectJob : public ::testing::Test {
  stan::math::vector_v shared_params_v;
  std::vector<stan::math::vector_v> job_params_v;
  stan::math::vector_v shared_params_v2;
  std::vector<stan::math::vector_v> job_params_v2;
  const std::size_t N = 10;
  std::vector<std::vector<double> > x_r = std::vector<std::vector<double>>(N, std::vector<double>(1,1.0));
  std::vector<std::vector<int> > x_i = std::vector<std::vector<int>>(N, std::vector<int>(1,0));

   virtual void SetUp() {
     shared_params_v.resize(2);
     shared_params_v << 2, 0;
     shared_params_v2.resize(2);
     shared_params_v2 << 2, 0;
     
     for(std::size_t n = 0; n != N; ++n) {
       x_i[n][0] = n;
       stan::math::vector_v job_v(2);
       job_v << n+1.0, n * n;
       job_params_v.push_back(job_v);

       stan::math::vector_v job_v2(2);
       job_v2 << n+1.0, n * n;
       job_params_v2.push_back(job_v2);
     }
   }
  
};

TEST_F(MapRectJob, hard_work_vv) {
  
  stan::math::vector_v result_serial = stan::math::map_rect_serial<0,hard_work>(shared_params_v, job_params_v, x_r, x_i, 0);

}

TEST_F(MapRectJob, always_faulty_functor_vv) {

  stan::math::vector_v result;

  EXPECT_NO_THROW((result = stan::math::map_rect<1,faulty_functor>(shared_params_v, job_params_v, x_r, x_i)));

  // faulty functor throws on theta(0) being -1.0
  // throwing during the first evaluation is quite severe and will
  // lead to a respective runtime error
  job_params_v[0](0) = -1;

  // upon the second evaluation throwing is handled internally different
  EXPECT_ANY_THROW((result = stan::math::map_rect<1,faulty_functor>(shared_params_v, job_params_v, x_r, x_i)));

  // throwing on the very first evaluation 
  EXPECT_ANY_THROW((result = stan::math::map_rect<2,faulty_functor>(shared_params_v, job_params_v, x_r, x_i)));
}
