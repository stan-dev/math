#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>

TEST(ErrorHandlingMatrix, checkStdVectorIndexMatrix) {
  std::vector<double> y;
  y.push_back(5);
  y.push_back(5);
  y.push_back(5);
  y.push_back(5);
  size_t i;
  
  i=2;
  y.resize(3);
  EXPECT_NO_THROW(stan::math::check_std_vector_index("checkStdVectorIndexMatrix",
                                                     "i", y, i));
  i=3;
  EXPECT_NO_THROW(stan::math::check_std_vector_index("checkStdVectorIndexMatrix",
                                                     "i", y, i));

  y.resize(2);
  EXPECT_THROW(stan::math::check_std_vector_index("checkStdVectorIndexMatrix",
                                                  "i", y, i),
               std::out_of_range);

  i=0;
  EXPECT_THROW(stan::math::check_std_vector_index("checkStdVectorIndexMatrix",
                                                  "i", y, i),
               std::out_of_range);
}

TEST(ErrorHandlingMatrix, checkStdVectorIndexMatrix_nan) {
  std::vector<double> y;
  double nan = std::numeric_limits<double>::quiet_NaN();
  y.push_back(nan);
  y.push_back(nan);
  y.push_back(nan);
  y.push_back(nan);
  size_t i;
  
  i=2;
  y.resize(3);
  EXPECT_NO_THROW(stan::math::check_std_vector_index("checkStdVectorIndexMatrix",
                                                     "i", y, i));
  i=3;
  EXPECT_NO_THROW(stan::math::check_std_vector_index("checkStdVectorIndexMatrix",
                                                     "i", y, i));

  y.resize(2);
  EXPECT_THROW(stan::math::check_std_vector_index("checkStdVectorIndexMatrix",
                                                  "i", y, i),
               std::out_of_range);

  i=0;
  EXPECT_THROW(stan::math::check_std_vector_index("checkStdVectorIndexMatrix",
                                                  "i", y, i),
               std::out_of_range);
}
