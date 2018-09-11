#include <stan/math/prim/arr.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathFunctions, asc_lin_seq_by_array_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_THROW(stan::math::asc_lin_seq_by_array(nan, 1, 1),
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(nan, 1, 1, true),
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(nan, 1, .1), 
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(nan, 1, .1, true),
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(nan, .1, 1),
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(nan, .1, .1, true),
               std::domain_error);

  EXPECT_THROW(stan::math::asc_lin_seq_by_array(0, nan, 1),
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(0, nan, 1, true),
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(0, nan, .1),
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(0, nan, .1, true),
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(1.1, nan, .1),
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(1.1, nan, .1, true), 
               std::domain_error);

  EXPECT_THROW(stan::math::asc_lin_seq_by_array(0, 1, nan), 
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(0, 1, nan, true), 
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(0, 1.1, nan), 
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(0, 1.1, nan, true), 
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(-.1, 1.1, nan), 
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(-.1, 1.1, nan, true), 
               std::domain_error);
}

TEST(MathFunctions, asc_lin_seq_by_array_inf) {
  double inf = std::numeric_limits<double>::infinity();

  EXPECT_THROW(stan::math::asc_lin_seq_by_array(-inf, 1, 2), 
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(-inf, 1, 2, true), 
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(-inf, 1.1, 2), 
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(-inf, 1.1, 2, true), 
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(-inf, 1.1, 2.1), 
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(-inf, 1.1, 2.1, true), 
               std::domain_error);

  EXPECT_THROW(stan::math::asc_lin_seq_by_array(1, inf, 1), 
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(1, inf, 1, true), 
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(.1, inf, 1), 
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(.1, inf, 1, true), 
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(1, inf, .1), 
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(1, inf, .1, true), 
               std::domain_error);

  EXPECT_THROW(stan::math::asc_lin_seq_by_array(1, 2, inf), 
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(1, 2, inf, true), 
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(1.1, 2, inf), 
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(1.1, 2, inf, true), 
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(1.1, 2.1, inf), 
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(1.1, 2.1, inf, true), 
               std::domain_error);
}

TEST(MathFunctions, asc_lin_seq_by_array_err) {
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(1, 2, -2), 
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(1, 2, -2, true),
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(1.1, 2, -2),
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(1.1, 2, -2, true), 
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(1, 2.1, -2), 
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(1, 2.1, -2, true), 
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(1, 2, -2.1), 
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(1, 2, -2.1, true), 
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(1.1, 2.1, -2), 
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(1.1, 2.1, -2, true), 
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(1, 2.1, -2), 
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(1, 2.1, -2.1, true),
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(1.1, 2.1, -2), 
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(1.1, 2.1, -2.1, true),
               std::domain_error);
 
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(1, -1, 1), 
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(1, -1, 1, true), 
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(.1, -1, 1), 
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(.1, -1, 1, true),
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(1, -.1, 1),
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(1, -.1, 1, true),
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(1, -1, .1), 
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(1, -1, .1, true),
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(.1, -.1, 1),
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(.1, -.1, 1, true), 
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(.1, -1, .1), 
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(.1, -1, .1, true), 
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(1, -.1, .1), 
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(1, -.1, .1, true), 
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(.1, -.1, .1), 
               std::domain_error);
  EXPECT_THROW(stan::math::asc_lin_seq_by_array(.1, -.1, .1, true), 
               std::domain_error);
}

TEST(MathFunctions, asc_lin_seq_by_array_int_size) {
  using std::vector;

  vector<int> test_v1 = stan::math::asc_lin_seq_by_array(1, 9, 2); 
  vector<int> test_v2 = stan::math::asc_lin_seq_by_array(1, 9, 2, false); 
  vector<int> test_v3 = stan::math::asc_lin_seq_by_array(1, 9, 2, true);
  vector<int> test_v4 = stan::math::asc_lin_seq_by_array(1, 8, 2);
  vector<int> test_v5 = stan::math::asc_lin_seq_by_array(1, 8, 2, true);
 
  EXPECT_EQ(4, test_v1.size());
  EXPECT_EQ(4, test_v2.size());
  EXPECT_EQ(5, test_v3.size());
  EXPECT_EQ(4, test_v4.size());
  EXPECT_EQ(5, test_v5.size());
}

TEST(MathFunctions, asc_lin_seq_by_array_double_size) {
  using std::vector;

  vector<double> test_v1 = stan::math::asc_lin_seq_by_array(1.1, 9, 2); 
  vector<double> test_v2 = stan::math::asc_lin_seq_by_array(1.1, 9.2, 2); 
  vector<double> test_v3 = stan::math::asc_lin_seq_by_array(1.1, 9.2, 2.0); 
  vector<double> test_v5 = stan::math::asc_lin_seq_by_array(1.1, 9,
                                                            2, true);
  vector<double> test_v6 = stan::math::asc_lin_seq_by_array(1, 8.3, 2);
  vector<double> test_v7 = stan::math::asc_lin_seq_by_array(1, 8.3, 2.0);
  vector<double> test_v8 = stan::math::asc_lin_seq_by_array(1.1, 8.3,
                                                            2, true);
  vector<double> test_v9 = stan::math::asc_lin_seq_by_array(1.1, 8.3,
                                                            2.0, true);
  EXPECT_EQ(4, test_v1.size());
  EXPECT_EQ(5, test_v2.size());
  EXPECT_EQ(5, test_v3.size());
  EXPECT_EQ(5, test_v5.size());
  EXPECT_EQ(4, test_v6.size());
  EXPECT_EQ(4, test_v7.size());
  EXPECT_EQ(5, test_v8.size());
  EXPECT_EQ(5, test_v9.size());
}

TEST(MathFunctions, asc_lin_seq_by_array_int) {
  using std::vector;

  vector<int> v1 = {1, 3, 5, 7};
  vector<int> v2 = {1, 3, 5, 7, 9};
  vector<int> v3 = {1, 3, 5, 7, 8};

  vector<int> test_v1 = stan::math::asc_lin_seq_by_array(1, 9, 2); 
  vector<int> test_v2 = stan::math::asc_lin_seq_by_array(1, 9, 2, true);
  vector<int> test_v3a = stan::math::asc_lin_seq_by_array(1, 8, 2);
  vector<int> test_v3b = stan::math::asc_lin_seq_by_array(1, 8, 2, true);

  for (size_t i = 0; i < v1.size(); ++i) {
    EXPECT_EQ(v1[i], test_v1[i]);
    EXPECT_EQ(v1[i], test_v3a[i]);
  }
  for (size_t i = 0; i < v2.size(); ++i) {
    EXPECT_EQ(v2[i], test_v2[i]);
  }
  for (size_t i = 0; i < v3.size(); ++i) {
    EXPECT_EQ(v3[i], test_v3b[i]);
  }
}

TEST(MathFunctions, asc_lin_seq_by_array_double) {
  using std::vector;

  vector<double> v1 = {1.1, 3.1, 5.1, 7.1};
  vector<double> v2 = {1.0, 3.0, 5.0, 7.0};
  vector<double> v3 = {1.1, 3.1, 5.1, 7.1, 9};
  vector<double> v4 = {1, 3, 5, 7, 8.5};
  vector<double> v5 = {1.0, 3.0, 5, 7, 8};
  vector<double> v6 = {1.1, 3.1, 5.1, 7.1, 9.1};

  vector<double> test_v1a = stan::math::asc_lin_seq_by_array(1.1, 9, 2); 
  vector<double> test_v1b = stan::math::asc_lin_seq_by_array(1.1, 9, 2.0); 
  vector<double> test_v1c = stan::math::asc_lin_seq_by_array(1.1, 8.5, 2); 
  vector<double> test_v1d = stan::math::asc_lin_seq_by_array(1.1, 8.5,
                                                             2.0); 
  vector<double> test_v2a = stan::math::asc_lin_seq_by_array(1, 8.5, 2); 
  vector<double> test_v2b = stan::math::asc_lin_seq_by_array(1, 8.5, 2.0); 
  vector<double> test_v2c = stan::math::asc_lin_seq_by_array(1, 9, 2.0); 
  vector<double> test_v3a = stan::math::asc_lin_seq_by_array(1.1, 9,
                                                             2, true);
  vector<double> test_v3b = stan::math::asc_lin_seq_by_array(1.1, 9,
                                                             2.0, true);
  vector<double> test_v4a = stan::math::asc_lin_seq_by_array(1, 8.5,
                                                             2, true);
  vector<double> test_v4b = stan::math::asc_lin_seq_by_array(1, 8.5,
                                                             2.0, true);
  vector<double> test_v5 = stan::math::asc_lin_seq_by_array(1, 8,
                                                            2.0, true);
  vector<double> test_v6 = stan::math::asc_lin_seq_by_array(1.1, 9.1,
                                                            2.0, true);
  for (size_t i = 0; i < v1.size(); ++i) {
    EXPECT_EQ(v1[i], test_v1a[i]);
    EXPECT_EQ(v1[i], test_v1b[i]);
    EXPECT_EQ(v1[i], test_v1c[i]);
    EXPECT_EQ(v1[i], test_v1d[i]);
  }
  for (size_t i = 0; i < v2.size(); ++i) {
    EXPECT_EQ(v2[i], test_v2a[i]);
    EXPECT_EQ(v2[i], test_v2b[i]);
    EXPECT_EQ(v2[i], test_v2c[i]);
  }
  for (size_t i = 0; i < v3.size(); ++i) {
    EXPECT_EQ(v3[i], test_v3a[i]);
    EXPECT_EQ(v3[i], test_v3b[i]);
  }
  for (size_t i = 0; i < v4.size(); ++i) {
    EXPECT_EQ(v4[i], test_v4a[i]);
    EXPECT_EQ(v4[i], test_v4b[i]);
  }
  for (size_t i = 0; i < v5.size(); ++i) {
    EXPECT_EQ(v5[i], test_v5[i]);
  }
  for (size_t i = 0; i < v6.size(); ++i) {
    EXPECT_EQ(v6[i], test_v6[i]);
  }
}
