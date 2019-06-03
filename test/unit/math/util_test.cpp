#include <test/unit/math/util.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(test_unit_math_test_ad, to_std_vector) {
  Eigen::VectorXd u(0);
  EXPECT_EQ(0, stan::test::to_std_vector(u).size());

  Eigen::VectorXd x(3);
  x << 1, 2, 3;
  std::vector<double> y = stan::test::to_std_vector(x);
  EXPECT_EQ(3, y.size());
  for (int i = 0; i < 3; ++i)
    EXPECT_EQ(x(i), y[i]);
}

TEST(test_unit_math_test_ad, to_eigen_vector) {
  std::vector<double> u;
  EXPECT_EQ(0, stan::test::to_eigen_vector(u).size());

  std::vector<double> v{1, 2, 3};
  Eigen::VectorXd vv = stan::test::to_eigen_vector(v);
  EXPECT_EQ(3, vv.size());
  for (int i = 0; i < 3; ++i)
    EXPECT_EQ(v[i], vv(i));
}

TEST(test_unit_math_util, serializer_deserializer) {
  stan::test::serializer<double> s;

  s.write(3.2);
  s.write(-1);
  s.write(std::vector<double>{10, 20, 30});

  Eigen::VectorXd a(2);
  a << -10, -20;
  s.write(a);

  Eigen::RowVectorXd b(3);
  b << 101, 102, 103;
  s.write(b);

  Eigen::MatrixXd c(3, 2);
  c << 1, 2, 3, 4, 5, 6;  // << is row major; index order is col major
  s.write(c);

  std::vector<double> expected{3.2, -1,  10, 20, 30, -10, -20, 101,
                               102, 103, 1,  3,  5,  2,   4,   6};
  for (size_t i = 0; i < expected.size(); ++i)
    EXPECT_EQ(expected[i], s.vals_[i]);

  stan::test::deserializer<double> d = stan::test::to_deserializer(s.vals_);

  EXPECT_EQ(3.2, d.read(0.0));
  EXPECT_EQ(-1, d.read(0.0));
  std::vector<double> x = d.read(std::vector<double>{0, 0, 0});
  EXPECT_EQ(10, x[0]);
  EXPECT_EQ(20, x[1]);
  EXPECT_EQ(30, x[2]);
  Eigen::VectorXd y = d.read(Eigen::VectorXd(2));
  EXPECT_EQ(-10, y(0));
  EXPECT_EQ(-20, y(1));
  Eigen::RowVectorXd z = d.read(Eigen::RowVectorXd(3));
  EXPECT_EQ(-10, y(0));
  EXPECT_EQ(-20, y(1));
  Eigen::MatrixXd u = d.read(Eigen::MatrixXd(3, 2));
  EXPECT_EQ(1, u(0, 0));
  EXPECT_EQ(2, u(0, 1));
  EXPECT_EQ(3, u(1, 0));
  EXPECT_EQ(4, u(1, 1));
  EXPECT_EQ(5, u(2, 0));
  EXPECT_EQ(6, u(2, 1));
}

TEST(test_unit_math_util, serialize) {
  std::vector<double> xs = stan::test::serialize<double>();
  EXPECT_EQ(0, xs.size());

  double a = 2;
  std::vector<double> b{3, 4, 5};
  Eigen::MatrixXd c(2, 3);
  c << -1, -2, -3, -4, -5, -6;
  std::vector<double> ys = stan::test::serialize<double>(a, b, c);

  std::vector<double> expected{2, 3, 4, 5, -1, -4, -2, -5, -3, -6};
  for (size_t i = 0; i < expected.size(); ++i)
    EXPECT_EQ(expected[i], ys[i]);
}
