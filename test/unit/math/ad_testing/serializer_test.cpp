#include <stan/math/prim/fun/serializer.hpp>
#include <gtest/gtest.h>
#include <complex>
#include <vector>

TEST(testUnitMathSerializer, serializer_deserializer) {
  stan::math::serializer<double> s;

  s.write(3.2);
  s.write(-1);

  s.write(std::complex<double>(-1.5, 15.75));

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

  std::vector<double> expected{3.2, -1,  -1.5, 15.75, 10, 20, 30, -10, -20,
                               101, 102, 103,  1,     3,  5,  2,  4,   6};
  for (size_t i = 0; i < expected.size(); ++i)
    EXPECT_EQ(expected[i], s.vals_[i]);

  stan::math::deserializer<double> d = stan::math::to_deserializer(s.vals_);

  EXPECT_EQ(3.2, d.read(0.0));
  EXPECT_EQ(-1, d.read(0.0));
  std::complex<double> w = d.read(std::complex<double>(1, 1));
  EXPECT_EQ(-1.5, w.real());
  EXPECT_EQ(15.75, w.imag());
  std::vector<double> x = d.read(std::vector<double>{0, 0, 0});
  EXPECT_EQ(10, x[0]);
  EXPECT_EQ(20, x[1]);
  EXPECT_EQ(30, x[2]);
  Eigen::VectorXd y = d.read(Eigen::VectorXd(2));
  EXPECT_EQ(-10, y(0));
  EXPECT_EQ(-20, y(1));
  Eigen::RowVectorXd z = d.read(Eigen::RowVectorXd(3));
  EXPECT_EQ(101, z(0));
  EXPECT_EQ(102, z(1));
  EXPECT_EQ(103, z(2));
  Eigen::MatrixXd u = d.read(Eigen::MatrixXd(3, 2));
  EXPECT_EQ(1, u(0, 0));
  EXPECT_EQ(2, u(0, 1));
  EXPECT_EQ(3, u(1, 0));
  EXPECT_EQ(4, u(1, 1));
  EXPECT_EQ(5, u(2, 0));
  EXPECT_EQ(6, u(2, 1));
}

TEST(testUnitMathSerializer, serialize) {
  std::vector<double> xs = stan::math::serialize<double>();
  EXPECT_EQ(0, xs.size());

  double a = 2;
  std::complex<double> b(1, 2);
  std::vector<double> c{3, 4, 5};
  Eigen::MatrixXd d(2, 3);
  d << -1, -2, -3, -4, -5, -6;
  std::vector<double> ys = stan::math::serialize<double>(a, b, c, d);

  std::vector<double> expected{2, 1, 2, 3, 4, 5, -1, -4, -2, -5, -3, -6};
  for (size_t i = 0; i < expected.size(); ++i)
    EXPECT_EQ(expected[i], ys[i]);
}
