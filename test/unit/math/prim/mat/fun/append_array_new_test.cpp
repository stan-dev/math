#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>

using namespace Eigen;

/**
 * Generate a new random variable, cast it to type T1, and store it in z. n0
 * is ignored in this specialization
 *
 * @tparam T1 Scalar type
 * @param n0 Ignored
 * @param z Output variable
 */
template<typename T1>
void build(int n0, T1 &z) {
  z = T1(rand());
}

/**
 * Generate a new random numbers to fill z (casting them to type T1).
 *
 * If z is a column vector, first resize to have n0 rows
 *
 * If z is a row vector, first resize to have n0 columns
 *
 * If z is a matrix, first resize to have n0 rows and n0 columns
 *
 * @tparam T1 Scalar type
 * @param n0 Dimension
 * @param z Output variable
 */
template<typename T1, int R, int C>
void build(int n0, Eigen::Matrix<T1, R, C>& z) {
  int nrow = n0, ncol = n0;

  if(R == 1)
    nrow = 1;

  if(C == 1)
    ncol = 1;

  z = Eigen::Matrix<T1, R, C>::Random(nrow, ncol);
}

/**
 * Make a new std::vector of length n1 and recursively call build on
 * each element of it.
 *
 * @tparam T1 Element type
 * @param n1 Dimension of vector
 * @param n0 Dimension of last element
 * @param z Output variable
 */
template<typename T1>
void build(int n1, int n0, std::vector<T1> &z) {
  z.resize(n1);
  for(size_t i = 0; i < n1; i++) {
    build(n0, z[i]);
  }
}

/**
 * Make a new std::vector of length n2 and recursively call build on
 * each element of it.
 *
 * @tparam T1 Element type
 * @param n2 Dimension of std::vector
 * @param n1 Dimension of child std::vector
 * @param n0 Dimension of last element
 * @param z Output variable
 */
template<typename T1>
void build(int n2, int n1, int n0, std::vector<T1> &z) {
  z.resize(n2);
  for(size_t i = 0; i < n2; i++) {
    build(n1, n0, z[i]);
  }
}

/**
 * Check if variables are equal via floating point macro
 *
 * @tparam T1 Type of first argument
 * @tparam T2 Type of second argument
 * @param z1 First argument
 * @param z2 Second argument
 */
template<typename T1, typename T2>
void check_eq(const T1& z1, const T2& z2) {
  EXPECT_FLOAT_EQ(z1, z2);
}

/**
 * Check if two integers are equal
 *
 * @param z1 First integer
 * @param z2 Second integer
 */
template<>
void check_eq(const int& z1, const int& z2) {
  EXPECT_EQ(z1, z2);
}

/**
 * Check if elements of two matrices are equal via floating point macro
 *
 * @tparam T1 Scalar type of first matrix
 * @tparam T2 Scalar type of second matrix
 * @param z1 First matrix
 * @param z2 Second matrix
 */
template<typename T1, typename T2, int R, int C>
void check_eq(const Eigen::Matrix<T1, R, C>& z1, const Eigen::Matrix<T2, R, C>& z2) {
  EXPECT_EQ(z1.rows(), z2.rows());
  EXPECT_EQ(z1.cols(), z2.cols());

  for(size_t i = 0; i < z1.rows(); i++)
    for(size_t j = 0; j < z1.cols(); j++)
      check_eq(z1(i, j), z2(i, j));
}

/**
 * Recursively check if elements of two std::vectors are equal
 *
 * @tparam T1 Element type of first std::vector
 * @tparam T2 Element type of second std::vector
 * @param z1 First std::vector
 * @param z2 Second std::vector
 */
template<typename T1, typename T2>
void check_eq(const std::vector<T1>& z1, const std::vector<T2>& z2) {
  EXPECT_EQ(z1.size(), z2.size());
  for(size_t i = 0; i < z1.size(); i++)
    check_eq(z1[i], z2[i]);
}

/**
 * Randomly create std::vectors of type T1 and T2, append them together, assign
 * them to a std::vector of type T3, and then check that the copy and cast
 * happened correctly.
 *
 * @tparam T1 Element type of first std::vector
 * @tparam T2 Element type of second std::vector
 * @tparam T3 Element type of return std::vector
 */
template<typename T1, typename T2, typename T3>
void checkv() {
  std::vector<T1> x;
  std::vector<T2> y;
  std::vector<T3> result;

  int r1 = rand() % 5 + 1,
    r2 = rand() % 5 + 3,
    r3 = rand() % 5 + 3;

  build(r1, r2, x);
  build(r3, r2, y);
  EXPECT_NO_THROW(result = stan::math::append_array(x, y));
  EXPECT_EQ(x.size() + y.size(), result.size());
  for(size_t i = 0; i < x.size(); i++)
    check_eq(result[i], x[i]);
  for(size_t i = 0; i < y.size(); i++)
    check_eq(result[x.size() + i], y[i]);
}

/**
 * Randomly create std::vector<std::vector>s of type T1 and T2, append them
 * together, assign them to a std::vector of type T3, and then check that the
 * copy and cast happened correctly.
 *
 * @tparam T1 Element type of first std::vector
 * @tparam T2 Element type of second std::vector
 * @tparam T3 Element type of return std::vector
 */
template<typename T1, typename T2, typename T3>
void checkvv() {
  std::vector<std::vector<T1> > x;
  std::vector<std::vector<T2> > y;
  std::vector<std::vector<T3> > result;

  int r1 = rand() % 5 + 1,
    r2 = rand() % 5 + 3,
    r3 = rand() % 5 + 3,
    r4 = rand() % 5 + 3;

  build(r1, r2, r3, x);
  build(r4, r2, r3, y);
  EXPECT_NO_THROW(result = stan::math::append_array(x, y));
  EXPECT_EQ(x.size() + y.size(), result.size());
  for(size_t i = 0; i < x.size(); i++)
    check_eq(result[i], x[i]);
  for(size_t i = 0; i < y.size(); i++)
    check_eq(result[x.size() + i], y[i]);
}

TEST(MathFunctions, append_array_promote_matrix) {
  typedef Matrix<double, Dynamic, 1> V;
  typedef Matrix<double, 1, Dynamic> RV;
  typedef Matrix<double, Dynamic, Dynamic> M;
  checkv<V, V, V>();
  checkv<RV, RV, RV>();
  checkv<M, M, M>();
  checkv<double, double, double>();
  checkv<double, int, double>();
  checkv<int, double, double>();

  checkvv<V, V, V>();
  checkvv<RV, RV, RV>();
  checkvv<M, M, M>();
  checkvv<double, double, double>();
  checkvv<double, int, double>();
  checkvv<int, double, double>();
}
