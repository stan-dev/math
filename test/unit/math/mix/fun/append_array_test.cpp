#include <stan/math/mix.hpp>
#include <gtest/gtest.h>
#include <vector>
#include <random>

using stan::math::fvar;
using stan::math::var;

typedef fvar<double> fd;
typedef fvar<var> fv;

typedef fvar<fvar<double> > ffd;
typedef fvar<fvar<var> > ffv;

typedef Eigen::Matrix<double, Eigen::Dynamic, 1> V;
typedef Eigen::Matrix<double, 1, Eigen::Dynamic> RV;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> M;

typedef Eigen::Matrix<var, Eigen::Dynamic, 1> Vv;
typedef Eigen::Matrix<var, 1, Eigen::Dynamic> RVv;
typedef Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> Mv;

typedef Eigen::Matrix<fd, Eigen::Dynamic, 1> Vfd;
typedef Eigen::Matrix<fd, 1, Eigen::Dynamic> RVfd;
typedef Eigen::Matrix<fd, Eigen::Dynamic, Eigen::Dynamic> Mfd;

typedef Eigen::Matrix<ffd, Eigen::Dynamic, 1> Vffd;
typedef Eigen::Matrix<ffd, 1, Eigen::Dynamic> RVffd;
typedef Eigen::Matrix<ffd, Eigen::Dynamic, Eigen::Dynamic> Mffd;

typedef Eigen::Matrix<fv, Eigen::Dynamic, 1> Vfv;
typedef Eigen::Matrix<fv, 1, Eigen::Dynamic> RVfv;
typedef Eigen::Matrix<fv, Eigen::Dynamic, Eigen::Dynamic> Mfv;

typedef Eigen::Matrix<ffv, Eigen::Dynamic, 1> Vffv;
typedef Eigen::Matrix<ffv, 1, Eigen::Dynamic> RVffv;
typedef Eigen::Matrix<ffv, Eigen::Dynamic, Eigen::Dynamic> Mffv;

/**
 * Generate a new random variable, cast it to type T1, and store it in z.
 *
 * n1 and n0 are ignored in this specialization.
 *
 * @tparam T1 Scalar type
 * @param n1 Ignored
 * @param n0 Ignored
 * @param z Output variable
 */
template <typename T1>
void build(int n1, int n0, T1& z) {
  std::random_device rd;
  std::mt19937 mt(rd());
  z = T1(rd());
}

/**
 * Generate new random numbers to fill z (casting them to type T1).
 *
 * If z is a column vector, first resize to have n1 rows
 *
 * If z is a row vector, first resize to have n0 columns
 *
 * If z is a matrix, first resize to have n1 rows and n0 columns
 *
 * @tparam T1 Scalar type
 * @param n1 Rows of z (ignored if R == 1)
 * @param n0 Columns of z (ignored if C == 1)
 * @param z Output variable
 */
template <typename T1, int R, int C>
void build(int n1, int n0, Eigen::Matrix<T1, R, C>& z) {
  z = Eigen::Matrix<T1, R, C>(R == 1 ? 1 : n1, C == 1 ? 1 : n0);

  for (int i = 0; i < z.rows(); i++) {
    for (int j = 0; j < z.cols(); j++) {
      build(n1, n0, z(i, j));
    }
  }
}

/**
 * Make a new std::vector of length n2 and recursively call build on
 * each element of it.
 *
 * @tparam T1 Element type
 * @param n2 Dimension of vector
 * @param n1 First dimension of last element
 * @param n0 Second dimension of last element
 * @param z Output variable
 */
template <typename T1>
void build(int n2, int n1, int n0, std::vector<T1>& z) {
  z.resize(n2);
  for (int i = 0; i < n2; i++) {
    build(n1, n0, z[i]);
  }
}

/**
 * Make a new std::vector of length n3 and recursively call build on
 * each element of it.
 *
 * @tparam T1 Element type
 * @param n3 Dimension of std::vector
 * @param n2 Dimension of child std::vector
 * @param n1 First dimension of last element
 * @param n0 Second dimension of last element
 * @param z Output variable
 */
template <typename T1>
void build(int n3, int n2, int n1, int n0, std::vector<T1>& z) {
  z.resize(n3);
  for (int i = 0; i < n3; i++) {
    build(n2, n1, n0, z[i]);
  }
}

/**
 * Get value of variable
 *
 * @param z1 Argument
 */
double get_value(const double& z1) { return z1; }

/**
 * Get value of var.
 *
 * @param z1 Argument
 */
double get_value(const var& z1) { return z1.val(); }

/**
 * Get value of fd
 *
 * @param z1 Argument
 */
double get_value(const fd& z1) { return z1.val(); }

/**
 * Get value of fvar<var>
 *
 * @param z1 Argument
 */
double get_value(const fv& z1) { return z1.val().val(); }

/**
 * Get value of fvar<fvar<double> >
 *
 * @param z1 Argument
 */
double get_value(const ffd& z1) { return z1.val().val(); }

/**
 * Get value of fvar<fvar<var> >
 *
 * @param z1 Argument
 */
double get_value(const ffv& z1) { return z1.val().val().val(); }

/**
 * Check if variables are equal via floating point macro
 *
 * @tparam T1 Type of first argument
 * @tparam T2 Type of second argument
 * @param z1 First argument
 * @param z2 Second argument
 */
template <typename T1, typename T2>
void check_eq(const T1& z1, const T2& z2) {
  EXPECT_FLOAT_EQ(get_value(z1), get_value(z2));
}

/**
 * Check if two integers are equal
 *
 * @param z1 First integer
 * @param z2 Second integer
 */
template <>
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
template <typename T1, typename T2, int R, int C>
void check_eq(const Eigen::Matrix<T1, R, C>& z1,
              const Eigen::Matrix<T2, R, C>& z2) {
  EXPECT_EQ(z1.rows(), z2.rows());
  EXPECT_EQ(z1.cols(), z2.cols());

  for (int i = 0; i < z1.rows(); i++)
    for (int j = 0; j < z1.cols(); j++)
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
template <typename T1, typename T2>
void check_eq(const std::vector<T1>& z1, const std::vector<T2>& z2) {
  EXPECT_EQ(z1.size(), z2.size());
  for (size_t i = 0; i < z1.size(); i++)
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
template <typename T1, typename T2, typename T3>
void checkv() {
  std::vector<T1> x;
  std::vector<T2> y;
  std::vector<T3> result;

  std::random_device rd;
  std::mt19937 mt(rd());

  int r1 = rd() % 5, r2 = rd() % 5 + 1, r3 = rd() % 5 + 1, r4 = rd() % 5;

  build(r1, r2, r3, x);
  build(r4, r2, r3, y);
  EXPECT_NO_THROW(result = stan::math::append_array(x, y));
  EXPECT_EQ(x.size() + y.size(), result.size());
  for (size_t i = 0; i < x.size(); i++)
    check_eq(result[i], x[i]);
  for (size_t i = 0; i < y.size(); i++)
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
template <typename T1, typename T2, typename T3>
void checkvv() {
  std::vector<std::vector<T1> > x;
  std::vector<std::vector<T2> > y;
  std::vector<std::vector<T3> > result;

  std::random_device rd;
  std::mt19937 mt(rd());

  int r1 = rd() % 5, r2 = rd() % 5 + 1, r3 = rd() % 5 + 1, r4 = rd() % 5 + 1,
      r5 = rd() % 5;

  build(r1, r2, r3, r4, x);
  build(r5, r2, r3, r4, y);
  EXPECT_NO_THROW(result = stan::math::append_array(x, y));
  EXPECT_EQ(x.size() + y.size(), result.size());
  for (size_t i = 0; i < x.size(); i++)
    check_eq(result[i], x[i]);
  for (size_t i = 0; i < y.size(); i++)
    check_eq(result[x.size() + i], y[i]);
}

/**
 * Call checkv and checkvv<T1, T2, T3>
 *
 * @tparam T1 Element type of first container
 * @tparam T2 Element type of second container
 * @tparam T3 Element type of third container
 */
template <typename T1, typename T2, typename T3>
void check() {
  // repeat the checks a few times since they're random
  for (int i = 0; i < 3; i++) {
    checkv<T1, T2, T3>();
    checkvv<T1, T2, T3>();
  }
}

TEST(MathFunctions, append_array_prim) {
  check<int, int, int>();
  check<double, double, double>();
  check<double, int, double>();
  check<int, double, double>();
  check<V, V, V>();
  check<RV, RV, RV>();
  check<M, M, M>();
}

TEST(MathFunctions, append_array_rev) {
  check<int, var, var>();
  check<var, int, var>();
  check<double, var, var>();
  check<var, double, var>();
  check<var, var, var>();
  check<V, Vv, Vv>();
  check<Vv, V, Vv>();
  check<Vv, Vv, Vv>();
  check<RV, RVv, RVv>();
  check<RVv, RV, RVv>();
  check<RVv, RVv, RVv>();
  check<M, Mv, Mv>();
  check<Mv, M, Mv>();
  check<Mv, Mv, Mv>();
}

TEST(MathFunctions, append_array_fwd) {
  check<int, fd, fd>();
  check<fd, int, fd>();
  check<double, fd, fd>();
  check<fd, double, fd>();
  check<fd, fd, fd>();
  check<V, Vfd, Vfd>();
  check<Vfd, V, Vfd>();
  check<Vfd, Vfd, Vfd>();
  check<RV, RVfd, RVfd>();
  check<RVfd, RV, RVfd>();
  check<RVfd, RVfd, RVfd>();
  check<M, Mfd, Mfd>();
  check<Mfd, M, Mfd>();
  check<Mfd, Mfd, Mfd>();

  check<int, ffd, ffd>();
  check<ffd, int, ffd>();
  check<double, ffd, ffd>();
  check<ffd, double, ffd>();
  check<ffd, ffd, ffd>();
  check<V, Vffd, Vffd>();
  check<Vffd, V, Vffd>();
  check<Vffd, Vffd, Vffd>();
  check<RV, RVffd, RVffd>();
  check<RVffd, RV, RVffd>();
  check<RVffd, RVffd, RVffd>();
  check<M, Mffd, Mffd>();
  check<Mffd, M, Mffd>();
  check<Mffd, Mffd, Mffd>();
}

TEST(MathFunctions, append_array_mix) {
  check<int, fv, fv>();
  check<fv, int, fv>();
  check<double, fv, fv>();
  check<fv, double, fv>();
  check<fv, fv, fv>();
  check<V, Vfv, Vfv>();
  check<Vfv, V, Vfv>();
  check<Vfv, Vfv, Vfv>();
  check<RV, RVfv, RVfv>();
  check<RVfv, RV, RVfv>();
  check<RVfv, RVfv, RVfv>();
  check<M, Mfv, Mfv>();
  check<Mfv, M, Mfv>();
  check<Mfv, Mfv, Mfv>();

  check<int, ffv, ffv>();
  check<ffv, int, ffv>();
  check<double, ffv, ffv>();
  check<ffv, double, ffv>();
  check<ffv, ffv, ffv>();
  check<V, Vffv, Vffv>();
  check<Vffv, V, Vffv>();
  check<Vffv, Vffv, Vffv>();
  check<RV, RVffv, RVffv>();
  check<RVffv, RV, RVffv>();
  check<RVffv, RVffv, RVffv>();
  check<M, Mffv, Mffv>();
  check<Mffv, M, Mffv>();
  check<Mffv, Mffv, Mffv>();
}
