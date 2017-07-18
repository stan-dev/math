#include <stan/math/mix/mat.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <gtest/gtest.h>

using namespace Eigen;
using stan::math::var;
using stan::math::fvar;

TEST(AgradMix, append_array_fvar_var_double) {
  std::vector<double> x(2);
  std::vector<fvar<var> > y(3), result;

  x[0] = 1.0;
  x[1] = 2.0;
  y[0] = 5.0;
  y[0].d_ = 1.5;
  y[1] = 6.0;
  y[1].d_ = -2.5;
  y[2] = 7.0;
  y[2].d_ = -3.5;

  EXPECT_NO_THROW(result = stan::math::append_array(y, x));
  EXPECT_EQ(5, result.size());
  EXPECT_FLOAT_EQ(5.0, result[0].val().val());
  EXPECT_FLOAT_EQ(6.0, result[1].val().val());
  EXPECT_FLOAT_EQ(7.0, result[2].val().val());
  EXPECT_FLOAT_EQ(1.0, result[3].val().val());
  EXPECT_FLOAT_EQ(2.0, result[4].val().val());

  EXPECT_FLOAT_EQ(1.5, result[0].tangent().val());
  EXPECT_FLOAT_EQ(-2.5, result[1].tangent().val());
  EXPECT_FLOAT_EQ(-3.5, result[2].tangent().val());
  EXPECT_FLOAT_EQ(0.0, result[3].tangent().val());
  EXPECT_FLOAT_EQ(0.0, result[4].tangent().val());

  stan::math::set_zero_all_adjoints();
  result[0].val().grad();
  EXPECT_FLOAT_EQ(1.0, y[0].val().adj());

  stan::math::set_zero_all_adjoints();
  result[0].tangent().grad();
  EXPECT_FLOAT_EQ(0.0, y[1].tangent().adj());
}

TEST(AgradMix, append_array_double_fvar_var) {
  std::vector<double> x(2);
  std::vector<fvar<var> > y(3), result;

  x[0] = 1.0;
  x[1] = 2.0;
  y[0] = 5.0;
  y[0].d_ = 1.5;
  y[1] = 6.0;
  y[1].d_ = -2.5;
  y[2] = 7.0;
  y[2].d_ = -3.5;

  EXPECT_NO_THROW(result = stan::math::append_array(x, y));
  EXPECT_EQ(5, result.size());
  EXPECT_FLOAT_EQ(1.0, result[0].val().val());
  EXPECT_FLOAT_EQ(2.0, result[1].val().val());
  EXPECT_FLOAT_EQ(5.0, result[2].val().val());
  EXPECT_FLOAT_EQ(6.0, result[3].val().val());
  EXPECT_FLOAT_EQ(7.0, result[4].val().val());

  EXPECT_FLOAT_EQ(0.0, result[0].tangent().val());
  EXPECT_FLOAT_EQ(0.0, result[1].tangent().val());
  EXPECT_FLOAT_EQ(1.5, result[2].tangent().val());
  EXPECT_FLOAT_EQ(-2.5, result[3].tangent().val());
  EXPECT_FLOAT_EQ(-3.5, result[4].tangent().val());

  stan::math::set_zero_all_adjoints();
  result[3].val().grad();
  EXPECT_FLOAT_EQ(1.0, y[1].val().adj());

  stan::math::set_zero_all_adjoints();
  result[0].tangent().grad();
  EXPECT_FLOAT_EQ(0.0, y[0].tangent().adj());
}

TEST(AgradMix, append_array_fvar_var_fvar_var) {
  std::vector<fvar<var> > x(2), y(3), result;

  x[0] = 1.0;
  x[0].d_ = 2.5;
  x[1] = 2.0;
  x[1].d_ = 3.5;
  y[0] = 5.0;
  y[0].d_ = 1.5;
  y[1] = 6.0;
  y[1].d_ = -2.5;
  y[2] = 7.0;
  y[2].d_ = -3.5;

  EXPECT_NO_THROW(result = stan::math::append_array(x, y));
  EXPECT_EQ(5, result.size());
  EXPECT_FLOAT_EQ(1.0, result[0].val().val());
  EXPECT_FLOAT_EQ(2.0, result[1].val().val());
  EXPECT_FLOAT_EQ(5.0, result[2].val().val());
  EXPECT_FLOAT_EQ(6.0, result[3].val().val());
  EXPECT_FLOAT_EQ(7.0, result[4].val().val());

  EXPECT_FLOAT_EQ(2.5, result[0].tangent().val());
  EXPECT_FLOAT_EQ(3.5, result[1].tangent().val());
  EXPECT_FLOAT_EQ(1.5, result[2].tangent().val());
  EXPECT_FLOAT_EQ(-2.5, result[3].tangent().val());
  EXPECT_FLOAT_EQ(-3.5, result[4].tangent().val());

  stan::math::set_zero_all_adjoints();
  result[2].val().grad();
  EXPECT_FLOAT_EQ(1.0, y[0].val().adj());

  stan::math::set_zero_all_adjoints();
  result[1].tangent().grad();
  EXPECT_FLOAT_EQ(1.0, x[1].tangent().adj());
}

TEST(AgradMix, append_array_fvar_fvar_var_double) {
  std::vector<double> x(2);
  std::vector<fvar<fvar<var> > > y(3), result;

  x[0] = 1.0;
  x[1] = 2.0;
  y[0] = 5.0;
  y[0].d_ = 1.5;
  y[1] = 6.0;
  y[1].d_ = -2.5;
  y[2] = 7.0;
  y[2].d_ = -3.5;

  EXPECT_NO_THROW(result = stan::math::append_array(y, x));
  EXPECT_EQ(5, result.size());
  EXPECT_FLOAT_EQ(5.0, result[0].val().val().val());
  EXPECT_FLOAT_EQ(6.0, result[1].val().val().val());
  EXPECT_FLOAT_EQ(7.0, result[2].val().val().val());
  EXPECT_FLOAT_EQ(1.0, result[3].val().val().val());
  EXPECT_FLOAT_EQ(2.0, result[4].val().val().val());

  EXPECT_FLOAT_EQ(1.5, result[0].tangent().val().val());
  EXPECT_FLOAT_EQ(-2.5, result[1].tangent().val().val());
  EXPECT_FLOAT_EQ(-3.5, result[2].tangent().val().val());
  EXPECT_FLOAT_EQ(0.0, result[3].tangent().val().val());
  EXPECT_FLOAT_EQ(0.0, result[4].tangent().val().val());

  for (size_t i = 0; i < result.size(); i++) {
    EXPECT_FLOAT_EQ(0.0, result[i].val().tangent().val());
    EXPECT_FLOAT_EQ(0.0, result[i].tangent().tangent().val());
  }

  stan::math::set_zero_all_adjoints();
  result[0].val().val().grad();
  EXPECT_FLOAT_EQ(1.0, y[0].val().val().adj());

  stan::math::set_zero_all_adjoints();
  result[0].val().tangent().grad();
  EXPECT_FLOAT_EQ(1.0, y[0].val().tangent().adj());
}

TEST(AgradMix, append_array_double_fvar_fvar_var) {
  std::vector<double> x(2);
  std::vector<fvar<fvar<var> > > y(3), result;

  x[0] = 1.0;
  x[1] = 2.0;
  y[0] = 5.0;
  y[0].d_ = 1.5;
  y[1] = 6.0;
  y[1].d_ = -2.5;
  y[2] = 7.0;
  y[2].d_ = -3.5;

  EXPECT_NO_THROW(result = stan::math::append_array(x, y));
  EXPECT_EQ(5, result.size());
  EXPECT_FLOAT_EQ(1.0, result[0].val().val().val());
  EXPECT_FLOAT_EQ(2.0, result[1].val().val().val());
  EXPECT_FLOAT_EQ(5.0, result[2].val().val().val());
  EXPECT_FLOAT_EQ(6.0, result[3].val().val().val());
  EXPECT_FLOAT_EQ(7.0, result[4].val().val().val());

  EXPECT_FLOAT_EQ(0.0, result[0].tangent().val().val());
  EXPECT_FLOAT_EQ(0.0, result[1].tangent().val().val());
  EXPECT_FLOAT_EQ(1.5, result[2].tangent().val().val());
  EXPECT_FLOAT_EQ(-2.5, result[3].tangent().val().val());
  EXPECT_FLOAT_EQ(-3.5, result[4].tangent().val().val());

  for (size_t i = 0; i < result.size(); i++) {
    EXPECT_FLOAT_EQ(0.0, result[i].val().tangent().val());
    EXPECT_FLOAT_EQ(0.0, result[i].tangent().tangent().val());
  }

  stan::math::set_zero_all_adjoints();
  result[2].val().val().grad();
  EXPECT_FLOAT_EQ(1.0, y[0].val().val().adj());

  stan::math::set_zero_all_adjoints();
  result[2].val().tangent().grad();
  EXPECT_FLOAT_EQ(1.0, y[0].val().tangent().adj());
}

TEST(AgradMix, append_array_fvar_fvar_var_fvar_fvar_var) {
  std::vector<fvar<fvar<var> > > x(2), y(3), result;

  x[0] = 1.0;
  x[0].d_ = 2.5;
  x[1] = 2.0;
  x[1].d_ = 3.5;
  y[0] = 5.0;
  y[0].d_ = 1.5;
  y[1] = 6.0;
  y[1].d_ = -2.5;
  y[2] = 7.0;
  y[2].d_ = -3.5;

  EXPECT_NO_THROW(result = stan::math::append_array(x, y));
  EXPECT_EQ(5, result.size());
  EXPECT_FLOAT_EQ(1.0, result[0].val().val().val());
  EXPECT_FLOAT_EQ(2.0, result[1].val().val().val());
  EXPECT_FLOAT_EQ(5.0, result[2].val().val().val());
  EXPECT_FLOAT_EQ(6.0, result[3].val().val().val());
  EXPECT_FLOAT_EQ(7.0, result[4].val().val().val());

  EXPECT_FLOAT_EQ(2.5, result[0].tangent().val().val());
  EXPECT_FLOAT_EQ(3.5, result[1].tangent().val().val());
  EXPECT_FLOAT_EQ(1.5, result[2].tangent().val().val());
  EXPECT_FLOAT_EQ(-2.5, result[3].tangent().val().val());
  EXPECT_FLOAT_EQ(-3.5, result[4].tangent().val().val());

  for (size_t i = 0; i < result.size(); i++) {
    EXPECT_FLOAT_EQ(0.0, result[i].tangent().tangent().val());
    EXPECT_FLOAT_EQ(0.0, result[i].val().tangent().val());
  }

  stan::math::set_zero_all_adjoints();
  result[0].tangent().val().grad();
  EXPECT_FLOAT_EQ(1.0, x[0].tangent().val().adj());

  stan::math::set_zero_all_adjoints();
  result[2].tangent().tangent().grad();
  EXPECT_FLOAT_EQ(1.0, y[0].tangent().tangent().adj());
}

TEST(AgradMix, append_array_matrix_double_matrix_fvar_var) {
  std::vector<Matrix<double, Dynamic, Dynamic> > x;
  std::vector<Matrix<fvar<var>, Dynamic, Dynamic> > y, result;

  for (int i = 0; i < 3; i++)
    x.push_back(Matrix<double, Dynamic, Dynamic>::Zero(3, 3));
  for (int i = 0; i < 2; i++)
    y.push_back(Matrix<fvar<var>, Dynamic, Dynamic>::Zero(3, 3));

  x[0](0, 0) = 1.0;
  y[1](2, 1) = 2.0;
  y[1](2, 1).d_ = 3.0;

  EXPECT_NO_THROW(result = stan::math::append_array(x, y));
  EXPECT_EQ(5, result.size());
  EXPECT_FLOAT_EQ(1.0, result[0](0, 0).val().val());
  EXPECT_FLOAT_EQ(2.0, result[4](2, 1).val().val());
  EXPECT_FLOAT_EQ(3.0, result[4](2, 1).tangent().val());
  EXPECT_FLOAT_EQ(0.0, result[4](2, 2).val().val());

  stan::math::set_zero_all_adjoints();
  result[3](2, 1).val().grad();
  EXPECT_FLOAT_EQ(1.0, y[0](2, 1).val().adj());

  stan::math::set_zero_all_adjoints();
  result[3](2, 1).tangent().grad();
  EXPECT_FLOAT_EQ(1.0, y[0](2, 1).tangent().adj());

  for (int i = 0; i < 2; i++)
    y[i] = Matrix<fvar<var>, Dynamic, Dynamic>::Zero(2, 1);

  EXPECT_THROW(result = stan::math::append_array(x, y), std::invalid_argument);
}

TEST(AgradMix, append_array_matrix_double_matrix_fvar_fvar_var) {
  std::vector<Matrix<double, Dynamic, Dynamic> > x;
  std::vector<Matrix<fvar<fvar<var> >, Dynamic, Dynamic> > y, result;

  for (int i = 0; i < 3; i++)
    x.push_back(Matrix<double, Dynamic, Dynamic>::Zero(3, 3));
  for (int i = 0; i < 2; i++)
    y.push_back(Matrix<fvar<fvar<var> >, Dynamic, Dynamic>::Zero(3, 3));

  x[0](0, 0) = 1.0;
  y[1](2, 1).val_ = 2.0;
  y[1](2, 1).val_.d_ = 3.0;
  y[1](2, 1).d_.val_ = 4.0;
  y[1](2, 1).d_.d_ = 5.0;

  EXPECT_NO_THROW(result = stan::math::append_array(x, y));
  EXPECT_EQ(5, result.size());
  EXPECT_FLOAT_EQ(1.0, result[0](0, 0).val().val().val());
  EXPECT_FLOAT_EQ(2.0, result[4](2, 1).val().val().val());
  EXPECT_FLOAT_EQ(3.0, result[4](2, 1).val().tangent().val());
  EXPECT_FLOAT_EQ(4.0, result[4](2, 1).tangent().val().val());
  EXPECT_FLOAT_EQ(5.0, result[4](2, 1).tangent().tangent().val());
  EXPECT_FLOAT_EQ(0.0, result[4](2, 2).val().val().val());

  stan::math::set_zero_all_adjoints();
  result[3](2, 1).val().val().grad();
  EXPECT_FLOAT_EQ(1.0, y[0](2, 1).val().val().adj());

  stan::math::set_zero_all_adjoints();
  result[3](2, 1).tangent().val().grad();
  EXPECT_FLOAT_EQ(1.0, y[0](2, 1).tangent().val().adj());

  for (int i = 0; i < 2; i++)
    y[i] = Matrix<fvar<fvar<var> >, Dynamic, Dynamic>::Zero(1, 5);

  EXPECT_THROW(result = stan::math::append_array(x, y), std::invalid_argument);
}

TEST(AgradMix, append_array_matrix_fvar_fvar_var_matrix_fvar_fvar_var) {
  std::vector<Matrix<fvar<fvar<var> >, Dynamic, Dynamic> > x, y, result;

  for (int i = 0; i < 3; i++)
    x.push_back(Matrix<fvar<fvar<var> >, Dynamic, Dynamic>::Zero(3, 3));
  for (int i = 0; i < 2; i++)
    y.push_back(Matrix<fvar<fvar<var> >, Dynamic, Dynamic>::Zero(3, 3));

  x[0](0, 0).val_ = 1.0;
  x[0](0, 0).val_.d_ = 6.0;
  x[0](0, 0).d_.val_ = 7.0;
  x[0](0, 0).d_.d_ = 8.0;
  y[1](2, 1).val_ = 2.0;
  y[1](2, 1).val_.d_ = 3.0;
  y[1](2, 1).d_.val_ = 4.0;
  y[1](2, 1).d_.d_ = 5.0;

  EXPECT_NO_THROW(result = stan::math::append_array(x, y));
  EXPECT_EQ(5, result.size());
  EXPECT_FLOAT_EQ(1.0, result[0](0, 0).val().val().val());
  EXPECT_FLOAT_EQ(6.0, result[0](0, 0).val().tangent().val());
  EXPECT_FLOAT_EQ(7.0, result[0](0, 0).tangent().val().val());
  EXPECT_FLOAT_EQ(8.0, result[0](0, 0).tangent().tangent().val());
  EXPECT_FLOAT_EQ(2.0, result[4](2, 1).val().val().val());
  EXPECT_FLOAT_EQ(3.0, result[4](2, 1).val().tangent().val());
  EXPECT_FLOAT_EQ(4.0, result[4](2, 1).tangent().val().val());
  EXPECT_FLOAT_EQ(5.0, result[4](2, 1).tangent().tangent().val());
  EXPECT_FLOAT_EQ(0.0, result[4](2, 2).val().val().val());

  stan::math::set_zero_all_adjoints();
  result[0](0, 0).val().val().grad();
  EXPECT_FLOAT_EQ(1.0, x[0](0, 0).val().val().adj());

  stan::math::set_zero_all_adjoints();
  result[0](0, 0).val().tangent().grad();
  EXPECT_FLOAT_EQ(1.0, x[0](0, 0).val().tangent().adj());

  stan::math::set_zero_all_adjoints();
  result[3](2, 1).val().val().grad();
  EXPECT_FLOAT_EQ(1.0, y[0](2, 1).val().val().adj());

  stan::math::set_zero_all_adjoints();
  result[3](2, 1).tangent().val().grad();
  EXPECT_FLOAT_EQ(1.0, y[0](2, 1).tangent().val().adj());

  for (int i = 0; i < 2; i++)
    y[i] = Matrix<fvar<fvar<var> >, Dynamic, Dynamic>::Zero(1, 5);

  EXPECT_THROW(result = stan::math::append_array(x, y), std::invalid_argument);
}

typedef fvar<double> fd;
typedef fvar<var> fv;

typedef fvar<fvar<double> > ffd;
typedef fvar<fvar<var> > ffv;

typedef Matrix<double, Dynamic, 1> V;
typedef Matrix<double, 1, Dynamic> RV;
typedef Matrix<double, Dynamic, Dynamic> M;

typedef Matrix<var, Dynamic, 1> Vv;
typedef Matrix<var, 1, Dynamic> RVv;
typedef Matrix<var, Dynamic, Dynamic> Mv;

typedef Matrix<fd, Dynamic, 1> Vfd;
typedef Matrix<fd, 1, Dynamic> RVfd;
typedef Matrix<fd, Dynamic, Dynamic> Mfd;

typedef Matrix<ffd, Dynamic, 1> Vffd;
typedef Matrix<ffd, 1, Dynamic> RVffd;
typedef Matrix<ffd, Dynamic, Dynamic> Mffd;

typedef Matrix<fv, Dynamic, 1> Vfv;
typedef Matrix<fv, 1, Dynamic> RVfv;
typedef Matrix<fv, Dynamic, Dynamic> Mfv;

typedef Matrix<ffv, Dynamic, 1> Vffv;
typedef Matrix<ffv, 1, Dynamic> RVffv;
typedef Matrix<ffv, Dynamic, Dynamic> Mffv;

/**
 * Generate a new random variable, cast it to type T1, and store it in z.
 *
 * n0 is ignored in this specialization.
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
 * Generate new random numbers to fill z (casting them to type T1).
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
  z = Eigen::Matrix<T1, R, C>(R == 1 ? 1 : n0, C == 1 ? 1 : n0);

  for (int i = 0; i < z.rows(); i++) {
    for (int j = 0; j < z.cols(); j++) {
      build(n0, z(i, j));
    }
  }
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
  for (int i = 0; i < n1; i++) {
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
  for (int i = 0; i < n2; i++) {
    build(n1, n0, z[i]);
  }
}

/**
 * Get value of variable.
 *
 * @param z1 Argument
 */
double get_value(const double& z1) {
  return z1;
}

/**
 * Get value of var.
 *
 * @param z1 Argument
 */
double get_value(const var& z1) {
  return z1.val();
}

/**
 * Get value of fd
 *
 * @param z1 Argument
 */
double get_value(const fd& z1) {
  return z1.val();
}

/**
 * Get value of fvar<var>
 *
 * @param z1 Argument
 */
double get_value(const fv& z1) {
  return z1.val().val();
}

/**
 * Get value of fvar<fvar<double> >
 *
 * @param z1 Argument
 */
double get_value(const ffd& z1) {
  return z1.val().val();
}

/**
 * Get value of fvar<fvar<var> >
 *
 * @param z1 Argument
 */
double get_value(const ffv& z1) {
  return z1.val().val().val();
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
  EXPECT_FLOAT_EQ(get_value(z1), get_value(z2));
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
template<typename T1, typename T2>
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
template<typename T1, typename T2, typename T3>
void check() {
  checkv<T1, T2, T3>();
  checkvv<T1, T2, T3>();
}

TEST(MathFunctions, append_array_prim) {
  check<double, double, double>();
  check<double, int, double>();
  check<int, double, double>();
  check<V, V, V>();
  check<RV, RV, RV>();
  check<M, M, M>();
}

TEST(MathFunctions, append_array_rev) {
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
