#ifndef TEST_UNIT_UTIL_HPP
#define TEST_UNIT_UTIL_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/typeof/typeof.hpp>
#include <test/unit/pretty_print_types.hpp>
#include <gtest/gtest.h>
#include <type_traits>
#include <string>

/**
 * Tests for  exact elementwise equality of the input matrices
 * with the EXPECT_EQ macro from GTest.
 *
 * @param A first input matrix to compare
 * @param B second input matrix to compare
 */
#define EXPECT_MATRIX_EQ(A, B)                                        \
  {                                                                   \
    using T_A = std::decay_t<decltype(A)>;                            \
    using T_B = std::decay_t<decltype(B)>;                            \
    const Eigen::Matrix<typename T_A::Scalar, T_A::RowsAtCompileTime, \
                        T_A::ColsAtCompileTime>                       \
        A_eval = A;                                                   \
    const Eigen::Matrix<typename T_B::Scalar, T_B::RowsAtCompileTime, \
                        T_B::ColsAtCompileTime>                       \
        B_eval = B;                                                   \
    EXPECT_EQ(A_eval.rows(), B_eval.rows());                          \
    EXPECT_EQ(A_eval.cols(), B_eval.cols());                          \
    for (int i = 0; i < A_eval.size(); i++)                           \
      EXPECT_EQ(A_eval(i), B_eval(i));                                \
  }

/**
 * Tests for elementwise equality of the input matrices
 * of doubles with the EXPECT_FLOAT_EQ macro from
 * GTest.
 *
 * @param A first input matrix to compare
 * @param B second input matrix to compare
 */
#define EXPECT_MATRIX_FLOAT_EQ(A, B)         \
  {                                          \
    const Eigen::MatrixXd& A_eval = A;       \
    const Eigen::MatrixXd& B_eval = B;       \
    EXPECT_EQ(A_eval.rows(), B_eval.rows()); \
    EXPECT_EQ(A_eval.cols(), B_eval.cols()); \
    for (int i = 0; i < A_eval.size(); i++)  \
      EXPECT_FLOAT_EQ(A_eval(i), B_eval(i)); \
  }

/**
 * Tests for elementwise equality of the input matrices
 * of std::complex<double>s with the EXPECT_FLOAT_EQ macro
 * from GTest.
 *
 * @param A first input matrix to compare
 * @param B second input matrix to compare
 */
#define EXPECT_MATRIX_COMPLEX_FLOAT_EQ(A, B)               \
  {                                                        \
    const Eigen::MatrixXcd& A_eval = A;                    \
    const Eigen::MatrixXcd& B_eval = B;                    \
    EXPECT_EQ(A_eval.rows(), B_eval.rows());               \
    EXPECT_EQ(A_eval.cols(), B_eval.cols());               \
    for (int i = 0; i < A_eval.size(); i++) {              \
      EXPECT_FLOAT_EQ(A_eval(i).real(), B_eval(i).real()); \
      EXPECT_FLOAT_EQ(A_eval(i).imag(), B_eval(i).imag()); \
    }                                                      \
  }

/**
 * Tests for elementwise equality of the input std::vectors
 * of any type with the EXPECT_FLOAT_EQ macro from GTest.
 *
 * @param A first input vector to compare
 * @param B second input vector to compare
 */
#define EXPECT_STD_VECTOR_FLOAT_EQ(A, B) \
  EXPECT_EQ(A.size(), B.size());         \
  for (int i = 0; i < A.size(); ++i)     \
    EXPECT_FLOAT_EQ(A[i], B[i]);

/**
 * Tests for elementwise equality of the input std::vectors
 * of any type with the EXPECT_EQ macro from GTest.
 *
 * @param A first input vector of integers to compare
 * @param B second input vector of integers to compare
 */
#define EXPECT_STD_VECTOR_EQ(A, B)   \
  EXPECT_EQ(A.size(), B.size());     \
  for (int i = 0; i < A.size(); ++i) \
    EXPECT_EQ(A[i], B[i]);

/**
 * Tests if any elementwise difference of the input matrices
 * of doubles is greater than DELTA. This uses the
 * EXPECT_NEAR macro from GTest.
 *
 * @param A first input matrix to compare
 * @param B second input matrix to compare
 * @param DELTA the maximum allowed difference
 */
#define EXPECT_MATRIX_NEAR(A, B, DELTA)                               \
  {                                                                   \
    using T_A = std::decay_t<decltype(A)>;                            \
    using T_B = std::decay_t<decltype(B)>;                            \
    const Eigen::Matrix<typename T_A::Scalar, T_A::RowsAtCompileTime, \
                        T_A::ColsAtCompileTime>                       \
        A_eval = A;                                                   \
    const Eigen::Matrix<typename T_B::Scalar, T_B::RowsAtCompileTime, \
                        T_B::ColsAtCompileTime>                       \
        B_eval = B;                                                   \
    EXPECT_EQ(A_eval.rows(), B_eval.rows());                          \
    EXPECT_EQ(A_eval.cols(), B_eval.cols());                          \
    for (int i = 0; i < A_eval.size(); i++)                           \
      EXPECT_NEAR(A_eval(i), B_eval(i), DELTA);                       \
  }

/**
 * Tests for elementwise equality of the input matrices
 * of std::complex<double>s with the EXPECT_FLOAT_EQ macro
 * from GTest.
 *
 * @param A first input matrix to compare
 * @param B second input matrix to compare
 */
#define EXPECT_MATRIX_COMPLEX_NEAR(A, B, DELTA)               \
  {                                                           \
    const Eigen::MatrixXcd& A_eval = A;                       \
    const Eigen::MatrixXcd& B_eval = B;                       \
    EXPECT_EQ(A_eval.rows(), B_eval.rows());                  \
    EXPECT_EQ(A_eval.cols(), B_eval.cols());                  \
    for (int i = 0; i < A_eval.size(); i++) {                 \
      EXPECT_NEAR(A_eval(i).real(), B_eval(i).real(), DELTA); \
      EXPECT_NEAR(A_eval(i).imag(), B_eval(i).imag(), DELTA); \
    }                                                         \
  }

/**
 * Tests if given types are the same type.
 *
 * @param a first type
 * @param b second type (code for this one can contain commas)
 **/
#define EXPECT_SAME_TYPE(a, ...)                                           \
  EXPECT_TRUE((std::is_same<a, __VA_ARGS__>::value))                       \
      << "Type a is" << stan::math::test::type_name<a>() << ". Type b is " \
      << stan::math::test::type_name<__VA_ARGS__>();

/**
 * Tests if given value is of given type.
 *
 * @param type type
 * @param value value (code for this one can contain commas)
 **/
#define EXPECT_TYPE(type, ...) EXPECT_SAME_TYPE(type, decltype(__VA_ARGS__))

/**
 * Count the number of times a substring is found in
 * a supplied string.
 *
 * @param target substring to match in s
 * @param s string to match count occurrences
 * @return number of found occurrences of target in s
 */
int count_matches(const std::string& target, const std::string& s) {
  if (target.size() == 0)
    return -1;  // error
  int count = 0;
  for (size_t pos = 0; (pos = s.find(target, pos)) != std::string::npos;
       pos += target.size())
    ++count;
  return count;
}

/**
 * Tests if the expression throws the expected
 * exception with a specific number of occurrences of
 * the expected message in the throw message.
 *
 * @param expr expression to test
 * @param T_e type of exception
 * @param msg expected message
 * @param count the expected number of occurrences of msg
 * in the message of the exception
 */
#define EXPECT_THROW_MSG_WITH_COUNT(expr, T_e, msg, count) \
  EXPECT_THROW(expr, T_e);                                 \
  try {                                                    \
    expr;                                                  \
  } catch (const T_e& e) {                                 \
    EXPECT_EQ(count, count_matches(msg, e.what()))         \
        << "expected message: " << msg << std::endl        \
        << "found message:    " << e.what();               \
  }

/**
 * Tests if the expression throws the expected
 * exception with the expected message.
 *
 * @param expr expression to test
 * @param T_e type of exception
 * @param msg expected message
 */
#define EXPECT_THROW_MSG(expr, T_e, msg) \
  EXPECT_THROW_MSG_WITH_COUNT(expr, T_e, msg, 1)

namespace stan {
namespace test {

auto make_sparse_matrix_random(int rows, int cols) {
  using eigen_triplet = Eigen::Triplet<double>;
  boost::mt19937 gen;
  boost::random::uniform_real_distribution<double> dist(0.0, 1.0);
  std::vector<eigen_triplet> tripletList;
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      auto v_ij = dist(gen);
      if (v_ij < 0.1) {
        tripletList.push_back(eigen_triplet(i, j, v_ij));
      }
    }
  }
  Eigen::SparseMatrix<double> mat(rows, cols);
  mat.setFromTriplets(tripletList.begin(), tripletList.end());
  return mat;
}

}  // namespace test
}  // namespace stan

#endif
