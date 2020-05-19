#ifndef TEST_UNIT_UTIL_HPP
#define TEST_UNIT_UTIL_HPP

#include <boost/typeof/typeof.hpp>
#include <gtest/gtest.h>
#include <type_traits>
#include <string>

#define EXPECT_MATRIX_FLOAT_EQ(A, B) \
  EXPECT_EQ(A.rows(), B.rows());        \
  EXPECT_EQ(A.cols(), B.cols());        \
  for (int i = 0; i < A.size(); i++)    \
    EXPECT_FLOAT_EQ(A(i), B(i));

#define EXPECT_MATRIX_NEAR(A, B, DELTA) \
  EXPECT_EQ(A.rows(), B.rows());        \
  EXPECT_EQ(A.cols(), B.cols());        \
  for (int i = 0; i < A.size(); i++)    \
    EXPECT_NEAR(A(i), B(i), DELTA);

#define EXPECT_THROW_MSG_WITH_COUNT(expr, T_e, msg, count) \
  EXPECT_THROW(expr, T_e);                                 \
  try {                                                    \
    expr;                                                  \
  } catch (const T_e& e) {                                 \
    EXPECT_EQ(count, count_matches(msg, e.what()))         \
        << "expected message: " << msg << std::endl        \
        << "found message:    " << e.what();               \
  }

#define EXPECT_THROW_MSG(expr, T_e, msg) \
  EXPECT_THROW_MSG_WITH_COUNT(expr, T_e, msg, 1)

int count_matches(const std::string& target, const std::string& s) {
  if (target.size() == 0)
    return -1;  // error
  int count = 0;
  for (size_t pos = 0; (pos = s.find(target, pos)) != std::string::npos;
       pos += target.size())
    ++count;
  return count;
}

namespace test {
template <typename T1, typename T2>
void expect_same_type() {
  bool b = std::is_same<T1, T2>::value;
  EXPECT_TRUE(b);
}
}  // namespace test

#endif
