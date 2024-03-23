#ifndef TEST_UNIT_MATH_REV_UTIL_HPP
#define TEST_UNIT_MATH_REV_UTIL_HPP

#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <vector>

struct AgradRev : public testing::Test {
  void SetUp() {
    // make sure memory's clean before starting each test
    stan::math::recover_memory();
  }
  void TearDown() {
    // make sure memory's clean after each test
    stan::math::recover_memory();
  }
};

namespace stan {
namespace math {
namespace test {
template <bool UseVarMat>
using cond_var_matrix_t = std::conditional_t<
    UseVarMat,
    stan::math::var_value<
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>,
    Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic>>;
template <bool UseVarMat>
using cond_var_vector_t = std::conditional_t<
    UseVarMat, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>>;
template <bool UseVarMat>
using cond_var_row_vector_t = std::conditional_t<
    UseVarMat, stan::math::var_value<Eigen::Matrix<double, 1, Eigen::Dynamic>>,
    Eigen::Matrix<stan::math::var, 1, Eigen::Dynamic>>;

template <bool UseVarMat>
struct var_matrix_types {
  using matrix_v = stan::math::test::cond_var_matrix_t<UseVarMat>;
  using row_vector_v = stan::math::test::cond_var_row_vector_t<UseVarMat>;
  using vector_v = stan::math::test::cond_var_vector_t<UseVarMat>;
};

template <class T>
class VarMatrixTypedTests : public testing::Test {
 public:
  using matrix_v = typename T::matrix_v;
  using row_vector_v = typename T::row_vector_v;
  using vector_v = typename T::vector_v;
  virtual ~VarMatrixTypedTests() { stan::math::recover_memory(); }
};
using VarMatImpls = testing::Types<stan::math::test::var_matrix_types<false>,
                                   stan::math::test::var_matrix_types<true>>;
}  // namespace test
}  // namespace math
}  // namespace stan

namespace test {

void check_varis_on_stack(const stan::math::var& x) {
  EXPECT_TRUE(stan::math::ChainableStack::instance_->memalloc_.in_stack(x.vi_))
      << "not on the stack";
}

void check_varis_on_stack(const std::vector<stan::math::var>& x) {
  for (size_t n = 0; n < x.size(); ++n)
    EXPECT_TRUE(
        stan::math::ChainableStack::instance_->memalloc_.in_stack(x[n].vi_))
        << n << " is not on the stack";
}

template <int R, int C>
void check_varis_on_stack(const Eigen::Matrix<stan::math::var, R, C>& x) {
  for (int j = 0; j < x.cols(); ++j)
    for (int i = 0; i < x.rows(); ++i)
      EXPECT_TRUE(stan::math::ChainableStack::instance_->memalloc_.in_stack(
          x(i, j).vi_))
          << i << ", " << j << " is not on the stack";
}

}  // namespace test
#endif
