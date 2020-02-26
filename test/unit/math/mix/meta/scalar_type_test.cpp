#include <stan/math/mix.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <vector>

template <typename R, typename T>
void expect_base() {
  using stan::scalar_type_t;
  static_assert(std::is_same<R, scalar_type_t<T>>::value, "NOT SAME");
  static_assert(std::is_same<R, scalar_type_t<T&>>::value, "NOT SAME");
  test::expect_same_type<R, scalar_type_t<T>>();
  test::expect_same_type<R, scalar_type_t<T>>();
  test::expect_same_type<R, scalar_type_t<T&>>();
  test::expect_same_type<R, scalar_type_t<T&>>();
  test::expect_same_type<R, scalar_type_t<const T&>>();
  test::expect_same_type<R, scalar_type_t<const T&>>();
  test::expect_same_type<R, scalar_type_t<const T>>();
  test::expect_same_type<R, scalar_type_t<const T>>();
}

template <typename T>
void test_base() {
  using Eigen::Matrix;
  using std::complex;
  using std::vector;
  expect_base<T, T>();
  expect_base<T, vector<T>>();
  // std::vector value type has to be copyable _and_ movable, const not movable
  //  expect_base<T, vector<const T>>();

  /**
   *
 (Steve): Idt we have to worry about std vector ref either either
    Unless we start using std::reference_wrapper<T>.

 test/unit/math/mix/meta/scalar_type_test.cpp:57:18:   required from here
 /usr/include/c++/8/ext/new_allocator.h:63:26: error: forming pointer to
 reference type ‘double&’ typedef _Tp*       pointer;
                          ^~~~~~~
 /usr/include/c++/8/ext/new_allocator.h:64:26: error: forming pointer to
 reference type ‘double&’ typedef const _Tp* const_pointer;
                          ^~~~~~~~~~~~~
   */
  // expect_base<T, vector<T&>>();
  expect_base<T, vector<vector<T>>>();
  expect_base<T, Matrix<T, -1, -1>>();
  expect_base<T, Matrix<T, 1, -1>>();
  expect_base<T, Matrix<T, -1, 1>>();
  expect_base<T, vector<Matrix<T, -1, -1>>>();
  expect_base<T, vector<Matrix<T, 1, -1>>>();
  expect_base<T, vector<Matrix<T, -1, 1>>>();

  expect_base<complex<T>, complex<T>>();
}

TEST(MathMetaPrim, scalarType) {
  using stan::base_type;
  using stan::math::fvar;
  using stan::math::var;
  using d_t = double;
  using v_t = var;
  using fd_t = fvar<double>;
  using ffd_t = fvar<fd_t>;
  using fv_t = fvar<var>;
  using ffv_t = fvar<fv_t>;
  test::expect_same_type<int, stan::scalar_type<int>::type>();
  test::expect_same_type<int, stan::scalar_type<std::vector<int>>::type>();

  test_base<d_t>();
  test_base<d_t const*>();
  test_base<v_t>();
  test_base<v_t const*>();
  test_base<fd_t>();
  test_base<ffd_t>();
  test_base<fv_t>();
  test_base<ffv_t>();
}

TEST(MathMetaPrim, ScalarTypeArrayConstConst) {
  using stan::scalar_type;
  using std::vector;
  // These ones work because pointers are always copied
  test::expect_same_type<double const*,
                         scalar_type<const vector<double const*>>::type>();
  test::expect_same_type<int const*,
                         scalar_type<const vector<int const*>>::type>();
  test::expect_same_type<
      double const*, scalar_type<const vector<vector<double const*>>>::type>();
}
