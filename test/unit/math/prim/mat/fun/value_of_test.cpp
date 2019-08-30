#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <type_traits>

#include <type_traits>
#include <typeinfo>
#ifndef _MSC_VER
#   include <cxxabi.h>
#endif
#include <memory>
#include <string>
#include <cstdlib>

template <class T>
std::string
type_name()
{
    typedef typename std::remove_reference<T>::type TR;
    std::unique_ptr<char, void(*)(void*)> own
           (
#ifndef _MSC_VER
                abi::__cxa_demangle(typeid(TR).name(), nullptr,
                                           nullptr, nullptr),
#else
                nullptr,
#endif
                std::free
           );
    std::string r = own != nullptr ? own.get() : typeid(TR).name();
    if (std::is_const<TR>::value)
        r += " const";
    if (std::is_volatile<TR>::value)
        r += " volatile";
    if (std::is_lvalue_reference<T>::value)
        r += "&";
    else if (std::is_rvalue_reference<T>::value)
        r += "&&";
    return r;
}

TEST(MathMatrix, value_of) {
  using stan::math::value_of;

  Eigen::Matrix<double, 2, 5> a;
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 5; j++)
      a(i, j) = i * 5 + j;
  Eigen::Matrix<double, 5, 1> b;
  for (int i = 0; i < 5; i++)
    for (int j = 0; j < 1; j++)
      b(i, j) = 10 + i * 5 + j;

  Eigen::MatrixXd d_a = value_of(a);
  Eigen::VectorXd d_b = value_of(b);

  for (int i = 0; i < 5; ++i)
    EXPECT_FLOAT_EQ(b(i), d_b(i));

  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 5; ++j)
      EXPECT_FLOAT_EQ(a(i, j), d_a(i, j));
}

TEST(MathFunctions, value_of_return_type_short_circuit_vector_xd) {
  Eigen::Matrix<double, Eigen::Dynamic, 1> a(5);
  EXPECT_FALSE((std::is_same<decltype(stan::math::value_of(a)),
                             Eigen::Matrix<double, Eigen::Dynamic, 1>>::value));
  EXPECT_FALSE(
      (std::is_same<decltype(stan::math::value_of(a)),
                    Eigen::Matrix<double, Eigen::Dynamic, 1>>::value));
std::cout << "decltype(stan::math::value_of(a)) is " << type_name<decltype(stan::math::value_of(a))>() << '\n';                    
  EXPECT_FALSE(
      (std::is_same<decltype(stan::math::value_of(a)),
                    const Eigen::Matrix<double, Eigen::Dynamic, 1>>::value));
}

TEST(MathFunctions, value_of_return_type_short_circuit_row_vector_xd) {
  Eigen::Matrix<double, 1, Eigen::Dynamic> a(5);
  EXPECT_FALSE((std::is_same<decltype(stan::math::value_of(a)),
                             Eigen::Matrix<double, 1, Eigen::Dynamic>>::value));
  EXPECT_FALSE(
      (std::is_same<decltype(stan::math::value_of(a)),
                    Eigen::Matrix<double, 1, Eigen::Dynamic>>::value));
  EXPECT_FALSE(
      (std::is_same<decltype(stan::math::value_of(a)),
                    const Eigen::Matrix<double, 1, Eigen::Dynamic>>::value));
  EXPECT_TRUE(
      (std::is_same<decltype(stan::math::value_of(a)),
                    const Eigen::Matrix<double, 1, Eigen::Dynamic>>::value));
}

TEST(MathFunctions, value_of_return_type_short_circuit_matrix_xd) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> a(5, 4);
  EXPECT_FALSE((std::is_same<
                decltype(stan::math::value_of(a)),
                Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>::value));
  EXPECT_FALSE(
      (std::is_same<
          decltype(stan::math::value_of(a)),
          Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>::value));
  EXPECT_FALSE(
      (std::is_same<
          decltype(stan::math::value_of(a)),
          const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>::value));
  EXPECT_TRUE((std::is_same<decltype(stan::math::value_of(a)),
                            const Eigen::Matrix<double, Eigen::Dynamic,
                                                Eigen::Dynamic>>::value));
}

TEST(MathFunctions, value_of_return_type_short_circuit_static_sized_matrix) {
  Eigen::Matrix<double, 5, 4> a;
  EXPECT_FALSE((std::is_same<decltype(stan::math::value_of(a)),
                             Eigen::Matrix<double, 5, 4>>::value));
  EXPECT_FALSE((std::is_same<decltype(stan::math::value_of(a)),
                             Eigen::Matrix<double, 5, 4>>::value));
  EXPECT_FALSE((std::is_same<decltype(stan::math::value_of(a)),
                             const Eigen::Matrix<double, 5, 4>>::value));
  EXPECT_TRUE((std::is_same<decltype(stan::math::value_of(a)),
                            const Eigen::Matrix<double, 5, 4>>::value));
}
