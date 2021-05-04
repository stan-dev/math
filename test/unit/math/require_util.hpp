#ifndef TEST_UNIT_MATH_REQUIRE_UTIL_HPP
#define TEST_UNIT_MATH_REQUIRE_UTIL_HPP

#include <gtest/gtest.h>
#include <type_traits>
#include <string>

namespace stan {
namespace test {
namespace internal {
/** See the link below for the definition of void_t
 * https://en.cppreference.com/w/cpp/types/void_t
 */
template <typename... Ts>
struct make_void {
  typedef void type;
};
template <typename... Ts>
using void_t = typename make_void<Ts...>::type;

/** See the link below for the definition of conjunction
 * https://en.cppreference.com/w/cpp/types/conjunction
 * This is the same as stan::math::conjunction but pulled in here
 *  so we do not need to pull in all of meta.hpp to get it.
 */
template <class...>
struct conjunction : std::true_type {};
template <class B1>
struct conjunction<B1> : B1 {};
template <class B1, class... Bn>
struct conjunction<B1, Bn...>
    : std::conditional_t<bool(B1::value), conjunction<Bn...>, B1> {};

template <class...>
struct disjunction : std::true_type {};
template <class B1>
struct disjunction<B1> : B1 {};
template <class B1, class... Bn>
struct disjunction<B1, Bn...>
    : std::conditional_t<bool(B1::value), B1, disjunction<Bn...>> {};
}  // namespace internal

/**
 * Used for checking an enable_if alias, when the enable_if goes off
 * the struct holds a /c value of true, else false.
 * See the stackoverflow post here for more info.
 * https://stackoverflow.com/q/58189867/2269255
 * @tparam Check The enable_if alias.
 * @tparam T1 In this base impl, the type that failed the /c Check.
 */
template <template <class...> class Check, typename T1, typename = void>
struct unary_require_tester_impl : std::false_type {};

/**
 * Used for checking an enable_if alias, when the enable_if goes off
 * the struct holds a /c value of true, else false.
 * @tparam Check The enable_if alias.
 * @tparam T1 In this base impl, the type that passed the /c Check
 */
template <template <class...> class Check, typename T1>
struct unary_require_tester_impl<Check, T1, internal::void_t<Check<T1>>>
    : std::true_type {};

/**
 * Used for checking an enable_if alias, when the enable_if goes off
 * the struct holds a /c value of true, else false.
 * @tparam Check The enable_if alias.
 * @tparam Types Parameter pack of types that will be checked individually.
 * @note Because of /c conjunction, the member /c value will only be true
 *  if all of the types pass.
 */
template <template <class...> class Check, typename... Types>
struct unary_require_tester
    : internal::conjunction<unary_require_tester_impl<Check, Types>...> {};

// For the end of recursion
template <template <class...> class Check>
struct unary_require_tester<Check, void> : std::true_type {};

/**
 * Used for checking an enable_if alias, when the enable_if goes off
 * the struct holds a /c value of true, else false.
 * @tparam Check The enable_if alias.
 * @tparam Types Parameter pack of types that will be checked individually.
 * @note Because of /c disjunction, the member /c value will only be true
 *  if at least one of the types passes.
 */
template <template <class...> class Check, typename... Types>
struct unary_not_require_tester
    : internal::disjunction<unary_require_tester_impl<Check, Types>...> {};

// For the end of recursion
template <template <class...> class Check>
struct unary_not_require_tester<Check, void> : std::true_type {};

/**
 * This is the same as above, but instead all types are checked at once
 * since this is used for the require_any/all type checks
 */
template <template <class...> class Check, typename Enable, typename... Types>
struct variadic_require_impl : std::false_type {};

template <template <class...> class Check, typename... Types>
struct variadic_require_impl<Check, internal::void_t<Check<Types...>>, Types...>
    : std::true_type {};

template <template <class...> class Check, typename... Types>
struct require_variadic_checker : variadic_require_impl<Check, void, Types...> {
};

// All the checks below are the same for each autodiff alias
template <template <class...> class checker, typename... Types>
struct require_scal_checker {
  struct dummy {};
  static void unary() {
    EXPECT_FALSE((unary_require_tester<checker, dummy>::value));
    EXPECT_TRUE((unary_require_tester<checker, Types...>::value));
  }
  static void not_unary() {
    EXPECT_TRUE((unary_not_require_tester<checker, dummy>::value));
    EXPECT_FALSE((unary_not_require_tester<checker, Types...>::value));
  }
  static void all() {
    EXPECT_FALSE((require_variadic_checker<checker, dummy>::value));
    EXPECT_FALSE((require_variadic_checker<checker, dummy, Types...>::value));
    EXPECT_FALSE((require_variadic_checker<checker, Types..., dummy>::value));
    EXPECT_TRUE((require_variadic_checker<checker, Types...>::value));
  }
  static void all_not() {
    EXPECT_TRUE((require_variadic_checker<checker, dummy>::value));
    EXPECT_FALSE((require_variadic_checker<checker, dummy, Types...>::value));
    EXPECT_FALSE((require_variadic_checker<checker, Types..., dummy>::value));
    EXPECT_FALSE((require_variadic_checker<checker, Types...>::value));
  }
  static void any() {
    EXPECT_FALSE((require_variadic_checker<checker, dummy>::value));
    EXPECT_TRUE((require_variadic_checker<checker, dummy, Types...>::value));
    EXPECT_TRUE((require_variadic_checker<checker, Types..., dummy>::value));
    EXPECT_TRUE((require_variadic_checker<checker, Types...>::value));
  }
  static void any_not() {
    EXPECT_TRUE((require_variadic_checker<checker, dummy>::value));
    EXPECT_TRUE((require_variadic_checker<checker, dummy, Types...>::value));
    EXPECT_TRUE((require_variadic_checker<checker, Types..., dummy>::value));
    EXPECT_FALSE((require_variadic_checker<checker, Types...>::value));
  }
};

// Base impl to check container and underlying type
template <template <template <class...> class TypeCheck, class...> class Check,
          template <class...> class TypeCheck, typename T1, typename = void>
struct unary_container_require_tester_impl : std::false_type {};

template <template <template <class...> class TypeCheck, class...> class Check,
          template <class...> class TypeCheck, typename T1>
struct unary_container_require_tester_impl<
    Check, TypeCheck, T1, stan::test::internal::void_t<Check<TypeCheck, T1>>>
    : std::true_type {};

template <template <template <class...> class TypeCheck, class...> class Check,
          template <class...> class TypeCheck, typename... Types>
struct unary_container_require_tester
    : stan::test::internal::conjunction<
          unary_container_require_tester_impl<Check, TypeCheck, Types>...> {};

// For the end of recursion
template <template <template <class...> class TypeCheck, class...> class Check,
          template <class...> class TypeCheck>
struct unary_container_require_tester<Check, TypeCheck, void> : std::true_type {
};

/**
 * This is the same as above, but instead all types are checked at once
 * since this is used for the require_any/all type checks
 */
template <template <template <class...> class TypeCheck, class...> class Check,
          template <class...> class TypeCheck, typename Enable,
          typename... Types>
struct variadic_container_require_impl : std::false_type {};

template <template <template <class...> class TypeCheck, class...> class Check,
          template <class...> class TypeCheck, typename... Types>
struct variadic_container_require_impl<
    Check, TypeCheck, internal::void_t<Check<TypeCheck, Types...>>, Types...>
    : std::true_type {};

template <template <template <class...> class TypeCheck, class...> class Check,
          template <class...> class TypeCheck, typename... Types>
struct variadic_container_require_tester
    : variadic_container_require_impl<Check, TypeCheck, void, Types...> {};

// All the checks below are the same for each container alias
template <template <template <class...> class TypeCheck, class...>
          class checker,
          template <class...> class Container>
struct require_container_checker {
  template <template <class...> class ContainerCheck>
  static void unary() {
    EXPECT_FALSE((unary_container_require_tester<checker, ContainerCheck,
                                                 double>::value));
    EXPECT_TRUE((unary_container_require_tester<checker, ContainerCheck,
                                                Container<double>>::value));
    EXPECT_TRUE((unary_container_require_tester<checker, ContainerCheck,
                                                Container<float>>::value));
  }

  template <template <class...> class ContainerCheck>
  static void not_unary() {
    EXPECT_TRUE((unary_container_require_tester<checker, ContainerCheck,
                                                double>::value));
    EXPECT_FALSE((unary_container_require_tester<checker, ContainerCheck,
                                                 Container<double>>::value));
    EXPECT_FALSE((unary_container_require_tester<checker, ContainerCheck,
                                                 Container<float>>::value));
  }

  template <template <class...> class ContainerCheck>
  static void all() {
    EXPECT_TRUE((variadic_container_require_tester<checker, ContainerCheck,
                                                   Container<double>,
                                                   Container<double>>::value));
    EXPECT_TRUE((variadic_container_require_tester<checker, ContainerCheck,
                                                   Container<float>,
                                                   Container<double>>::value));
    EXPECT_FALSE(
        (variadic_container_require_tester<checker, ContainerCheck, double,
                                           Container<double>>::value));
    EXPECT_FALSE(
        (variadic_container_require_tester<checker, ContainerCheck,
                                           Container<double>, double>::value));
    EXPECT_FALSE((variadic_container_require_tester<checker, ContainerCheck,
                                                    int, std::string>::value));
    EXPECT_FALSE((variadic_container_require_tester<checker, ContainerCheck,
                                                    double, double>::value));
  }

  template <template <class...> class ContainerCheck>
  static void all_not() {
    EXPECT_FALSE((variadic_container_require_tester<checker, ContainerCheck,
                                                    Container<double>,
                                                    Container<double>>::value));
    EXPECT_FALSE(
        (variadic_container_require_tester<checker, ContainerCheck, double,
                                           Container<double>>::value));
    EXPECT_FALSE(
        (variadic_container_require_tester<checker, ContainerCheck,
                                           Container<double>, double>::value));
    EXPECT_FALSE((variadic_container_require_tester<checker, ContainerCheck,
                                                    Container<std::string>,
                                                    Container<double>>::value));
    EXPECT_TRUE((variadic_container_require_tester<checker, ContainerCheck, int,
                                                   float>::value));
    EXPECT_TRUE((variadic_container_require_tester<checker, ContainerCheck,
                                                   double, double>::value));
  }

  template <template <class...> class ContainerCheck>
  static void any() {
    EXPECT_TRUE((variadic_container_require_tester<checker, ContainerCheck,
                                                   Container<double>,
                                                   Container<double>>::value));
    EXPECT_TRUE(
        (variadic_container_require_tester<checker, ContainerCheck, double,
                                           Container<double>>::value));
    EXPECT_TRUE(
        (variadic_container_require_tester<checker, ContainerCheck,
                                           Container<double>, double>::value));
    EXPECT_TRUE((variadic_container_require_tester<checker, ContainerCheck,
                                                   Container<float>,
                                                   Container<double>>::value));
    EXPECT_FALSE((variadic_container_require_tester<checker, ContainerCheck,
                                                    int, std::string>::value));
    EXPECT_FALSE((variadic_container_require_tester<checker, ContainerCheck,
                                                    double, double>::value));
  }

  template <template <class...> class ContainerCheck>
  static void any_not() {
    EXPECT_FALSE((variadic_container_require_tester<checker, ContainerCheck,
                                                    Container<double>,
                                                    Container<double>>::value));
    EXPECT_FALSE((variadic_container_require_tester<checker, ContainerCheck,
                                                    Container<float>,
                                                    Container<double>>::value));
    EXPECT_TRUE(
        (variadic_container_require_tester<checker, ContainerCheck, double,
                                           Container<double>>::value));
    EXPECT_TRUE(
        (variadic_container_require_tester<checker, ContainerCheck,
                                           Container<double>, double>::value));
    EXPECT_TRUE((variadic_container_require_tester<checker, ContainerCheck, int,
                                                   std::string>::value));
    EXPECT_TRUE((variadic_container_require_tester<checker, ContainerCheck,
                                                   double, double>::value));
  }
};

}  // namespace test
}  // namespace stan

#endif
