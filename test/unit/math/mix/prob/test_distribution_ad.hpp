#ifndef TEST_UNIT_MATH_TEST_DISTRIBUTION_AD_HPP
#define TEST_UNIT_MATH_TEST_DISTRIBUTION_AD_HPP

#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>
namespace stan {
  namespace test {
    template <typename Cast, typename T, stan::require_eigen_t<Cast>* = nullptr>
    Cast do_cast(T&& x) {
        return x;
    }
    template <typename Cast, typename T, stan::require_std_vector_t<Cast>* = nullptr>
    Cast do_cast(T&& x) {
        return Cast(x.data(), x.data() + x.size());
    }
    template <typename Cast, typename T, stan::require_arithmetic_t<Cast>* = nullptr>
    Cast do_cast(T&& x) {
        return x(0);
    }

  using unary_base_set = std::tuple<std::tuple<double>,
    std::tuple<std::vector<double>>,
    std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>>>;

  using binary_base_set = std::tuple<std::tuple<double, double>,
  std::tuple<double, std::vector<double>>,
  std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>>,
  std::tuple<std::vector<double>, double>,
  std::tuple<std::vector<double>, std::vector<double>>,
  std::tuple<std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>>,
  std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double>,
  std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>>,
  std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Matrix<double, Eigen::Dynamic, 1>>>;


  using ternary_base_sets = std::tuple<std::tuple<double, double, double>,
  std::tuple<double, double, std::vector<double>>,
  std::tuple<double, double, Eigen::Matrix<double, Eigen::Dynamic, 1>>,
  std::tuple<double, std::vector<double>, double>,
  std::tuple<double, std::vector<double>, std::vector<double>>,
  std::tuple<double, std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>>,
  std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>, double>,
  std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>>,
  std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Matrix<double, Eigen::Dynamic, 1>>,
  std::tuple<std::vector<double>, double, double>,
  std::tuple<std::vector<double>, double, std::vector<double>>,
  std::tuple<std::vector<double>, double, Eigen::Matrix<double, Eigen::Dynamic, 1>>,
  std::tuple<std::vector<double>, std::vector<double>, double>,
  std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>,
  std::tuple<std::vector<double>, std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>>,
  std::tuple<std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, double>,
  std::tuple<std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>>,
  std::tuple<std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Matrix<double, Eigen::Dynamic, 1>>,
  std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double, double>,
  std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double, std::vector<double>>,
  std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double, Eigen::Matrix<double, Eigen::Dynamic, 1>>,
  std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, double>,
  std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, std::vector<double>>,
  std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>>,
  std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Matrix<double, Eigen::Dynamic, 1>, double>,
  std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>>,
  std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Matrix<double, Eigen::Dynamic, 1>>>;


using quad_base_set = std::tuple<std::tuple<double, double, double, double>,
  std::tuple<double, double, double, std::vector<double>>,
  std::tuple<double, double, double, Eigen::Matrix<double, Eigen::Dynamic, 1>>,
  std::tuple<double, double, std::vector<double>, double>,
  std::tuple<double, double, std::vector<double>, std::vector<double>>,
  std::tuple<double, double, std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>>,
  std::tuple<double, double, Eigen::Matrix<double, Eigen::Dynamic, 1>, double>,
  std::tuple<double, double, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>>,
  std::tuple<double, double, Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Matrix<double, Eigen::Dynamic, 1>>,
  std::tuple<double, std::vector<double>, double, double>,
  std::tuple<double, std::vector<double>, double, std::vector<double>>,
  std::tuple<double, std::vector<double>, double, Eigen::Matrix<double, Eigen::Dynamic, 1>>,
  std::tuple<double, std::vector<double>, std::vector<double>, double>,
  std::tuple<double, std::vector<double>, std::vector<double>, std::vector<double>>,
  std::tuple<double, std::vector<double>, std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>>,
  std::tuple<double, std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, double>,
  std::tuple<double, std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>>,
  std::tuple<double, std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Matrix<double, Eigen::Dynamic, 1>>,
  std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>, double, double>,
  std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>, double, std::vector<double>>,
  std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>, double, Eigen::Matrix<double, Eigen::Dynamic, 1>>,
  std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, double>,
  std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, std::vector<double>>,
  std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>>,
  std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Matrix<double, Eigen::Dynamic, 1>, double>,
  std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>>,
  std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Matrix<double, Eigen::Dynamic, 1>>,
  std::tuple<std::vector<double>, double, double, double>,
  std::tuple<std::vector<double>, double, double, std::vector<double>>,
  std::tuple<std::vector<double>, double, double, Eigen::Matrix<double, Eigen::Dynamic, 1>>,
  std::tuple<std::vector<double>, double, std::vector<double>, double>,
  std::tuple<std::vector<double>, double, std::vector<double>, std::vector<double>>,
  std::tuple<std::vector<double>, double, std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>>,
  std::tuple<std::vector<double>, double, Eigen::Matrix<double, Eigen::Dynamic, 1>, double>,
  std::tuple<std::vector<double>, double, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>>,
  std::tuple<std::vector<double>, double, Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Matrix<double, Eigen::Dynamic, 1>>,
  std::tuple<std::vector<double>, std::vector<double>, double, double>,
  std::tuple<std::vector<double>, std::vector<double>, double, std::vector<double>>,
  std::tuple<std::vector<double>, std::vector<double>, double, Eigen::Matrix<double, Eigen::Dynamic, 1>>,
  std::tuple<std::vector<double>, std::vector<double>, std::vector<double>, double>,
  std::tuple<std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>>,
  std::tuple<std::vector<double>, std::vector<double>, std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>>,
  std::tuple<std::vector<double>, std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, double>,
  std::tuple<std::vector<double>, std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>>,
  std::tuple<std::vector<double>, std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Matrix<double, Eigen::Dynamic, 1>>,
  std::tuple<std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, double, double>,
  std::tuple<std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, double, std::vector<double>>,
  std::tuple<std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, double, Eigen::Matrix<double, Eigen::Dynamic, 1>>,
  std::tuple<std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, double>,
  std::tuple<std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, std::vector<double>>,
  std::tuple<std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>>,
  std::tuple<std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Matrix<double, Eigen::Dynamic, 1>, double>,
  std::tuple<std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>>,
  std::tuple<std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Matrix<double, Eigen::Dynamic, 1>>,
  std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double, double, double>,
  std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double, double, std::vector<double>>,
  std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double, double, Eigen::Matrix<double, Eigen::Dynamic, 1>>,
  std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double, std::vector<double>, double>,
  std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double, std::vector<double>, std::vector<double>>,
  std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double, std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>>,
  std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double, Eigen::Matrix<double, Eigen::Dynamic, 1>, double>,
  std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>>,
  std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double, Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Matrix<double, Eigen::Dynamic, 1>>,
  std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, double, double>,
  std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, double, std::vector<double>>,
  std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, double, Eigen::Matrix<double, Eigen::Dynamic, 1>>,
  std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, std::vector<double>, double>,
  std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, std::vector<double>, std::vector<double>>,
  std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>>,
  std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, double>,
  std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>>,
  std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Matrix<double, Eigen::Dynamic, 1>>,
  std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Matrix<double, Eigen::Dynamic, 1>, double, double>,
  std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Matrix<double, Eigen::Dynamic, 1>, double, std::vector<double>>,
  std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Matrix<double, Eigen::Dynamic, 1>, double, Eigen::Matrix<double, Eigen::Dynamic, 1>>,
  std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, double>,
  std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, std::vector<double>>,
  std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>>,
  std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Matrix<double, Eigen::Dynamic, 1>, double>,
  std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>>,
  std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Matrix<double, Eigen::Dynamic, 1>>>;



template <size_t sizeof_args>
using base_set_t = std::conditional_t<sizeof_args == 1, unary_base_set,
 std::conditional_t<sizeof_args == 2, binary_base_set,
 std::conditional_t<sizeof_args == 3, ternary_base_sets,
 std::conditional_t<sizeof_args == 4, quad_base_set, quad_base_set>>>>;



    template <typename F, typename... EigVecs>
    void expect_ad_distribution(F&& f, const EigVecs&... args) {
      // Need to test all combinations of scalars, Eigen vectors, and std vectors
      using base_set = base_set_t<sizeof...(args)>;
      // TODO Write a version of for each that just takes in a template
      // and not an object
      stan::math::for_each([&f, &args...](auto&& type_tuple) {
        stan::math::apply([&f, &args...](auto&&... types) {
          stan::test::expect_ad(ad_tolerances{}, f, do_cast<std::decay_t<decltype(types)>>(args)...);
        }, type_tuple);
      }, base_set{});
    }

  }
}
#endif
