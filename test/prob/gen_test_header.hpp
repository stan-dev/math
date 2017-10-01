#include <gtest/gtest.h>
#include <boost/mpl/vector.hpp>
#include <test/prob/test_fixture_distr.hpp>
#include <test/prob/test_fixture_cdf.hpp>
#include <test/prob/test_fixture_cdf_log.hpp>
#include <test/prob/test_fixture_ccdf_log.hpp>

template <typename T>
using type_0 = boost::mpl::vector<int, T, empty, empty, empty, empty>;
template <typename T>
using type_1 = boost::mpl::vector<int, std::vector<T >, empty, empty, empty, empty>;
template <typename T>
using type_2 = boost::mpl::vector<int, Eigen::Matrix<T, Eigen::Dynamic, 1>, empty, empty, empty, empty>;
template <typename T>
using type_3 = boost::mpl::vector<int, Eigen::Matrix<T, 1, Eigen::Dynamic>, empty, empty, empty, empty>;
template <typename T>
using type_4 = boost::mpl::vector<std::vector<int>, T, empty, empty, empty, empty>;
template <typename T>
using type_5 = boost::mpl::vector<std::vector<int>, std::vector<T >, empty, empty, empty, empty>;
template <typename T>
using type_6 = boost::mpl::vector<std::vector<int>, Eigen::Matrix<T, Eigen::Dynamic, 1>, empty, empty, empty, empty>;
template <typename T>
using type_7 = boost::mpl::vector<std::vector<int>, Eigen::Matrix<T, 1, Eigen::Dynamic>, empty, empty, empty, empty>;
template <typename T>
using type_8 = boost::mpl::vector<Eigen::Matrix<int, Eigen::Dynamic, 1>, T, empty, empty, empty, empty>;
template <typename T>
using type_9 = boost::mpl::vector<Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<T >, empty, empty, empty, empty>;
template <typename T>
using type_10 = boost::mpl::vector<Eigen::Matrix<int, Eigen::Dynamic, 1>, Eigen::Matrix<T, Eigen::Dynamic, 1>, empty, empty, empty, empty>;
template <typename T>
using type_11 = boost::mpl::vector<Eigen::Matrix<int, Eigen::Dynamic, 1>, Eigen::Matrix<T, 1, Eigen::Dynamic>, empty, empty, empty, empty>;
template <typename T>
using type_12 = boost::mpl::vector<Eigen::Matrix<int, 1, Eigen::Dynamic>, T, empty, empty, empty, empty>;
template <typename T>
using type_13 = boost::mpl::vector<Eigen::Matrix<int, 1, Eigen::Dynamic>, std::vector<T >, empty, empty, empty, empty>;
template <typename T>
using type_14 = boost::mpl::vector<Eigen::Matrix<int, 1, Eigen::Dynamic>, Eigen::Matrix<T, Eigen::Dynamic, 1>, empty, empty, empty, empty>;
template <typename T>
using type_15 = boost::mpl::vector<Eigen::Matrix<int, 1, Eigen::Dynamic>, Eigen::Matrix<T, 1, Eigen::Dynamic>, empty, empty, empty, empty>;
