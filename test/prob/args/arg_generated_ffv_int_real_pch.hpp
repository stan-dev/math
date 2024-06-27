#include <gtest/gtest.h>
#include <tuple>
#include <test/prob/test_fixture_distr.hpp>
#include <test/prob/test_fixture_cdf.hpp>
#include <test/prob/test_fixture_cdf_log.hpp>
#include <test/prob/test_fixture_ccdf_log.hpp>

typedef std::tuple<int, fvar<fvar<var>>, empty, empty, empty, empty>
    type_ffv_int_real_0;
typedef std::tuple<int, std::vector<fvar<fvar<var>>>, empty, empty, empty,
                   empty>
    type_ffv_int_real_1;
typedef std::tuple<int, Eigen::Matrix<fvar<fvar<var>>, Eigen::Dynamic, 1>,
                   empty, empty, empty, empty>
    type_ffv_int_real_2;
typedef std::tuple<std::vector<int>, fvar<fvar<var>>, empty, empty, empty,
                   empty>
    type_ffv_int_real_3;
typedef std::tuple<std::vector<int>, std::vector<fvar<fvar<var>>>, empty, empty,
                   empty, empty>
    type_ffv_int_real_4;
typedef std::tuple<std::vector<int>,
                   Eigen::Matrix<fvar<fvar<var>>, Eigen::Dynamic, 1>, empty,
                   empty, empty, empty>
    type_ffv_int_real_5;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, fvar<fvar<var>>,
                   empty, empty, empty, empty>
    type_ffv_int_real_6;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<fvar<fvar<var>>>, empty, empty, empty, empty>
    type_ffv_int_real_7;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<fvar<fvar<var>>, Eigen::Dynamic, 1>, empty,
                   empty, empty, empty>
    type_ffv_int_real_8;
