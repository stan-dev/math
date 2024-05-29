#include <gtest/gtest.h>
#include <tuple>
#include <test/prob/test_fixture_distr.hpp>
#include <test/prob/test_fixture_cdf.hpp>
#include <test/prob/test_fixture_cdf_log.hpp>
#include <test/prob/test_fixture_ccdf_log.hpp>

typedef std::tuple<int, int, var, empty, empty, empty> type_vv_int_int_real_0;
typedef std::tuple<int, int, std::vector<var>, empty, empty, empty>
    type_vv_int_int_real_1;
typedef std::tuple<
    int, int, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    empty, empty, empty>
    type_vv_int_int_real_2;
typedef std::tuple<int, std::vector<int>, var, empty, empty, empty>
    type_vv_int_int_real_3;
typedef std::tuple<int, std::vector<int>, std::vector<var>, empty, empty, empty>
    type_vv_int_int_real_4;
typedef std::tuple<
    int, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty, empty>
    type_vv_int_int_real_5;
typedef std::tuple<int, Eigen::Matrix<int, Eigen::Dynamic, 1>, var, empty,
                   empty, empty>
    type_vv_int_int_real_6;
typedef std::tuple<int, Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   empty, empty, empty>
    type_vv_int_int_real_7;
typedef std::tuple<
    int, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty, empty>
    type_vv_int_int_real_8;
typedef std::tuple<std::vector<int>, int, var, empty, empty, empty>
    type_vv_int_int_real_9;
typedef std::tuple<std::vector<int>, int, std::vector<var>, empty, empty, empty>
    type_vv_int_int_real_10;
typedef std::tuple<
    std::vector<int>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty, empty>
    type_vv_int_int_real_11;
typedef std::tuple<std::vector<int>, std::vector<int>, var, empty, empty, empty>
    type_vv_int_int_real_12;
typedef std::tuple<std::vector<int>, std::vector<int>, std::vector<var>, empty,
                   empty, empty>
    type_vv_int_int_real_13;
typedef std::tuple<
    std::vector<int>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty, empty>
    type_vv_int_int_real_14;
typedef std::tuple<std::vector<int>, Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
                   empty, empty, empty>
    type_vv_int_int_real_15;
typedef std::tuple<std::vector<int>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<var>, empty, empty, empty>
    type_vv_int_int_real_16;
typedef std::tuple<
    std::vector<int>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty, empty>
    type_vv_int_int_real_17;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, int, var, empty,
                   empty, empty>
    type_vv_int_int_real_18;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, int, std::vector<var>,
                   empty, empty, empty>
    type_vv_int_int_real_19;
typedef std::tuple<
    Eigen::Matrix<int, Eigen::Dynamic, 1>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty, empty>
    type_vv_int_int_real_20;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<int>, var,
                   empty, empty, empty>
    type_vv_int_int_real_21;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<int>,
                   std::vector<var>, empty, empty, empty>
    type_vv_int_int_real_22;
typedef std::tuple<
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty, empty>
    type_vv_int_int_real_23;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var, empty, empty,
                   empty>
    type_vv_int_int_real_24;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   empty, empty, empty>
    type_vv_int_int_real_25;
typedef std::tuple<
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty, empty>
    type_vv_int_int_real_26;
