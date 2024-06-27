#include <gtest/gtest.h>
#include <tuple>
#include <test/prob/test_fixture_distr.hpp>
#include <test/prob/test_fixture_cdf.hpp>
#include <test/prob/test_fixture_cdf_log.hpp>
#include <test/prob/test_fixture_ccdf_log.hpp>

typedef std::tuple<int, double, empty, empty, empty, empty> type_v_int_real_0;
typedef std::tuple<int, std::vector<double>, empty, empty, empty, empty>
    type_v_int_real_1;
typedef std::tuple<int, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty,
                   empty, empty>
    type_v_int_real_2;
typedef std::tuple<int, var, empty, empty, empty, empty> type_v_int_real_3;
typedef std::tuple<int, std::vector<var>, empty, empty, empty, empty>
    type_v_int_real_4;
typedef std::tuple<int, Eigen::Matrix<var, Eigen::Dynamic, 1>, empty, empty,
                   empty, empty>
    type_v_int_real_5;
typedef std::tuple<std::vector<int>, double, empty, empty, empty, empty>
    type_v_int_real_6;
typedef std::tuple<std::vector<int>, std::vector<double>, empty, empty, empty,
                   empty>
    type_v_int_real_7;
typedef std::tuple<std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   empty, empty, empty, empty>
    type_v_int_real_8;
typedef std::tuple<std::vector<int>, var, empty, empty, empty, empty>
    type_v_int_real_9;
typedef std::tuple<std::vector<int>, std::vector<var>, empty, empty, empty,
                   empty>
    type_v_int_real_10;
typedef std::tuple<std::vector<int>, Eigen::Matrix<var, Eigen::Dynamic, 1>,
                   empty, empty, empty, empty>
    type_v_int_real_11;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, double, empty, empty,
                   empty, empty>
    type_v_int_real_12;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   empty, empty, empty, empty>
    type_v_int_real_13;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty,
                   empty, empty>
    type_v_int_real_14;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, var, empty, empty,
                   empty, empty>
    type_v_int_real_15;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   empty, empty, empty, empty>
    type_v_int_real_16;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<var, Eigen::Dynamic, 1>, empty, empty, empty,
                   empty>
    type_v_int_real_17;
