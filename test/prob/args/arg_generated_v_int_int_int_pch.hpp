#include <gtest/gtest.h>
#include <tuple>
#include <test/prob/test_fixture_distr.hpp>
#include <test/prob/test_fixture_cdf.hpp>
#include <test/prob/test_fixture_cdf_log.hpp>
#include <test/prob/test_fixture_ccdf_log.hpp>

typedef std::tuple<int, int, int, empty, empty, empty> type_v_int_int_int_0;
typedef std::tuple<int, int, std::vector<int>, empty, empty, empty>
    type_v_int_int_int_1;
typedef std::tuple<int, int, Eigen::Matrix<int, Eigen::Dynamic, 1>, empty,
                   empty, empty>
    type_v_int_int_int_2;
typedef std::tuple<int, std::vector<int>, int, empty, empty, empty>
    type_v_int_int_int_3;
typedef std::tuple<int, std::vector<int>, std::vector<int>, empty, empty, empty>
    type_v_int_int_int_4;
typedef std::tuple<int, std::vector<int>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   empty, empty, empty>
    type_v_int_int_int_5;
typedef std::tuple<int, Eigen::Matrix<int, Eigen::Dynamic, 1>, int, empty,
                   empty, empty>
    type_v_int_int_int_6;
typedef std::tuple<int, Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<int>,
                   empty, empty, empty>
    type_v_int_int_int_7;
typedef std::tuple<int, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, empty, empty, empty>
    type_v_int_int_int_8;
typedef std::tuple<std::vector<int>, int, int, empty, empty, empty>
    type_v_int_int_int_9;
typedef std::tuple<std::vector<int>, int, std::vector<int>, empty, empty, empty>
    type_v_int_int_int_10;
typedef std::tuple<std::vector<int>, int, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   empty, empty, empty>
    type_v_int_int_int_11;
typedef std::tuple<std::vector<int>, std::vector<int>, int, empty, empty, empty>
    type_v_int_int_int_12;
typedef std::tuple<std::vector<int>, std::vector<int>, std::vector<int>, empty,
                   empty, empty>
    type_v_int_int_int_13;
typedef std::tuple<std::vector<int>, std::vector<int>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, empty, empty, empty>
    type_v_int_int_int_14;
typedef std::tuple<std::vector<int>, Eigen::Matrix<int, Eigen::Dynamic, 1>, int,
                   empty, empty, empty>
    type_v_int_int_int_15;
typedef std::tuple<std::vector<int>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<int>, empty, empty, empty>
    type_v_int_int_int_16;
typedef std::tuple<std::vector<int>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, empty, empty, empty>
    type_v_int_int_int_17;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, int, int, empty,
                   empty, empty>
    type_v_int_int_int_18;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, int, std::vector<int>,
                   empty, empty, empty>
    type_v_int_int_int_19;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, int,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, empty, empty, empty>
    type_v_int_int_int_20;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<int>, int,
                   empty, empty, empty>
    type_v_int_int_int_21;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<int>,
                   std::vector<int>, empty, empty, empty>
    type_v_int_int_int_22;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<int>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, empty, empty, empty>
    type_v_int_int_int_23;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, int, empty, empty,
                   empty>
    type_v_int_int_int_24;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<int>,
                   empty, empty, empty>
    type_v_int_int_int_25;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, empty, empty, empty>
    type_v_int_int_int_26;
