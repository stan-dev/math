#include <gtest/gtest.h>
#include <tuple>
#include <test/prob/test_fixture_distr.hpp>
#include <test/prob/test_fixture_cdf.hpp>
#include <test/prob/test_fixture_cdf_log.hpp>
#include <test/prob/test_fixture_ccdf_log.hpp>

typedef std::tuple<double, var, empty, empty, empty, empty> type_vv_real_real_0;
typedef std::tuple<double, std::vector<var>, empty, empty, empty, empty> type_vv_real_real_1;
typedef std::tuple<double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty, empty, empty, empty> type_vv_real_real_2;
typedef std::tuple<std::vector<double>, var, empty, empty, empty, empty> type_vv_real_real_3;
typedef std::tuple<std::vector<double>, std::vector<var>, empty, empty, empty, empty> type_vv_real_real_4;
typedef std::tuple<std::vector<double>, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty, empty, empty, empty> type_vv_real_real_5;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty, empty, empty, empty> type_vv_real_real_6;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, empty, empty, empty, empty> type_vv_real_real_7;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty, empty, empty, empty> type_vv_real_real_8;
typedef std::tuple<var, double, empty, empty, empty, empty> type_vv_real_real_9;
typedef std::tuple<var, std::vector<double>, empty, empty, empty, empty> type_vv_real_real_10;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty, empty, empty> type_vv_real_real_11;
typedef std::tuple<var, var, empty, empty, empty, empty> type_vv_real_real_12;
typedef std::tuple<var, std::vector<var>, empty, empty, empty, empty> type_vv_real_real_13;
typedef std::tuple<var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty, empty, empty, empty> type_vv_real_real_14;
typedef std::tuple<std::vector<var>, double, empty, empty, empty, empty> type_vv_real_real_15;
typedef std::tuple<std::vector<var>, std::vector<double>, empty, empty, empty, empty> type_vv_real_real_16;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty, empty, empty> type_vv_real_real_17;
typedef std::tuple<std::vector<var>, var, empty, empty, empty, empty> type_vv_real_real_18;
typedef std::tuple<std::vector<var>, std::vector<var>, empty, empty, empty, empty> type_vv_real_real_19;
typedef std::tuple<std::vector<var>, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty, empty, empty, empty> type_vv_real_real_20;
typedef std::tuple<stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double, empty, empty, empty, empty> type_vv_real_real_21;
typedef std::tuple<stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, std::vector<double>, empty, empty, empty, empty> type_vv_real_real_22;
typedef std::tuple<stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty, empty, empty> type_vv_real_real_23;
typedef std::tuple<stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty, empty, empty, empty> type_vv_real_real_24;
typedef std::tuple<stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, std::vector<var>, empty, empty, empty, empty> type_vv_real_real_25;
typedef std::tuple<stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty, empty, empty, empty> type_vv_real_real_26;

