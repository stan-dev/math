#include <gtest/gtest.h>
#include <tuple>
#include <test/prob/test_fixture_distr.hpp>
#include <test/prob/test_fixture_cdf.hpp>
#include <test/prob/test_fixture_cdf_log.hpp>
#include <test/prob/test_fixture_ccdf_log.hpp>

typedef std::tuple<int, int, double, empty, empty, empty> type_v_int_int_real_0;
typedef std::tuple<int, int, std::vector<double>, empty, empty, empty> type_v_int_int_real_1;
typedef std::tuple<int, int, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty, empty> type_v_int_int_real_2;
typedef std::tuple<int, int, var, empty, empty, empty> type_v_int_int_real_3;
typedef std::tuple<int, int, std::vector<var>, empty, empty, empty> type_v_int_int_real_4;
typedef std::tuple<int, int, Eigen::Matrix<var, Eigen::Dynamic, 1>, empty, empty, empty> type_v_int_int_real_5;
typedef std::tuple<int, std::vector<int>, double, empty, empty, empty> type_v_int_int_real_6;
typedef std::tuple<int, std::vector<int>, std::vector<double>, empty, empty, empty> type_v_int_int_real_7;
typedef std::tuple<int, std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty, empty> type_v_int_int_real_8;
typedef std::tuple<int, std::vector<int>, var, empty, empty, empty> type_v_int_int_real_9;
typedef std::tuple<int, std::vector<int>, std::vector<var>, empty, empty, empty> type_v_int_int_real_10;
typedef std::tuple<int, std::vector<int>, Eigen::Matrix<var, Eigen::Dynamic, 1>, empty, empty, empty> type_v_int_int_real_11;
typedef std::tuple<int, Eigen::Matrix<int, Eigen::Dynamic, 1>, double, empty, empty, empty> type_v_int_int_real_12;
typedef std::tuple<int, Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>, empty, empty, empty> type_v_int_int_real_13;
typedef std::tuple<int, Eigen::Matrix<int, Eigen::Dynamic, 1>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty, empty> type_v_int_int_real_14;
typedef std::tuple<int, Eigen::Matrix<int, Eigen::Dynamic, 1>, var, empty, empty, empty> type_v_int_int_real_15;
typedef std::tuple<int, Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>, empty, empty, empty> type_v_int_int_real_16;
typedef std::tuple<int, Eigen::Matrix<int, Eigen::Dynamic, 1>, Eigen::Matrix<var, Eigen::Dynamic, 1>, empty, empty, empty> type_v_int_int_real_17;
typedef std::tuple<std::vector<int>, int, double, empty, empty, empty> type_v_int_int_real_18;
typedef std::tuple<std::vector<int>, int, std::vector<double>, empty, empty, empty> type_v_int_int_real_19;
typedef std::tuple<std::vector<int>, int, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty, empty> type_v_int_int_real_20;
typedef std::tuple<std::vector<int>, int, var, empty, empty, empty> type_v_int_int_real_21;
typedef std::tuple<std::vector<int>, int, std::vector<var>, empty, empty, empty> type_v_int_int_real_22;
typedef std::tuple<std::vector<int>, int, Eigen::Matrix<var, Eigen::Dynamic, 1>, empty, empty, empty> type_v_int_int_real_23;
typedef std::tuple<std::vector<int>, std::vector<int>, double, empty, empty, empty> type_v_int_int_real_24;
typedef std::tuple<std::vector<int>, std::vector<int>, std::vector<double>, empty, empty, empty> type_v_int_int_real_25;
typedef std::tuple<std::vector<int>, std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty, empty> type_v_int_int_real_26;
typedef std::tuple<std::vector<int>, std::vector<int>, var, empty, empty, empty> type_v_int_int_real_27;
typedef std::tuple<std::vector<int>, std::vector<int>, std::vector<var>, empty, empty, empty> type_v_int_int_real_28;
typedef std::tuple<std::vector<int>, std::vector<int>, Eigen::Matrix<var, Eigen::Dynamic, 1>, empty, empty, empty> type_v_int_int_real_29;
typedef std::tuple<std::vector<int>, Eigen::Matrix<int, Eigen::Dynamic, 1>, double, empty, empty, empty> type_v_int_int_real_30;
typedef std::tuple<std::vector<int>, Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>, empty, empty, empty> type_v_int_int_real_31;
typedef std::tuple<std::vector<int>, Eigen::Matrix<int, Eigen::Dynamic, 1>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty, empty> type_v_int_int_real_32;
typedef std::tuple<std::vector<int>, Eigen::Matrix<int, Eigen::Dynamic, 1>, var, empty, empty, empty> type_v_int_int_real_33;
typedef std::tuple<std::vector<int>, Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>, empty, empty, empty> type_v_int_int_real_34;
typedef std::tuple<std::vector<int>, Eigen::Matrix<int, Eigen::Dynamic, 1>, Eigen::Matrix<var, Eigen::Dynamic, 1>, empty, empty, empty> type_v_int_int_real_35;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, int, double, empty, empty, empty> type_v_int_int_real_36;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, int, std::vector<double>, empty, empty, empty> type_v_int_int_real_37;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, int, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty, empty> type_v_int_int_real_38;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, int, var, empty, empty, empty> type_v_int_int_real_39;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, int, std::vector<var>, empty, empty, empty> type_v_int_int_real_40;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, int, Eigen::Matrix<var, Eigen::Dynamic, 1>, empty, empty, empty> type_v_int_int_real_41;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<int>, double, empty, empty, empty> type_v_int_int_real_42;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<int>, std::vector<double>, empty, empty, empty> type_v_int_int_real_43;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty, empty> type_v_int_int_real_44;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<int>, var, empty, empty, empty> type_v_int_int_real_45;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<int>, std::vector<var>, empty, empty, empty> type_v_int_int_real_46;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<int>, Eigen::Matrix<var, Eigen::Dynamic, 1>, empty, empty, empty> type_v_int_int_real_47;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, Eigen::Matrix<int, Eigen::Dynamic, 1>, double, empty, empty, empty> type_v_int_int_real_48;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>, empty, empty, empty> type_v_int_int_real_49;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, Eigen::Matrix<int, Eigen::Dynamic, 1>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty, empty> type_v_int_int_real_50;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, Eigen::Matrix<int, Eigen::Dynamic, 1>, var, empty, empty, empty> type_v_int_int_real_51;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>, empty, empty, empty> type_v_int_int_real_52;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, Eigen::Matrix<int, Eigen::Dynamic, 1>, Eigen::Matrix<var, Eigen::Dynamic, 1>, empty, empty, empty> type_v_int_int_real_53;

