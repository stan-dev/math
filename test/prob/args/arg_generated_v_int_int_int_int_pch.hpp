#include <gtest/gtest.h>
#include <tuple>
#include <test/prob/test_fixture_distr.hpp>
#include <test/prob/test_fixture_cdf.hpp>
#include <test/prob/test_fixture_cdf_log.hpp>
#include <test/prob/test_fixture_ccdf_log.hpp>

typedef std::tuple<int, int, int, int, empty, empty> type_v_int_int_int_int_0;
typedef std::tuple<int, int, int, std::vector<int>, empty, empty> type_v_int_int_int_int_1;
typedef std::tuple<int, int, int, Eigen::Matrix<int, Eigen::Dynamic, 1>, empty, empty> type_v_int_int_int_int_2;
typedef std::tuple<int, int, std::vector<int>, int, empty, empty> type_v_int_int_int_int_3;
typedef std::tuple<int, int, std::vector<int>, std::vector<int>, empty, empty> type_v_int_int_int_int_4;
typedef std::tuple<int, int, std::vector<int>, Eigen::Matrix<int, Eigen::Dynamic, 1>, empty, empty> type_v_int_int_int_int_5;
typedef std::tuple<int, int, Eigen::Matrix<int, Eigen::Dynamic, 1>, int, empty, empty> type_v_int_int_int_int_6;
typedef std::tuple<int, int, Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<int>, empty, empty> type_v_int_int_int_int_7;
typedef std::tuple<int, int, Eigen::Matrix<int, Eigen::Dynamic, 1>, Eigen::Matrix<int, Eigen::Dynamic, 1>, empty, empty> type_v_int_int_int_int_8;
typedef std::tuple<int, std::vector<int>, int, int, empty, empty> type_v_int_int_int_int_9;
typedef std::tuple<int, std::vector<int>, int, std::vector<int>, empty, empty> type_v_int_int_int_int_10;
typedef std::tuple<int, std::vector<int>, int, Eigen::Matrix<int, Eigen::Dynamic, 1>, empty, empty> type_v_int_int_int_int_11;
typedef std::tuple<int, std::vector<int>, std::vector<int>, int, empty, empty> type_v_int_int_int_int_12;
typedef std::tuple<int, std::vector<int>, std::vector<int>, std::vector<int>, empty, empty> type_v_int_int_int_int_13;
typedef std::tuple<int, std::vector<int>, std::vector<int>, Eigen::Matrix<int, Eigen::Dynamic, 1>, empty, empty> type_v_int_int_int_int_14;
typedef std::tuple<int, std::vector<int>, Eigen::Matrix<int, Eigen::Dynamic, 1>, int, empty, empty> type_v_int_int_int_int_15;
typedef std::tuple<int, std::vector<int>, Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<int>, empty, empty> type_v_int_int_int_int_16;
typedef std::tuple<int, std::vector<int>, Eigen::Matrix<int, Eigen::Dynamic, 1>, Eigen::Matrix<int, Eigen::Dynamic, 1>, empty, empty> type_v_int_int_int_int_17;
typedef std::tuple<int, Eigen::Matrix<int, Eigen::Dynamic, 1>, int, int, empty, empty> type_v_int_int_int_int_18;
typedef std::tuple<int, Eigen::Matrix<int, Eigen::Dynamic, 1>, int, std::vector<int>, empty, empty> type_v_int_int_int_int_19;
typedef std::tuple<int, Eigen::Matrix<int, Eigen::Dynamic, 1>, int, Eigen::Matrix<int, Eigen::Dynamic, 1>, empty, empty> type_v_int_int_int_int_20;
typedef std::tuple<int, Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<int>, int, empty, empty> type_v_int_int_int_int_21;
typedef std::tuple<int, Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<int>, std::vector<int>, empty, empty> type_v_int_int_int_int_22;
typedef std::tuple<int, Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<int>, Eigen::Matrix<int, Eigen::Dynamic, 1>, empty, empty> type_v_int_int_int_int_23;
typedef std::tuple<int, Eigen::Matrix<int, Eigen::Dynamic, 1>, Eigen::Matrix<int, Eigen::Dynamic, 1>, int, empty, empty> type_v_int_int_int_int_24;
typedef std::tuple<int, Eigen::Matrix<int, Eigen::Dynamic, 1>, Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<int>, empty, empty> type_v_int_int_int_int_25;
typedef std::tuple<int, Eigen::Matrix<int, Eigen::Dynamic, 1>, Eigen::Matrix<int, Eigen::Dynamic, 1>, Eigen::Matrix<int, Eigen::Dynamic, 1>, empty, empty> type_v_int_int_int_int_26;
typedef std::tuple<std::vector<int>, int, int, int, empty, empty> type_v_int_int_int_int_27;
typedef std::tuple<std::vector<int>, int, int, std::vector<int>, empty, empty> type_v_int_int_int_int_28;
typedef std::tuple<std::vector<int>, int, int, Eigen::Matrix<int, Eigen::Dynamic, 1>, empty, empty> type_v_int_int_int_int_29;
typedef std::tuple<std::vector<int>, int, std::vector<int>, int, empty, empty> type_v_int_int_int_int_30;
typedef std::tuple<std::vector<int>, int, std::vector<int>, std::vector<int>, empty, empty> type_v_int_int_int_int_31;
typedef std::tuple<std::vector<int>, int, std::vector<int>, Eigen::Matrix<int, Eigen::Dynamic, 1>, empty, empty> type_v_int_int_int_int_32;
typedef std::tuple<std::vector<int>, int, Eigen::Matrix<int, Eigen::Dynamic, 1>, int, empty, empty> type_v_int_int_int_int_33;
typedef std::tuple<std::vector<int>, int, Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<int>, empty, empty> type_v_int_int_int_int_34;
typedef std::tuple<std::vector<int>, int, Eigen::Matrix<int, Eigen::Dynamic, 1>, Eigen::Matrix<int, Eigen::Dynamic, 1>, empty, empty> type_v_int_int_int_int_35;
typedef std::tuple<std::vector<int>, std::vector<int>, int, int, empty, empty> type_v_int_int_int_int_36;
typedef std::tuple<std::vector<int>, std::vector<int>, int, std::vector<int>, empty, empty> type_v_int_int_int_int_37;
typedef std::tuple<std::vector<int>, std::vector<int>, int, Eigen::Matrix<int, Eigen::Dynamic, 1>, empty, empty> type_v_int_int_int_int_38;
typedef std::tuple<std::vector<int>, std::vector<int>, std::vector<int>, int, empty, empty> type_v_int_int_int_int_39;
typedef std::tuple<std::vector<int>, std::vector<int>, std::vector<int>, std::vector<int>, empty, empty> type_v_int_int_int_int_40;
typedef std::tuple<std::vector<int>, std::vector<int>, std::vector<int>, Eigen::Matrix<int, Eigen::Dynamic, 1>, empty, empty> type_v_int_int_int_int_41;
typedef std::tuple<std::vector<int>, std::vector<int>, Eigen::Matrix<int, Eigen::Dynamic, 1>, int, empty, empty> type_v_int_int_int_int_42;
typedef std::tuple<std::vector<int>, std::vector<int>, Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<int>, empty, empty> type_v_int_int_int_int_43;
typedef std::tuple<std::vector<int>, std::vector<int>, Eigen::Matrix<int, Eigen::Dynamic, 1>, Eigen::Matrix<int, Eigen::Dynamic, 1>, empty, empty> type_v_int_int_int_int_44;
typedef std::tuple<std::vector<int>, Eigen::Matrix<int, Eigen::Dynamic, 1>, int, int, empty, empty> type_v_int_int_int_int_45;
typedef std::tuple<std::vector<int>, Eigen::Matrix<int, Eigen::Dynamic, 1>, int, std::vector<int>, empty, empty> type_v_int_int_int_int_46;
typedef std::tuple<std::vector<int>, Eigen::Matrix<int, Eigen::Dynamic, 1>, int, Eigen::Matrix<int, Eigen::Dynamic, 1>, empty, empty> type_v_int_int_int_int_47;
typedef std::tuple<std::vector<int>, Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<int>, int, empty, empty> type_v_int_int_int_int_48;
typedef std::tuple<std::vector<int>, Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<int>, std::vector<int>, empty, empty> type_v_int_int_int_int_49;
typedef std::tuple<std::vector<int>, Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<int>, Eigen::Matrix<int, Eigen::Dynamic, 1>, empty, empty> type_v_int_int_int_int_50;
typedef std::tuple<std::vector<int>, Eigen::Matrix<int, Eigen::Dynamic, 1>, Eigen::Matrix<int, Eigen::Dynamic, 1>, int, empty, empty> type_v_int_int_int_int_51;
typedef std::tuple<std::vector<int>, Eigen::Matrix<int, Eigen::Dynamic, 1>, Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<int>, empty, empty> type_v_int_int_int_int_52;
typedef std::tuple<std::vector<int>, Eigen::Matrix<int, Eigen::Dynamic, 1>, Eigen::Matrix<int, Eigen::Dynamic, 1>, Eigen::Matrix<int, Eigen::Dynamic, 1>, empty, empty> type_v_int_int_int_int_53;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, int, int, int, empty, empty> type_v_int_int_int_int_54;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, int, int, std::vector<int>, empty, empty> type_v_int_int_int_int_55;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, int, int, Eigen::Matrix<int, Eigen::Dynamic, 1>, empty, empty> type_v_int_int_int_int_56;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, int, std::vector<int>, int, empty, empty> type_v_int_int_int_int_57;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, int, std::vector<int>, std::vector<int>, empty, empty> type_v_int_int_int_int_58;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, int, std::vector<int>, Eigen::Matrix<int, Eigen::Dynamic, 1>, empty, empty> type_v_int_int_int_int_59;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, int, Eigen::Matrix<int, Eigen::Dynamic, 1>, int, empty, empty> type_v_int_int_int_int_60;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, int, Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<int>, empty, empty> type_v_int_int_int_int_61;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, int, Eigen::Matrix<int, Eigen::Dynamic, 1>, Eigen::Matrix<int, Eigen::Dynamic, 1>, empty, empty> type_v_int_int_int_int_62;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<int>, int, int, empty, empty> type_v_int_int_int_int_63;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<int>, int, std::vector<int>, empty, empty> type_v_int_int_int_int_64;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<int>, int, Eigen::Matrix<int, Eigen::Dynamic, 1>, empty, empty> type_v_int_int_int_int_65;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<int>, std::vector<int>, int, empty, empty> type_v_int_int_int_int_66;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<int>, std::vector<int>, std::vector<int>, empty, empty> type_v_int_int_int_int_67;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<int>, std::vector<int>, Eigen::Matrix<int, Eigen::Dynamic, 1>, empty, empty> type_v_int_int_int_int_68;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<int>, Eigen::Matrix<int, Eigen::Dynamic, 1>, int, empty, empty> type_v_int_int_int_int_69;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<int>, Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<int>, empty, empty> type_v_int_int_int_int_70;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<int>, Eigen::Matrix<int, Eigen::Dynamic, 1>, Eigen::Matrix<int, Eigen::Dynamic, 1>, empty, empty> type_v_int_int_int_int_71;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, Eigen::Matrix<int, Eigen::Dynamic, 1>, int, int, empty, empty> type_v_int_int_int_int_72;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, Eigen::Matrix<int, Eigen::Dynamic, 1>, int, std::vector<int>, empty, empty> type_v_int_int_int_int_73;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, Eigen::Matrix<int, Eigen::Dynamic, 1>, int, Eigen::Matrix<int, Eigen::Dynamic, 1>, empty, empty> type_v_int_int_int_int_74;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<int>, int, empty, empty> type_v_int_int_int_int_75;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<int>, std::vector<int>, empty, empty> type_v_int_int_int_int_76;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<int>, Eigen::Matrix<int, Eigen::Dynamic, 1>, empty, empty> type_v_int_int_int_int_77;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, Eigen::Matrix<int, Eigen::Dynamic, 1>, Eigen::Matrix<int, Eigen::Dynamic, 1>, int, empty, empty> type_v_int_int_int_int_78;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, Eigen::Matrix<int, Eigen::Dynamic, 1>, Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<int>, empty, empty> type_v_int_int_int_int_79;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, Eigen::Matrix<int, Eigen::Dynamic, 1>, Eigen::Matrix<int, Eigen::Dynamic, 1>, Eigen::Matrix<int, Eigen::Dynamic, 1>, empty, empty> type_v_int_int_int_int_80;

