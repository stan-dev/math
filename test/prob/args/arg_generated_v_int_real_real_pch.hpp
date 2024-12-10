#include <gtest/gtest.h>
#include <tuple>
#include <test/prob/test_fixture_distr.hpp>
#include <test/prob/test_fixture_cdf.hpp>
#include <test/prob/test_fixture_cdf_log.hpp>
#include <test/prob/test_fixture_ccdf_log.hpp>

typedef std::tuple<int, double, double, empty, empty, empty>
    type_v_int_real_real_0;
typedef std::tuple<int, double, std::vector<double>, empty, empty, empty>
    type_v_int_real_real_1;
typedef std::tuple<int, double, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty,
                   empty, empty>
    type_v_int_real_real_2;
typedef std::tuple<int, double, var, empty, empty, empty>
    type_v_int_real_real_3;
typedef std::tuple<int, double, std::vector<var>, empty, empty, empty>
    type_v_int_real_real_4;
typedef std::tuple<int, double, Eigen::Matrix<var, Eigen::Dynamic, 1>, empty,
                   empty, empty>
    type_v_int_real_real_5;
typedef std::tuple<int, std::vector<double>, double, empty, empty, empty>
    type_v_int_real_real_6;
typedef std::tuple<int, std::vector<double>, std::vector<double>, empty, empty,
                   empty>
    type_v_int_real_real_7;
typedef std::tuple<int, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty,
                   empty>
    type_v_int_real_real_8;
typedef std::tuple<int, std::vector<double>, var, empty, empty, empty>
    type_v_int_real_real_9;
typedef std::tuple<int, std::vector<double>, std::vector<var>, empty, empty,
                   empty>
    type_v_int_real_real_10;
typedef std::tuple<int, std::vector<double>,
                   Eigen::Matrix<var, Eigen::Dynamic, 1>, empty, empty, empty>
    type_v_int_real_real_11;
typedef std::tuple<int, Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty,
                   empty, empty>
    type_v_int_real_real_12;
typedef std::tuple<int, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty, empty, empty>
    type_v_int_real_real_13;
typedef std::tuple<int, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty,
                   empty>
    type_v_int_real_real_14;
typedef std::tuple<int, Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty,
                   empty, empty>
    type_v_int_real_real_15;
typedef std::tuple<int, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<var>, empty, empty, empty>
    type_v_int_real_real_16;
typedef std::tuple<int, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<var, Eigen::Dynamic, 1>, empty, empty, empty>
    type_v_int_real_real_17;
typedef std::tuple<int, var, double, empty, empty, empty>
    type_v_int_real_real_18;
typedef std::tuple<int, var, std::vector<double>, empty, empty, empty>
    type_v_int_real_real_19;
typedef std::tuple<int, var, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty,
                   empty, empty>
    type_v_int_real_real_20;
typedef std::tuple<int, var, var, empty, empty, empty> type_v_int_real_real_21;
typedef std::tuple<int, var, std::vector<var>, empty, empty, empty>
    type_v_int_real_real_22;
typedef std::tuple<int, var, Eigen::Matrix<var, Eigen::Dynamic, 1>, empty,
                   empty, empty>
    type_v_int_real_real_23;
typedef std::tuple<int, std::vector<var>, double, empty, empty, empty>
    type_v_int_real_real_24;
typedef std::tuple<int, std::vector<var>, std::vector<double>, empty, empty,
                   empty>
    type_v_int_real_real_25;
typedef std::tuple<int, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty,
                   empty>
    type_v_int_real_real_26;
typedef std::tuple<int, std::vector<var>, var, empty, empty, empty>
    type_v_int_real_real_27;
typedef std::tuple<int, std::vector<var>, std::vector<var>, empty, empty, empty>
    type_v_int_real_real_28;
typedef std::tuple<int, std::vector<var>, Eigen::Matrix<var, Eigen::Dynamic, 1>,
                   empty, empty, empty>
    type_v_int_real_real_29;
typedef std::tuple<int, Eigen::Matrix<var, Eigen::Dynamic, 1>, double, empty,
                   empty, empty>
    type_v_int_real_real_30;
typedef std::tuple<int, Eigen::Matrix<var, Eigen::Dynamic, 1>,
                   std::vector<double>, empty, empty, empty>
    type_v_int_real_real_31;
typedef std::tuple<int, Eigen::Matrix<var, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty,
                   empty>
    type_v_int_real_real_32;
typedef std::tuple<int, Eigen::Matrix<var, Eigen::Dynamic, 1>, var, empty,
                   empty, empty>
    type_v_int_real_real_33;
typedef std::tuple<int, Eigen::Matrix<var, Eigen::Dynamic, 1>, std::vector<var>,
                   empty, empty, empty>
    type_v_int_real_real_34;
typedef std::tuple<int, Eigen::Matrix<var, Eigen::Dynamic, 1>,
                   Eigen::Matrix<var, Eigen::Dynamic, 1>, empty, empty, empty>
    type_v_int_real_real_35;
typedef std::tuple<std::vector<int>, double, double, empty, empty, empty>
    type_v_int_real_real_36;
typedef std::tuple<std::vector<int>, double, std::vector<double>, empty, empty,
                   empty>
    type_v_int_real_real_37;
typedef std::tuple<std::vector<int>, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty,
                   empty>
    type_v_int_real_real_38;
typedef std::tuple<std::vector<int>, double, var, empty, empty, empty>
    type_v_int_real_real_39;
typedef std::tuple<std::vector<int>, double, std::vector<var>, empty, empty,
                   empty>
    type_v_int_real_real_40;
typedef std::tuple<std::vector<int>, double,
                   Eigen::Matrix<var, Eigen::Dynamic, 1>, empty, empty, empty>
    type_v_int_real_real_41;
typedef std::tuple<std::vector<int>, std::vector<double>, double, empty, empty,
                   empty>
    type_v_int_real_real_42;
typedef std::tuple<std::vector<int>, std::vector<double>, std::vector<double>,
                   empty, empty, empty>
    type_v_int_real_real_43;
typedef std::tuple<std::vector<int>, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty,
                   empty>
    type_v_int_real_real_44;
typedef std::tuple<std::vector<int>, std::vector<double>, var, empty, empty,
                   empty>
    type_v_int_real_real_45;
typedef std::tuple<std::vector<int>, std::vector<double>, std::vector<var>,
                   empty, empty, empty>
    type_v_int_real_real_46;
typedef std::tuple<std::vector<int>, std::vector<double>,
                   Eigen::Matrix<var, Eigen::Dynamic, 1>, empty, empty, empty>
    type_v_int_real_real_47;
typedef std::tuple<std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   double, empty, empty, empty>
    type_v_int_real_real_48;
typedef std::tuple<std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty, empty, empty>
    type_v_int_real_real_49;
typedef std::tuple<std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty,
                   empty>
    type_v_int_real_real_50;
typedef std::tuple<std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   var, empty, empty, empty>
    type_v_int_real_real_51;
typedef std::tuple<std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<var>, empty, empty, empty>
    type_v_int_real_real_52;
typedef std::tuple<std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<var, Eigen::Dynamic, 1>, empty, empty, empty>
    type_v_int_real_real_53;
typedef std::tuple<std::vector<int>, var, double, empty, empty, empty>
    type_v_int_real_real_54;
typedef std::tuple<std::vector<int>, var, std::vector<double>, empty, empty,
                   empty>
    type_v_int_real_real_55;
typedef std::tuple<std::vector<int>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty,
                   empty>
    type_v_int_real_real_56;
typedef std::tuple<std::vector<int>, var, var, empty, empty, empty>
    type_v_int_real_real_57;
typedef std::tuple<std::vector<int>, var, std::vector<var>, empty, empty, empty>
    type_v_int_real_real_58;
typedef std::tuple<std::vector<int>, var, Eigen::Matrix<var, Eigen::Dynamic, 1>,
                   empty, empty, empty>
    type_v_int_real_real_59;
typedef std::tuple<std::vector<int>, std::vector<var>, double, empty, empty,
                   empty>
    type_v_int_real_real_60;
typedef std::tuple<std::vector<int>, std::vector<var>, std::vector<double>,
                   empty, empty, empty>
    type_v_int_real_real_61;
typedef std::tuple<std::vector<int>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty,
                   empty>
    type_v_int_real_real_62;
typedef std::tuple<std::vector<int>, std::vector<var>, var, empty, empty, empty>
    type_v_int_real_real_63;
typedef std::tuple<std::vector<int>, std::vector<var>, std::vector<var>, empty,
                   empty, empty>
    type_v_int_real_real_64;
typedef std::tuple<std::vector<int>, std::vector<var>,
                   Eigen::Matrix<var, Eigen::Dynamic, 1>, empty, empty, empty>
    type_v_int_real_real_65;
typedef std::tuple<std::vector<int>, Eigen::Matrix<var, Eigen::Dynamic, 1>,
                   double, empty, empty, empty>
    type_v_int_real_real_66;
typedef std::tuple<std::vector<int>, Eigen::Matrix<var, Eigen::Dynamic, 1>,
                   std::vector<double>, empty, empty, empty>
    type_v_int_real_real_67;
typedef std::tuple<std::vector<int>, Eigen::Matrix<var, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty,
                   empty>
    type_v_int_real_real_68;
typedef std::tuple<std::vector<int>, Eigen::Matrix<var, Eigen::Dynamic, 1>, var,
                   empty, empty, empty>
    type_v_int_real_real_69;
typedef std::tuple<std::vector<int>, Eigen::Matrix<var, Eigen::Dynamic, 1>,
                   std::vector<var>, empty, empty, empty>
    type_v_int_real_real_70;
typedef std::tuple<std::vector<int>, Eigen::Matrix<var, Eigen::Dynamic, 1>,
                   Eigen::Matrix<var, Eigen::Dynamic, 1>, empty, empty, empty>
    type_v_int_real_real_71;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, double, double, empty,
                   empty, empty>
    type_v_int_real_real_72;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
                   std::vector<double>, empty, empty, empty>
    type_v_int_real_real_73;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty,
                   empty>
    type_v_int_real_real_74;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, double, var, empty,
                   empty, empty>
    type_v_int_real_real_75;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
                   std::vector<var>, empty, empty, empty>
    type_v_int_real_real_76;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
                   Eigen::Matrix<var, Eigen::Dynamic, 1>, empty, empty, empty>
    type_v_int_real_real_77;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   double, empty, empty, empty>
    type_v_int_real_real_78;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   std::vector<double>, empty, empty, empty>
    type_v_int_real_real_79;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty,
                   empty>
    type_v_int_real_real_80;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   var, empty, empty, empty>
    type_v_int_real_real_81;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   std::vector<var>, empty, empty, empty>
    type_v_int_real_real_82;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   Eigen::Matrix<var, Eigen::Dynamic, 1>, empty, empty, empty>
    type_v_int_real_real_83;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty,
                   empty, empty>
    type_v_int_real_real_84;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty, empty, empty>
    type_v_int_real_real_85;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty,
                   empty>
    type_v_int_real_real_86;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty, empty,
                   empty>
    type_v_int_real_real_87;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty, empty, empty>
    type_v_int_real_real_88;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<var, Eigen::Dynamic, 1>, empty, empty, empty>
    type_v_int_real_real_89;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, var, double, empty,
                   empty, empty>
    type_v_int_real_real_90;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
                   std::vector<double>, empty, empty, empty>
    type_v_int_real_real_91;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty,
                   empty>
    type_v_int_real_real_92;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, var, var, empty,
                   empty, empty>
    type_v_int_real_real_93;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, var, std::vector<var>,
                   empty, empty, empty>
    type_v_int_real_real_94;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<var, Eigen::Dynamic, 1>, empty, empty, empty>
    type_v_int_real_real_95;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   double, empty, empty, empty>
    type_v_int_real_real_96;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<double>, empty, empty, empty>
    type_v_int_real_real_97;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty,
                   empty>
    type_v_int_real_real_98;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>, var,
                   empty, empty, empty>
    type_v_int_real_real_99;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<var>, empty, empty, empty>
    type_v_int_real_real_100;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   Eigen::Matrix<var, Eigen::Dynamic, 1>, empty, empty, empty>
    type_v_int_real_real_101;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<var, Eigen::Dynamic, 1>, double, empty, empty,
                   empty>
    type_v_int_real_real_102;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<var, Eigen::Dynamic, 1>, std::vector<double>,
                   empty, empty, empty>
    type_v_int_real_real_103;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<var, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty,
                   empty>
    type_v_int_real_real_104;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<var, Eigen::Dynamic, 1>, var, empty, empty,
                   empty>
    type_v_int_real_real_105;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<var, Eigen::Dynamic, 1>, std::vector<var>,
                   empty, empty, empty>
    type_v_int_real_real_106;
typedef std::tuple<Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<var, Eigen::Dynamic, 1>,
                   Eigen::Matrix<var, Eigen::Dynamic, 1>, empty, empty, empty>
    type_v_int_real_real_107;
