#include <gtest/gtest.h>
#include <tuple>
#include <test/prob/test_fixture_distr.hpp>
#include <test/prob/test_fixture_cdf.hpp>
#include <test/prob/test_fixture_cdf_log.hpp>
#include <test/prob/test_fixture_ccdf_log.hpp>

typedef std::tuple<double, double, int, double, var, empty>
    type_vv_real_real_int_real_real_0;
typedef std::tuple<double, double, int, double, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1;
typedef std::tuple<
    double, double, int, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2;
typedef std::tuple<double, double, int, std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_3;
typedef std::tuple<double, double, int, std::vector<double>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_4;
typedef std::tuple<
    double, double, int, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_5;
typedef std::tuple<double, double, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_6;
typedef std::tuple<double, double, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_7;
typedef std::tuple<
    double, double, int, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_8;
typedef std::tuple<double, double, int, var, double, empty>
    type_vv_real_real_int_real_real_9;
typedef std::tuple<double, double, int, var, std::vector<double>, empty>
    type_vv_real_real_int_real_real_10;
typedef std::tuple<double, double, int, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_11;
typedef std::tuple<double, double, int, var, var, empty>
    type_vv_real_real_int_real_real_12;
typedef std::tuple<double, double, int, var, std::vector<var>, empty>
    type_vv_real_real_int_real_real_13;
typedef std::tuple<
    double, double, int, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_14;
typedef std::tuple<double, double, int, std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_15;
typedef std::tuple<double, double, int, std::vector<var>, std::vector<double>,
                   empty>
    type_vv_real_real_int_real_real_16;
typedef std::tuple<double, double, int, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_17;
typedef std::tuple<double, double, int, std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_18;
typedef std::tuple<double, double, int, std::vector<var>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_19;
typedef std::tuple<
    double, double, int, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_20;
typedef std::tuple<
    double, double, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_21;
typedef std::tuple<
    double, double, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_22;
typedef std::tuple<
    double, double, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_23;
typedef std::tuple<
    double, double, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_24;
typedef std::tuple<
    double, double, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_25;
typedef std::tuple<
    double, double, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_26;
typedef std::tuple<double, double, std::vector<int>, double, var, empty>
    type_vv_real_real_int_real_real_27;
typedef std::tuple<double, double, std::vector<int>, double, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_28;
typedef std::tuple<
    double, double, std::vector<int>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_29;
typedef std::tuple<double, double, std::vector<int>, std::vector<double>, var,
                   empty>
    type_vv_real_real_int_real_real_30;
typedef std::tuple<double, double, std::vector<int>, std::vector<double>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_31;
typedef std::tuple<
    double, double, std::vector<int>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_32;
typedef std::tuple<double, double, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_33;
typedef std::tuple<double, double, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_34;
typedef std::tuple<
    double, double, std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_35;
typedef std::tuple<double, double, std::vector<int>, var, double, empty>
    type_vv_real_real_int_real_real_36;
typedef std::tuple<double, double, std::vector<int>, var, std::vector<double>,
                   empty>
    type_vv_real_real_int_real_real_37;
typedef std::tuple<double, double, std::vector<int>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_38;
typedef std::tuple<double, double, std::vector<int>, var, var, empty>
    type_vv_real_real_int_real_real_39;
typedef std::tuple<double, double, std::vector<int>, var, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_40;
typedef std::tuple<
    double, double, std::vector<int>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_41;
typedef std::tuple<double, double, std::vector<int>, std::vector<var>, double,
                   empty>
    type_vv_real_real_int_real_real_42;
typedef std::tuple<double, double, std::vector<int>, std::vector<var>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_43;
typedef std::tuple<double, double, std::vector<int>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_44;
typedef std::tuple<double, double, std::vector<int>, std::vector<var>, var,
                   empty>
    type_vv_real_real_int_real_real_45;
typedef std::tuple<double, double, std::vector<int>, std::vector<var>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_46;
typedef std::tuple<
    double, double, std::vector<int>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_47;
typedef std::tuple<
    double, double, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_48;
typedef std::tuple<
    double, double, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_49;
typedef std::tuple<
    double, double, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_50;
typedef std::tuple<
    double, double, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_51;
typedef std::tuple<
    double, double, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_52;
typedef std::tuple<
    double, double, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_53;
typedef std::tuple<double, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   double, var, empty>
    type_vv_real_real_int_real_real_54;
typedef std::tuple<double, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   double, std::vector<var>, empty>
    type_vv_real_real_int_real_real_55;
typedef std::tuple<
    double, double, Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_56;
typedef std::tuple<double, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_57;
typedef std::tuple<double, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<double>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_58;
typedef std::tuple<
    double, double, Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_59;
typedef std::tuple<double, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_60;
typedef std::tuple<double, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_61;
typedef std::tuple<
    double, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_62;
typedef std::tuple<double, double, Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
                   double, empty>
    type_vv_real_real_int_real_real_63;
typedef std::tuple<double, double, Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_64;
typedef std::tuple<double, double, Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_65;
typedef std::tuple<double, double, Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
                   var, empty>
    type_vv_real_real_int_real_real_66;
typedef std::tuple<double, double, Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_67;
typedef std::tuple<
    double, double, Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_68;
typedef std::tuple<double, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_69;
typedef std::tuple<double, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<var>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_70;
typedef std::tuple<double, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   empty>
    type_vv_real_real_int_real_real_71;
typedef std::tuple<double, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_72;
typedef std::tuple<double, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<var>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_73;
typedef std::tuple<
    double, double, Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_74;
typedef std::tuple<
    double, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_75;
typedef std::tuple<
    double, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_76;
typedef std::tuple<
    double, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_77;
typedef std::tuple<
    double, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_78;
typedef std::tuple<
    double, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_79;
typedef std::tuple<
    double, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_80;
typedef std::tuple<double, std::vector<double>, int, double, var, empty>
    type_vv_real_real_int_real_real_81;
typedef std::tuple<double, std::vector<double>, int, double, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_82;
typedef std::tuple<
    double, std::vector<double>, int, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_83;
typedef std::tuple<double, std::vector<double>, int, std::vector<double>, var,
                   empty>
    type_vv_real_real_int_real_real_84;
typedef std::tuple<double, std::vector<double>, int, std::vector<double>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_85;
typedef std::tuple<
    double, std::vector<double>, int, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_86;
typedef std::tuple<double, std::vector<double>, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_87;
typedef std::tuple<double, std::vector<double>, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_88;
typedef std::tuple<
    double, std::vector<double>, int, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_89;
typedef std::tuple<double, std::vector<double>, int, var, double, empty>
    type_vv_real_real_int_real_real_90;
typedef std::tuple<double, std::vector<double>, int, var, std::vector<double>,
                   empty>
    type_vv_real_real_int_real_real_91;
typedef std::tuple<double, std::vector<double>, int, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_92;
typedef std::tuple<double, std::vector<double>, int, var, var, empty>
    type_vv_real_real_int_real_real_93;
typedef std::tuple<double, std::vector<double>, int, var, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_94;
typedef std::tuple<
    double, std::vector<double>, int, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_95;
typedef std::tuple<double, std::vector<double>, int, std::vector<var>, double,
                   empty>
    type_vv_real_real_int_real_real_96;
typedef std::tuple<double, std::vector<double>, int, std::vector<var>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_97;
typedef std::tuple<double, std::vector<double>, int, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_98;
typedef std::tuple<double, std::vector<double>, int, std::vector<var>, var,
                   empty>
    type_vv_real_real_int_real_real_99;
typedef std::tuple<double, std::vector<double>, int, std::vector<var>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_100;
typedef std::tuple<
    double, std::vector<double>, int, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_101;
typedef std::tuple<
    double, std::vector<double>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_102;
typedef std::tuple<
    double, std::vector<double>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_103;
typedef std::tuple<
    double, std::vector<double>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_104;
typedef std::tuple<
    double, std::vector<double>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_105;
typedef std::tuple<
    double, std::vector<double>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_106;
typedef std::tuple<
    double, std::vector<double>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_107;
typedef std::tuple<double, std::vector<double>, std::vector<int>, double, var,
                   empty>
    type_vv_real_real_int_real_real_108;
typedef std::tuple<double, std::vector<double>, std::vector<int>, double,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_109;
typedef std::tuple<
    double, std::vector<double>, std::vector<int>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_110;
typedef std::tuple<double, std::vector<double>, std::vector<int>,
                   std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_111;
typedef std::tuple<double, std::vector<double>, std::vector<int>,
                   std::vector<double>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_112;
typedef std::tuple<
    double, std::vector<double>, std::vector<int>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_113;
typedef std::tuple<double, std::vector<double>, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_114;
typedef std::tuple<double, std::vector<double>, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_115;
typedef std::tuple<
    double, std::vector<double>, std::vector<int>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_116;
typedef std::tuple<double, std::vector<double>, std::vector<int>, var, double,
                   empty>
    type_vv_real_real_int_real_real_117;
typedef std::tuple<double, std::vector<double>, std::vector<int>, var,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_118;
typedef std::tuple<double, std::vector<double>, std::vector<int>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_119;
typedef std::tuple<double, std::vector<double>, std::vector<int>, var, var,
                   empty>
    type_vv_real_real_int_real_real_120;
typedef std::tuple<double, std::vector<double>, std::vector<int>, var,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_121;
typedef std::tuple<
    double, std::vector<double>, std::vector<int>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_122;
typedef std::tuple<double, std::vector<double>, std::vector<int>,
                   std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_123;
typedef std::tuple<double, std::vector<double>, std::vector<int>,
                   std::vector<var>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_124;
typedef std::tuple<double, std::vector<double>, std::vector<int>,
                   std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   empty>
    type_vv_real_real_int_real_real_125;
typedef std::tuple<double, std::vector<double>, std::vector<int>,
                   std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_126;
typedef std::tuple<double, std::vector<double>, std::vector<int>,
                   std::vector<var>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_127;
typedef std::tuple<
    double, std::vector<double>, std::vector<int>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_128;
typedef std::tuple<
    double, std::vector<double>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_129;
typedef std::tuple<
    double, std::vector<double>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_130;
typedef std::tuple<
    double, std::vector<double>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_131;
typedef std::tuple<
    double, std::vector<double>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_132;
typedef std::tuple<
    double, std::vector<double>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_133;
typedef std::tuple<
    double, std::vector<double>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_134;
typedef std::tuple<double, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double, var, empty>
    type_vv_real_real_int_real_real_135;
typedef std::tuple<double, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_136;
typedef std::tuple<
    double, std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_137;
typedef std::tuple<double, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   var, empty>
    type_vv_real_real_int_real_real_138;
typedef std::tuple<double, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_139;
typedef std::tuple<
    double, std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_140;
typedef std::tuple<double, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_141;
typedef std::tuple<
    double, std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_142;
typedef std::tuple<
    double, std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_143;
typedef std::tuple<double, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var, double, empty>
    type_vv_real_real_int_real_real_144;
typedef std::tuple<double, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_145;
typedef std::tuple<double, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_146;
typedef std::tuple<double, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var, var, empty>
    type_vv_real_real_int_real_real_147;
typedef std::tuple<double, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_148;
typedef std::tuple<
    double, std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_149;
typedef std::tuple<double, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   double, empty>
    type_vv_real_real_int_real_real_150;
typedef std::tuple<double, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_151;
typedef std::tuple<double, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_152;
typedef std::tuple<double, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>, var,
                   empty>
    type_vv_real_real_int_real_real_153;
typedef std::tuple<double, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_154;
typedef std::tuple<
    double, std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_155;
typedef std::tuple<
    double, std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_156;
typedef std::tuple<
    double, std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_157;
typedef std::tuple<
    double, std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_158;
typedef std::tuple<
    double, std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_159;
typedef std::tuple<
    double, std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_160;
typedef std::tuple<
    double, std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_161;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
                   double, var, empty>
    type_vv_real_real_int_real_real_162;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
                   double, std::vector<var>, empty>
    type_vv_real_real_int_real_real_163;
typedef std::tuple<
    double, Eigen::Matrix<double, Eigen::Dynamic, 1>, int, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_164;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
                   std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_165;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
                   std::vector<double>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_166;
typedef std::tuple<
    double, Eigen::Matrix<double, Eigen::Dynamic, 1>, int, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_167;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_168;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_169;
typedef std::tuple<
    double, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_170;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>, int, var,
                   double, empty>
    type_vv_real_real_int_real_real_171;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>, int, var,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_172;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>, int, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_173;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>, int, var,
                   var, empty>
    type_vv_real_real_int_real_real_174;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>, int, var,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_175;
typedef std::tuple<
    double, Eigen::Matrix<double, Eigen::Dynamic, 1>, int, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_176;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
                   std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_177;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
                   std::vector<var>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_178;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
                   std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   empty>
    type_vv_real_real_int_real_real_179;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
                   std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_180;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
                   std::vector<var>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_181;
typedef std::tuple<
    double, Eigen::Matrix<double, Eigen::Dynamic, 1>, int, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_182;
typedef std::tuple<
    double, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_183;
typedef std::tuple<
    double, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_184;
typedef std::tuple<
    double, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_185;
typedef std::tuple<
    double, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_186;
typedef std::tuple<
    double, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_187;
typedef std::tuple<
    double, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_188;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, double, var, empty>
    type_vv_real_real_int_real_real_189;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, double, std::vector<var>, empty>
    type_vv_real_real_int_real_real_190;
typedef std::tuple<
    double, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_191;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_192;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, std::vector<double>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_193;
typedef std::tuple<
    double, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_194;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   var, empty>
    type_vv_real_real_int_real_real_195;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_196;
typedef std::tuple<
    double, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_197;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, var, double, empty>
    type_vv_real_real_int_real_real_198;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, var, std::vector<double>, empty>
    type_vv_real_real_int_real_real_199;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_200;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, var, var, empty>
    type_vv_real_real_int_real_real_201;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, var, std::vector<var>, empty>
    type_vv_real_real_int_real_real_202;
typedef std::tuple<
    double, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_203;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_204;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, std::vector<var>, std::vector<double>,
                   empty>
    type_vv_real_real_int_real_real_205;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_206;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_207;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, std::vector<var>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_208;
typedef std::tuple<
    double, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_209;
typedef std::tuple<
    double, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_210;
typedef std::tuple<
    double, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_211;
typedef std::tuple<
    double, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_212;
typedef std::tuple<
    double, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_213;
typedef std::tuple<
    double, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_214;
typedef std::tuple<
    double, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_215;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double, var, empty>
    type_vv_real_real_int_real_real_216;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_217;
typedef std::tuple<
    double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_218;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   var, empty>
    type_vv_real_real_int_real_real_219;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_220;
typedef std::tuple<
    double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_221;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_222;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_223;
typedef std::tuple<
    double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_224;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var, double, empty>
    type_vv_real_real_int_real_real_225;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_226;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_227;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var, var, empty>
    type_vv_real_real_int_real_real_228;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_229;
typedef std::tuple<
    double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_230;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   double, empty>
    type_vv_real_real_int_real_real_231;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_232;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_233;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>, var,
                   empty>
    type_vv_real_real_int_real_real_234;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_235;
typedef std::tuple<
    double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_236;
typedef std::tuple<
    double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_237;
typedef std::tuple<
    double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_238;
typedef std::tuple<
    double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_239;
typedef std::tuple<
    double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_240;
typedef std::tuple<
    double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_241;
typedef std::tuple<
    double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_242;
typedef std::tuple<double, var, int, double, double, empty>
    type_vv_real_real_int_real_real_243;
typedef std::tuple<double, var, int, double, std::vector<double>, empty>
    type_vv_real_real_int_real_real_244;
typedef std::tuple<double, var, int, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_245;
typedef std::tuple<double, var, int, double, var, empty>
    type_vv_real_real_int_real_real_246;
typedef std::tuple<double, var, int, double, std::vector<var>, empty>
    type_vv_real_real_int_real_real_247;
typedef std::tuple<
    double, var, int, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_248;
typedef std::tuple<double, var, int, std::vector<double>, double, empty>
    type_vv_real_real_int_real_real_249;
typedef std::tuple<double, var, int, std::vector<double>, std::vector<double>,
                   empty>
    type_vv_real_real_int_real_real_250;
typedef std::tuple<double, var, int, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_251;
typedef std::tuple<double, var, int, std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_252;
typedef std::tuple<double, var, int, std::vector<double>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_253;
typedef std::tuple<
    double, var, int, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_254;
typedef std::tuple<double, var, int, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   double, empty>
    type_vv_real_real_int_real_real_255;
typedef std::tuple<double, var, int, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_256;
typedef std::tuple<double, var, int, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_257;
typedef std::tuple<double, var, int, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   var, empty>
    type_vv_real_real_int_real_real_258;
typedef std::tuple<double, var, int, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_259;
typedef std::tuple<
    double, var, int, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_260;
typedef std::tuple<double, var, int, var, double, empty>
    type_vv_real_real_int_real_real_261;
typedef std::tuple<double, var, int, var, std::vector<double>, empty>
    type_vv_real_real_int_real_real_262;
typedef std::tuple<double, var, int, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_263;
typedef std::tuple<double, var, int, var, var, empty>
    type_vv_real_real_int_real_real_264;
typedef std::tuple<double, var, int, var, std::vector<var>, empty>
    type_vv_real_real_int_real_real_265;
typedef std::tuple<
    double, var, int, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_266;
typedef std::tuple<double, var, int, std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_267;
typedef std::tuple<double, var, int, std::vector<var>, std::vector<double>,
                   empty>
    type_vv_real_real_int_real_real_268;
typedef std::tuple<double, var, int, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_269;
typedef std::tuple<double, var, int, std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_270;
typedef std::tuple<double, var, int, std::vector<var>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_271;
typedef std::tuple<
    double, var, int, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_272;
typedef std::tuple<
    double, var, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_273;
typedef std::tuple<
    double, var, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_274;
typedef std::tuple<
    double, var, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_275;
typedef std::tuple<
    double, var, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_276;
typedef std::tuple<
    double, var, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_277;
typedef std::tuple<
    double, var, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_278;
typedef std::tuple<double, var, std::vector<int>, double, double, empty>
    type_vv_real_real_int_real_real_279;
typedef std::tuple<double, var, std::vector<int>, double, std::vector<double>,
                   empty>
    type_vv_real_real_int_real_real_280;
typedef std::tuple<double, var, std::vector<int>, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_281;
typedef std::tuple<double, var, std::vector<int>, double, var, empty>
    type_vv_real_real_int_real_real_282;
typedef std::tuple<double, var, std::vector<int>, double, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_283;
typedef std::tuple<
    double, var, std::vector<int>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_284;
typedef std::tuple<double, var, std::vector<int>, std::vector<double>, double,
                   empty>
    type_vv_real_real_int_real_real_285;
typedef std::tuple<double, var, std::vector<int>, std::vector<double>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_286;
typedef std::tuple<double, var, std::vector<int>, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_287;
typedef std::tuple<double, var, std::vector<int>, std::vector<double>, var,
                   empty>
    type_vv_real_real_int_real_real_288;
typedef std::tuple<double, var, std::vector<int>, std::vector<double>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_289;
typedef std::tuple<
    double, var, std::vector<int>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_290;
typedef std::tuple<double, var, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_291;
typedef std::tuple<double, var, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_292;
typedef std::tuple<double, var, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_293;
typedef std::tuple<double, var, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_294;
typedef std::tuple<double, var, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_295;
typedef std::tuple<
    double, var, std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_296;
typedef std::tuple<double, var, std::vector<int>, var, double, empty>
    type_vv_real_real_int_real_real_297;
typedef std::tuple<double, var, std::vector<int>, var, std::vector<double>,
                   empty>
    type_vv_real_real_int_real_real_298;
typedef std::tuple<double, var, std::vector<int>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_299;
typedef std::tuple<double, var, std::vector<int>, var, var, empty>
    type_vv_real_real_int_real_real_300;
typedef std::tuple<double, var, std::vector<int>, var, std::vector<var>, empty>
    type_vv_real_real_int_real_real_301;
typedef std::tuple<
    double, var, std::vector<int>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_302;
typedef std::tuple<double, var, std::vector<int>, std::vector<var>, double,
                   empty>
    type_vv_real_real_int_real_real_303;
typedef std::tuple<double, var, std::vector<int>, std::vector<var>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_304;
typedef std::tuple<double, var, std::vector<int>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_305;
typedef std::tuple<double, var, std::vector<int>, std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_306;
typedef std::tuple<double, var, std::vector<int>, std::vector<var>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_307;
typedef std::tuple<
    double, var, std::vector<int>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_308;
typedef std::tuple<
    double, var, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_309;
typedef std::tuple<
    double, var, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_310;
typedef std::tuple<
    double, var, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_311;
typedef std::tuple<
    double, var, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_312;
typedef std::tuple<
    double, var, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_313;
typedef std::tuple<
    double, var, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_314;
typedef std::tuple<double, var, Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
                   double, empty>
    type_vv_real_real_int_real_real_315;
typedef std::tuple<double, var, Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_316;
typedef std::tuple<double, var, Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_317;
typedef std::tuple<double, var, Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
                   var, empty>
    type_vv_real_real_int_real_real_318;
typedef std::tuple<double, var, Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_319;
typedef std::tuple<
    double, var, Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_320;
typedef std::tuple<double, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<double>, double, empty>
    type_vv_real_real_int_real_real_321;
typedef std::tuple<double, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<double>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_322;
typedef std::tuple<double, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_323;
typedef std::tuple<double, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_324;
typedef std::tuple<double, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<double>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_325;
typedef std::tuple<
    double, var, Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_326;
typedef std::tuple<double, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_327;
typedef std::tuple<double, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_328;
typedef std::tuple<double, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_329;
typedef std::tuple<double, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_330;
typedef std::tuple<double, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_331;
typedef std::tuple<
    double, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_332;
typedef std::tuple<double, var, Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
                   double, empty>
    type_vv_real_real_int_real_real_333;
typedef std::tuple<double, var, Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_334;
typedef std::tuple<double, var, Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_335;
typedef std::tuple<double, var, Eigen::Matrix<int, Eigen::Dynamic, 1>, var, var,
                   empty>
    type_vv_real_real_int_real_real_336;
typedef std::tuple<double, var, Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_337;
typedef std::tuple<
    double, var, Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_338;
typedef std::tuple<double, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_339;
typedef std::tuple<double, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<var>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_340;
typedef std::tuple<double, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   empty>
    type_vv_real_real_int_real_real_341;
typedef std::tuple<double, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_342;
typedef std::tuple<double, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<var>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_343;
typedef std::tuple<
    double, var, Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_344;
typedef std::tuple<
    double, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_345;
typedef std::tuple<
    double, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_346;
typedef std::tuple<
    double, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_347;
typedef std::tuple<
    double, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_348;
typedef std::tuple<
    double, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_349;
typedef std::tuple<
    double, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_350;
typedef std::tuple<double, std::vector<var>, int, double, double, empty>
    type_vv_real_real_int_real_real_351;
typedef std::tuple<double, std::vector<var>, int, double, std::vector<double>,
                   empty>
    type_vv_real_real_int_real_real_352;
typedef std::tuple<double, std::vector<var>, int, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_353;
typedef std::tuple<double, std::vector<var>, int, double, var, empty>
    type_vv_real_real_int_real_real_354;
typedef std::tuple<double, std::vector<var>, int, double, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_355;
typedef std::tuple<
    double, std::vector<var>, int, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_356;
typedef std::tuple<double, std::vector<var>, int, std::vector<double>, double,
                   empty>
    type_vv_real_real_int_real_real_357;
typedef std::tuple<double, std::vector<var>, int, std::vector<double>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_358;
typedef std::tuple<double, std::vector<var>, int, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_359;
typedef std::tuple<double, std::vector<var>, int, std::vector<double>, var,
                   empty>
    type_vv_real_real_int_real_real_360;
typedef std::tuple<double, std::vector<var>, int, std::vector<double>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_361;
typedef std::tuple<
    double, std::vector<var>, int, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_362;
typedef std::tuple<double, std::vector<var>, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_363;
typedef std::tuple<double, std::vector<var>, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_364;
typedef std::tuple<double, std::vector<var>, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_365;
typedef std::tuple<double, std::vector<var>, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_366;
typedef std::tuple<double, std::vector<var>, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_367;
typedef std::tuple<
    double, std::vector<var>, int, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_368;
typedef std::tuple<double, std::vector<var>, int, var, double, empty>
    type_vv_real_real_int_real_real_369;
typedef std::tuple<double, std::vector<var>, int, var, std::vector<double>,
                   empty>
    type_vv_real_real_int_real_real_370;
typedef std::tuple<double, std::vector<var>, int, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_371;
typedef std::tuple<double, std::vector<var>, int, var, var, empty>
    type_vv_real_real_int_real_real_372;
typedef std::tuple<double, std::vector<var>, int, var, std::vector<var>, empty>
    type_vv_real_real_int_real_real_373;
typedef std::tuple<
    double, std::vector<var>, int, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_374;
typedef std::tuple<double, std::vector<var>, int, std::vector<var>, double,
                   empty>
    type_vv_real_real_int_real_real_375;
typedef std::tuple<double, std::vector<var>, int, std::vector<var>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_376;
typedef std::tuple<double, std::vector<var>, int, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_377;
typedef std::tuple<double, std::vector<var>, int, std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_378;
typedef std::tuple<double, std::vector<var>, int, std::vector<var>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_379;
typedef std::tuple<
    double, std::vector<var>, int, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_380;
typedef std::tuple<
    double, std::vector<var>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_381;
typedef std::tuple<
    double, std::vector<var>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_382;
typedef std::tuple<
    double, std::vector<var>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_383;
typedef std::tuple<
    double, std::vector<var>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_384;
typedef std::tuple<
    double, std::vector<var>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_385;
typedef std::tuple<
    double, std::vector<var>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_386;
typedef std::tuple<double, std::vector<var>, std::vector<int>, double, double,
                   empty>
    type_vv_real_real_int_real_real_387;
typedef std::tuple<double, std::vector<var>, std::vector<int>, double,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_388;
typedef std::tuple<double, std::vector<var>, std::vector<int>, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_389;
typedef std::tuple<double, std::vector<var>, std::vector<int>, double, var,
                   empty>
    type_vv_real_real_int_real_real_390;
typedef std::tuple<double, std::vector<var>, std::vector<int>, double,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_391;
typedef std::tuple<
    double, std::vector<var>, std::vector<int>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_392;
typedef std::tuple<double, std::vector<var>, std::vector<int>,
                   std::vector<double>, double, empty>
    type_vv_real_real_int_real_real_393;
typedef std::tuple<double, std::vector<var>, std::vector<int>,
                   std::vector<double>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_394;
typedef std::tuple<double, std::vector<var>, std::vector<int>,
                   std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_395;
typedef std::tuple<double, std::vector<var>, std::vector<int>,
                   std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_396;
typedef std::tuple<double, std::vector<var>, std::vector<int>,
                   std::vector<double>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_397;
typedef std::tuple<
    double, std::vector<var>, std::vector<int>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_398;
typedef std::tuple<double, std::vector<var>, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_399;
typedef std::tuple<double, std::vector<var>, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_400;
typedef std::tuple<double, std::vector<var>, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_401;
typedef std::tuple<double, std::vector<var>, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_402;
typedef std::tuple<double, std::vector<var>, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_403;
typedef std::tuple<
    double, std::vector<var>, std::vector<int>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_404;
typedef std::tuple<double, std::vector<var>, std::vector<int>, var, double,
                   empty>
    type_vv_real_real_int_real_real_405;
typedef std::tuple<double, std::vector<var>, std::vector<int>, var,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_406;
typedef std::tuple<double, std::vector<var>, std::vector<int>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_407;
typedef std::tuple<double, std::vector<var>, std::vector<int>, var, var, empty>
    type_vv_real_real_int_real_real_408;
typedef std::tuple<double, std::vector<var>, std::vector<int>, var,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_409;
typedef std::tuple<
    double, std::vector<var>, std::vector<int>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_410;
typedef std::tuple<double, std::vector<var>, std::vector<int>, std::vector<var>,
                   double, empty>
    type_vv_real_real_int_real_real_411;
typedef std::tuple<double, std::vector<var>, std::vector<int>, std::vector<var>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_412;
typedef std::tuple<double, std::vector<var>, std::vector<int>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_413;
typedef std::tuple<double, std::vector<var>, std::vector<int>, std::vector<var>,
                   var, empty>
    type_vv_real_real_int_real_real_414;
typedef std::tuple<double, std::vector<var>, std::vector<int>, std::vector<var>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_415;
typedef std::tuple<
    double, std::vector<var>, std::vector<int>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_416;
typedef std::tuple<
    double, std::vector<var>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_417;
typedef std::tuple<
    double, std::vector<var>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_418;
typedef std::tuple<
    double, std::vector<var>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_419;
typedef std::tuple<
    double, std::vector<var>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_420;
typedef std::tuple<
    double, std::vector<var>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_421;
typedef std::tuple<
    double, std::vector<var>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_422;
typedef std::tuple<double, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double, double, empty>
    type_vv_real_real_int_real_real_423;
typedef std::tuple<double, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_424;
typedef std::tuple<double, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_425;
typedef std::tuple<double, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double, var, empty>
    type_vv_real_real_int_real_real_426;
typedef std::tuple<double, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_427;
typedef std::tuple<
    double, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_428;
typedef std::tuple<double, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   double, empty>
    type_vv_real_real_int_real_real_429;
typedef std::tuple<double, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_430;
typedef std::tuple<double, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_431;
typedef std::tuple<double, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   var, empty>
    type_vv_real_real_int_real_real_432;
typedef std::tuple<double, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_433;
typedef std::tuple<
    double, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_434;
typedef std::tuple<double, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_435;
typedef std::tuple<
    double, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_436;
typedef std::tuple<double, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_437;
typedef std::tuple<double, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_438;
typedef std::tuple<
    double, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_439;
typedef std::tuple<
    double, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_440;
typedef std::tuple<double, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var, double, empty>
    type_vv_real_real_int_real_real_441;
typedef std::tuple<double, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_442;
typedef std::tuple<double, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_443;
typedef std::tuple<double, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var, var, empty>
    type_vv_real_real_int_real_real_444;
typedef std::tuple<double, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_445;
typedef std::tuple<
    double, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_446;
typedef std::tuple<double, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   double, empty>
    type_vv_real_real_int_real_real_447;
typedef std::tuple<double, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_448;
typedef std::tuple<double, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_449;
typedef std::tuple<double, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>, var,
                   empty>
    type_vv_real_real_int_real_real_450;
typedef std::tuple<double, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_451;
typedef std::tuple<
    double, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_452;
typedef std::tuple<
    double, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_453;
typedef std::tuple<
    double, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_454;
typedef std::tuple<
    double, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_455;
typedef std::tuple<
    double, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_456;
typedef std::tuple<
    double, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_457;
typedef std::tuple<
    double, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_458;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    int, double, double, empty>
    type_vv_real_real_int_real_real_459;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    int, double, std::vector<double>, empty>
    type_vv_real_real_int_real_real_460;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    int, double, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_461;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    int, double, var, empty>
    type_vv_real_real_int_real_real_462;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    int, double, std::vector<var>, empty>
    type_vv_real_real_int_real_real_463;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    int, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_464;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    int, std::vector<double>, double, empty>
    type_vv_real_real_int_real_real_465;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    int, std::vector<double>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_466;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    int, std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_467;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    int, std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_468;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    int, std::vector<double>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_469;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    int, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_470;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    int, Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_471;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    int, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_472;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    int, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_473;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    int, Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_474;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    int, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_475;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    int, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_476;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    int, var, double, empty>
    type_vv_real_real_int_real_real_477;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    int, var, std::vector<double>, empty>
    type_vv_real_real_int_real_real_478;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    int, var, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_479;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    int, var, var, empty>
    type_vv_real_real_int_real_real_480;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    int, var, std::vector<var>, empty>
    type_vv_real_real_int_real_real_481;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    int, var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    empty>
    type_vv_real_real_int_real_real_482;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    int, std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_483;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    int, std::vector<var>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_484;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    int, std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_485;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    int, std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_486;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    int, std::vector<var>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_487;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    int, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_488;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    int, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    double, empty>
    type_vv_real_real_int_real_real_489;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    int, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_490;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    int, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_491;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    int, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    empty>
    type_vv_real_real_int_real_real_492;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    int, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_493;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    int, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_494;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, double, double, empty>
    type_vv_real_real_int_real_real_495;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, double, std::vector<double>, empty>
    type_vv_real_real_int_real_real_496;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, double, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_497;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, double, var, empty>
    type_vv_real_real_int_real_real_498;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, double, std::vector<var>, empty>
    type_vv_real_real_int_real_real_499;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_500;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<double>, double, empty>
    type_vv_real_real_int_real_real_501;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<double>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_502;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<double>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_503;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_504;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<double>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_505;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_506;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_507;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_508;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_509;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_510;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_511;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_512;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, var, double, empty>
    type_vv_real_real_int_real_real_513;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, var, std::vector<double>, empty>
    type_vv_real_real_int_real_real_514;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, var, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_515;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, var, var, empty>
    type_vv_real_real_int_real_real_516;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, var, std::vector<var>, empty>
    type_vv_real_real_int_real_real_517;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_518;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_519;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<var>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_520;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<var>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_521;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_522;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<var>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_523;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_524;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_525;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_526;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_527;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_528;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_529;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_530;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double, double, empty>
    type_vv_real_real_int_real_real_531;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double, std::vector<double>, empty>
    type_vv_real_real_int_real_real_532;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_533;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double, var, empty>
    type_vv_real_real_int_real_real_534;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double, std::vector<var>, empty>
    type_vv_real_real_int_real_real_535;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_536;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>, double, empty>
    type_vv_real_real_int_real_real_537;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_538;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_539;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_540;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_541;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_542;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_543;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_544;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_545;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_546;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_547;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_548;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var, double, empty>
    type_vv_real_real_int_real_real_549;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var, std::vector<double>, empty>
    type_vv_real_real_int_real_real_550;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_551;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var, var, empty>
    type_vv_real_real_int_real_real_552;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var, std::vector<var>, empty>
    type_vv_real_real_int_real_real_553;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_554;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_555;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_556;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_557;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_558;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>, std::vector<var>,
    empty>
    type_vv_real_real_int_real_real_559;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_560;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_561;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_562;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_563;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_564;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_565;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_566;
typedef std::tuple<std::vector<double>, double, int, double, var, empty>
    type_vv_real_real_int_real_real_567;
typedef std::tuple<std::vector<double>, double, int, double, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_568;
typedef std::tuple<
    std::vector<double>, double, int, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_569;
typedef std::tuple<std::vector<double>, double, int, std::vector<double>, var,
                   empty>
    type_vv_real_real_int_real_real_570;
typedef std::tuple<std::vector<double>, double, int, std::vector<double>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_571;
typedef std::tuple<
    std::vector<double>, double, int, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_572;
typedef std::tuple<std::vector<double>, double, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_573;
typedef std::tuple<std::vector<double>, double, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_574;
typedef std::tuple<
    std::vector<double>, double, int, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_575;
typedef std::tuple<std::vector<double>, double, int, var, double, empty>
    type_vv_real_real_int_real_real_576;
typedef std::tuple<std::vector<double>, double, int, var, std::vector<double>,
                   empty>
    type_vv_real_real_int_real_real_577;
typedef std::tuple<std::vector<double>, double, int, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_578;
typedef std::tuple<std::vector<double>, double, int, var, var, empty>
    type_vv_real_real_int_real_real_579;
typedef std::tuple<std::vector<double>, double, int, var, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_580;
typedef std::tuple<
    std::vector<double>, double, int, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_581;
typedef std::tuple<std::vector<double>, double, int, std::vector<var>, double,
                   empty>
    type_vv_real_real_int_real_real_582;
typedef std::tuple<std::vector<double>, double, int, std::vector<var>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_583;
typedef std::tuple<std::vector<double>, double, int, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_584;
typedef std::tuple<std::vector<double>, double, int, std::vector<var>, var,
                   empty>
    type_vv_real_real_int_real_real_585;
typedef std::tuple<std::vector<double>, double, int, std::vector<var>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_586;
typedef std::tuple<
    std::vector<double>, double, int, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_587;
typedef std::tuple<
    std::vector<double>, double, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_588;
typedef std::tuple<
    std::vector<double>, double, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_589;
typedef std::tuple<
    std::vector<double>, double, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_590;
typedef std::tuple<
    std::vector<double>, double, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_591;
typedef std::tuple<
    std::vector<double>, double, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_592;
typedef std::tuple<
    std::vector<double>, double, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_593;
typedef std::tuple<std::vector<double>, double, std::vector<int>, double, var,
                   empty>
    type_vv_real_real_int_real_real_594;
typedef std::tuple<std::vector<double>, double, std::vector<int>, double,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_595;
typedef std::tuple<
    std::vector<double>, double, std::vector<int>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_596;
typedef std::tuple<std::vector<double>, double, std::vector<int>,
                   std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_597;
typedef std::tuple<std::vector<double>, double, std::vector<int>,
                   std::vector<double>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_598;
typedef std::tuple<
    std::vector<double>, double, std::vector<int>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_599;
typedef std::tuple<std::vector<double>, double, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_600;
typedef std::tuple<std::vector<double>, double, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_601;
typedef std::tuple<
    std::vector<double>, double, std::vector<int>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_602;
typedef std::tuple<std::vector<double>, double, std::vector<int>, var, double,
                   empty>
    type_vv_real_real_int_real_real_603;
typedef std::tuple<std::vector<double>, double, std::vector<int>, var,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_604;
typedef std::tuple<std::vector<double>, double, std::vector<int>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_605;
typedef std::tuple<std::vector<double>, double, std::vector<int>, var, var,
                   empty>
    type_vv_real_real_int_real_real_606;
typedef std::tuple<std::vector<double>, double, std::vector<int>, var,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_607;
typedef std::tuple<
    std::vector<double>, double, std::vector<int>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_608;
typedef std::tuple<std::vector<double>, double, std::vector<int>,
                   std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_609;
typedef std::tuple<std::vector<double>, double, std::vector<int>,
                   std::vector<var>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_610;
typedef std::tuple<std::vector<double>, double, std::vector<int>,
                   std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   empty>
    type_vv_real_real_int_real_real_611;
typedef std::tuple<std::vector<double>, double, std::vector<int>,
                   std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_612;
typedef std::tuple<std::vector<double>, double, std::vector<int>,
                   std::vector<var>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_613;
typedef std::tuple<
    std::vector<double>, double, std::vector<int>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_614;
typedef std::tuple<
    std::vector<double>, double, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_615;
typedef std::tuple<
    std::vector<double>, double, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_616;
typedef std::tuple<
    std::vector<double>, double, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_617;
typedef std::tuple<
    std::vector<double>, double, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_618;
typedef std::tuple<
    std::vector<double>, double, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_619;
typedef std::tuple<
    std::vector<double>, double, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_620;
typedef std::tuple<std::vector<double>, double,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double, var, empty>
    type_vv_real_real_int_real_real_621;
typedef std::tuple<std::vector<double>, double,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_622;
typedef std::tuple<
    std::vector<double>, double, Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_623;
typedef std::tuple<std::vector<double>, double,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   var, empty>
    type_vv_real_real_int_real_real_624;
typedef std::tuple<std::vector<double>, double,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_625;
typedef std::tuple<
    std::vector<double>, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_626;
typedef std::tuple<std::vector<double>, double,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_627;
typedef std::tuple<
    std::vector<double>, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_628;
typedef std::tuple<
    std::vector<double>, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_629;
typedef std::tuple<std::vector<double>, double,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var, double, empty>
    type_vv_real_real_int_real_real_630;
typedef std::tuple<std::vector<double>, double,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_631;
typedef std::tuple<std::vector<double>, double,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_632;
typedef std::tuple<std::vector<double>, double,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var, var, empty>
    type_vv_real_real_int_real_real_633;
typedef std::tuple<std::vector<double>, double,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_634;
typedef std::tuple<
    std::vector<double>, double, Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_635;
typedef std::tuple<std::vector<double>, double,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   double, empty>
    type_vv_real_real_int_real_real_636;
typedef std::tuple<std::vector<double>, double,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_637;
typedef std::tuple<std::vector<double>, double,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_638;
typedef std::tuple<std::vector<double>, double,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>, var,
                   empty>
    type_vv_real_real_int_real_real_639;
typedef std::tuple<std::vector<double>, double,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_640;
typedef std::tuple<
    std::vector<double>, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_641;
typedef std::tuple<
    std::vector<double>, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_642;
typedef std::tuple<
    std::vector<double>, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_643;
typedef std::tuple<
    std::vector<double>, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_644;
typedef std::tuple<
    std::vector<double>, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_645;
typedef std::tuple<
    std::vector<double>, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_646;
typedef std::tuple<
    std::vector<double>, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_647;
typedef std::tuple<std::vector<double>, std::vector<double>, int, double, var,
                   empty>
    type_vv_real_real_int_real_real_648;
typedef std::tuple<std::vector<double>, std::vector<double>, int, double,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_649;
typedef std::tuple<
    std::vector<double>, std::vector<double>, int, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_650;
typedef std::tuple<std::vector<double>, std::vector<double>, int,
                   std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_651;
typedef std::tuple<std::vector<double>, std::vector<double>, int,
                   std::vector<double>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_652;
typedef std::tuple<
    std::vector<double>, std::vector<double>, int, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_653;
typedef std::tuple<std::vector<double>, std::vector<double>, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_654;
typedef std::tuple<std::vector<double>, std::vector<double>, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_655;
typedef std::tuple<
    std::vector<double>, std::vector<double>, int,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_656;
typedef std::tuple<std::vector<double>, std::vector<double>, int, var, double,
                   empty>
    type_vv_real_real_int_real_real_657;
typedef std::tuple<std::vector<double>, std::vector<double>, int, var,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_658;
typedef std::tuple<std::vector<double>, std::vector<double>, int, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_659;
typedef std::tuple<std::vector<double>, std::vector<double>, int, var, var,
                   empty>
    type_vv_real_real_int_real_real_660;
typedef std::tuple<std::vector<double>, std::vector<double>, int, var,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_661;
typedef std::tuple<
    std::vector<double>, std::vector<double>, int, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_662;
typedef std::tuple<std::vector<double>, std::vector<double>, int,
                   std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_663;
typedef std::tuple<std::vector<double>, std::vector<double>, int,
                   std::vector<var>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_664;
typedef std::tuple<std::vector<double>, std::vector<double>, int,
                   std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   empty>
    type_vv_real_real_int_real_real_665;
typedef std::tuple<std::vector<double>, std::vector<double>, int,
                   std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_666;
typedef std::tuple<std::vector<double>, std::vector<double>, int,
                   std::vector<var>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_667;
typedef std::tuple<
    std::vector<double>, std::vector<double>, int, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_668;
typedef std::tuple<
    std::vector<double>, std::vector<double>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_669;
typedef std::tuple<
    std::vector<double>, std::vector<double>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_670;
typedef std::tuple<
    std::vector<double>, std::vector<double>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_671;
typedef std::tuple<
    std::vector<double>, std::vector<double>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_672;
typedef std::tuple<
    std::vector<double>, std::vector<double>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_673;
typedef std::tuple<
    std::vector<double>, std::vector<double>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_674;
typedef std::tuple<std::vector<double>, std::vector<double>, std::vector<int>,
                   double, var, empty>
    type_vv_real_real_int_real_real_675;
typedef std::tuple<std::vector<double>, std::vector<double>, std::vector<int>,
                   double, std::vector<var>, empty>
    type_vv_real_real_int_real_real_676;
typedef std::tuple<
    std::vector<double>, std::vector<double>, std::vector<int>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_677;
typedef std::tuple<std::vector<double>, std::vector<double>, std::vector<int>,
                   std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_678;
typedef std::tuple<std::vector<double>, std::vector<double>, std::vector<int>,
                   std::vector<double>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_679;
typedef std::tuple<
    std::vector<double>, std::vector<double>, std::vector<int>,
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_680;
typedef std::tuple<std::vector<double>, std::vector<double>, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_681;
typedef std::tuple<std::vector<double>, std::vector<double>, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_682;
typedef std::tuple<
    std::vector<double>, std::vector<double>, std::vector<int>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_683;
typedef std::tuple<std::vector<double>, std::vector<double>, std::vector<int>,
                   var, double, empty>
    type_vv_real_real_int_real_real_684;
typedef std::tuple<std::vector<double>, std::vector<double>, std::vector<int>,
                   var, std::vector<double>, empty>
    type_vv_real_real_int_real_real_685;
typedef std::tuple<std::vector<double>, std::vector<double>, std::vector<int>,
                   var, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_686;
typedef std::tuple<std::vector<double>, std::vector<double>, std::vector<int>,
                   var, var, empty>
    type_vv_real_real_int_real_real_687;
typedef std::tuple<std::vector<double>, std::vector<double>, std::vector<int>,
                   var, std::vector<var>, empty>
    type_vv_real_real_int_real_real_688;
typedef std::tuple<
    std::vector<double>, std::vector<double>, std::vector<int>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_689;
typedef std::tuple<std::vector<double>, std::vector<double>, std::vector<int>,
                   std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_690;
typedef std::tuple<std::vector<double>, std::vector<double>, std::vector<int>,
                   std::vector<var>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_691;
typedef std::tuple<std::vector<double>, std::vector<double>, std::vector<int>,
                   std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   empty>
    type_vv_real_real_int_real_real_692;
typedef std::tuple<std::vector<double>, std::vector<double>, std::vector<int>,
                   std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_693;
typedef std::tuple<std::vector<double>, std::vector<double>, std::vector<int>,
                   std::vector<var>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_694;
typedef std::tuple<
    std::vector<double>, std::vector<double>, std::vector<int>,
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_695;
typedef std::tuple<
    std::vector<double>, std::vector<double>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_696;
typedef std::tuple<
    std::vector<double>, std::vector<double>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_697;
typedef std::tuple<
    std::vector<double>, std::vector<double>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_698;
typedef std::tuple<
    std::vector<double>, std::vector<double>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_699;
typedef std::tuple<
    std::vector<double>, std::vector<double>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_700;
typedef std::tuple<
    std::vector<double>, std::vector<double>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_701;
typedef std::tuple<std::vector<double>, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double, var, empty>
    type_vv_real_real_int_real_real_702;
typedef std::tuple<std::vector<double>, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_703;
typedef std::tuple<
    std::vector<double>, std::vector<double>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_704;
typedef std::tuple<std::vector<double>, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   var, empty>
    type_vv_real_real_int_real_real_705;
typedef std::tuple<std::vector<double>, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_706;
typedef std::tuple<
    std::vector<double>, std::vector<double>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_707;
typedef std::tuple<std::vector<double>, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_708;
typedef std::tuple<std::vector<double>, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_709;
typedef std::tuple<
    std::vector<double>, std::vector<double>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_710;
typedef std::tuple<std::vector<double>, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var, double, empty>
    type_vv_real_real_int_real_real_711;
typedef std::tuple<std::vector<double>, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_712;
typedef std::tuple<std::vector<double>, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_713;
typedef std::tuple<std::vector<double>, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var, var, empty>
    type_vv_real_real_int_real_real_714;
typedef std::tuple<std::vector<double>, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_715;
typedef std::tuple<
    std::vector<double>, std::vector<double>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_716;
typedef std::tuple<std::vector<double>, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   double, empty>
    type_vv_real_real_int_real_real_717;
typedef std::tuple<std::vector<double>, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_718;
typedef std::tuple<std::vector<double>, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_719;
typedef std::tuple<std::vector<double>, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>, var,
                   empty>
    type_vv_real_real_int_real_real_720;
typedef std::tuple<std::vector<double>, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_721;
typedef std::tuple<
    std::vector<double>, std::vector<double>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_722;
typedef std::tuple<
    std::vector<double>, std::vector<double>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_723;
typedef std::tuple<
    std::vector<double>, std::vector<double>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_724;
typedef std::tuple<
    std::vector<double>, std::vector<double>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_725;
typedef std::tuple<
    std::vector<double>, std::vector<double>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_726;
typedef std::tuple<
    std::vector<double>, std::vector<double>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_727;
typedef std::tuple<
    std::vector<double>, std::vector<double>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_728;
typedef std::tuple<std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, int, double, var,
                   empty>
    type_vv_real_real_int_real_real_729;
typedef std::tuple<std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, int, double,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_730;
typedef std::tuple<
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, int, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_731;
typedef std::tuple<std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
                   std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_732;
typedef std::tuple<std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
                   std::vector<double>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_733;
typedef std::tuple<
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_734;
typedef std::tuple<std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_735;
typedef std::tuple<
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_736;
typedef std::tuple<
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_737;
typedef std::tuple<std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, int, var, double,
                   empty>
    type_vv_real_real_int_real_real_738;
typedef std::tuple<std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, int, var,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_739;
typedef std::tuple<std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, int, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_740;
typedef std::tuple<std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, int, var, var,
                   empty>
    type_vv_real_real_int_real_real_741;
typedef std::tuple<std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, int, var,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_742;
typedef std::tuple<
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, int, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_743;
typedef std::tuple<std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
                   std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_744;
typedef std::tuple<std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
                   std::vector<var>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_745;
typedef std::tuple<
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_746;
typedef std::tuple<std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
                   std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_747;
typedef std::tuple<std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
                   std::vector<var>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_748;
typedef std::tuple<
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_749;
typedef std::tuple<
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_750;
typedef std::tuple<
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_751;
typedef std::tuple<
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_752;
typedef std::tuple<
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_753;
typedef std::tuple<
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_754;
typedef std::tuple<
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_755;
typedef std::tuple<std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
                   double, var, empty>
    type_vv_real_real_int_real_real_756;
typedef std::tuple<std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
                   double, std::vector<var>, empty>
    type_vv_real_real_int_real_real_757;
typedef std::tuple<
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    std::vector<int>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_758;
typedef std::tuple<std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
                   std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_759;
typedef std::tuple<std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
                   std::vector<double>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_760;
typedef std::tuple<
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    std::vector<int>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_761;
typedef std::tuple<std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_762;
typedef std::tuple<std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_763;
typedef std::tuple<
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_764;
typedef std::tuple<std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
                   var, double, empty>
    type_vv_real_real_int_real_real_765;
typedef std::tuple<std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
                   var, std::vector<double>, empty>
    type_vv_real_real_int_real_real_766;
typedef std::tuple<std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
                   var, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_767;
typedef std::tuple<std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
                   var, var, empty>
    type_vv_real_real_int_real_real_768;
typedef std::tuple<std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
                   var, std::vector<var>, empty>
    type_vv_real_real_int_real_real_769;
typedef std::tuple<
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    std::vector<int>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_770;
typedef std::tuple<std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
                   std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_771;
typedef std::tuple<std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
                   std::vector<var>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_772;
typedef std::tuple<std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
                   std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   empty>
    type_vv_real_real_int_real_real_773;
typedef std::tuple<std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
                   std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_774;
typedef std::tuple<std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
                   std::vector<var>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_775;
typedef std::tuple<
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    std::vector<int>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_776;
typedef std::tuple<
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_777;
typedef std::tuple<
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_778;
typedef std::tuple<
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_779;
typedef std::tuple<
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_780;
typedef std::tuple<
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_781;
typedef std::tuple<
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_782;
typedef std::tuple<std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double, var, empty>
    type_vv_real_real_int_real_real_783;
typedef std::tuple<
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double, std::vector<var>, empty>
    type_vv_real_real_int_real_real_784;
typedef std::tuple<
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_785;
typedef std::tuple<
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_786;
typedef std::tuple<std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_787;
typedef std::tuple<
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_788;
typedef std::tuple<std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_789;
typedef std::tuple<
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_790;
typedef std::tuple<
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_791;
typedef std::tuple<std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var, double, empty>
    type_vv_real_real_int_real_real_792;
typedef std::tuple<
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var, std::vector<double>, empty>
    type_vv_real_real_int_real_real_793;
typedef std::tuple<std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_794;
typedef std::tuple<std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var, var, empty>
    type_vv_real_real_int_real_real_795;
typedef std::tuple<
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var, std::vector<var>, empty>
    type_vv_real_real_int_real_real_796;
typedef std::tuple<
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_797;
typedef std::tuple<
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_798;
typedef std::tuple<std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_799;
typedef std::tuple<std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_800;
typedef std::tuple<
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_801;
typedef std::tuple<std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_802;
typedef std::tuple<
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_803;
typedef std::tuple<
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_804;
typedef std::tuple<
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_805;
typedef std::tuple<
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_806;
typedef std::tuple<
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_807;
typedef std::tuple<
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_808;
typedef std::tuple<
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_809;
typedef std::tuple<std::vector<double>, var, int, double, double, empty>
    type_vv_real_real_int_real_real_810;
typedef std::tuple<std::vector<double>, var, int, double, std::vector<double>,
                   empty>
    type_vv_real_real_int_real_real_811;
typedef std::tuple<std::vector<double>, var, int, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_812;
typedef std::tuple<std::vector<double>, var, int, double, var, empty>
    type_vv_real_real_int_real_real_813;
typedef std::tuple<std::vector<double>, var, int, double, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_814;
typedef std::tuple<
    std::vector<double>, var, int, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_815;
typedef std::tuple<std::vector<double>, var, int, std::vector<double>, double,
                   empty>
    type_vv_real_real_int_real_real_816;
typedef std::tuple<std::vector<double>, var, int, std::vector<double>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_817;
typedef std::tuple<std::vector<double>, var, int, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_818;
typedef std::tuple<std::vector<double>, var, int, std::vector<double>, var,
                   empty>
    type_vv_real_real_int_real_real_819;
typedef std::tuple<std::vector<double>, var, int, std::vector<double>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_820;
typedef std::tuple<
    std::vector<double>, var, int, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_821;
typedef std::tuple<std::vector<double>, var, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_822;
typedef std::tuple<std::vector<double>, var, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_823;
typedef std::tuple<std::vector<double>, var, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_824;
typedef std::tuple<std::vector<double>, var, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_825;
typedef std::tuple<std::vector<double>, var, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_826;
typedef std::tuple<
    std::vector<double>, var, int, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_827;
typedef std::tuple<std::vector<double>, var, int, var, double, empty>
    type_vv_real_real_int_real_real_828;
typedef std::tuple<std::vector<double>, var, int, var, std::vector<double>,
                   empty>
    type_vv_real_real_int_real_real_829;
typedef std::tuple<std::vector<double>, var, int, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_830;
typedef std::tuple<std::vector<double>, var, int, var, var, empty>
    type_vv_real_real_int_real_real_831;
typedef std::tuple<std::vector<double>, var, int, var, std::vector<var>, empty>
    type_vv_real_real_int_real_real_832;
typedef std::tuple<
    std::vector<double>, var, int, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_833;
typedef std::tuple<std::vector<double>, var, int, std::vector<var>, double,
                   empty>
    type_vv_real_real_int_real_real_834;
typedef std::tuple<std::vector<double>, var, int, std::vector<var>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_835;
typedef std::tuple<std::vector<double>, var, int, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_836;
typedef std::tuple<std::vector<double>, var, int, std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_837;
typedef std::tuple<std::vector<double>, var, int, std::vector<var>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_838;
typedef std::tuple<
    std::vector<double>, var, int, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_839;
typedef std::tuple<
    std::vector<double>, var, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_840;
typedef std::tuple<
    std::vector<double>, var, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_841;
typedef std::tuple<
    std::vector<double>, var, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_842;
typedef std::tuple<
    std::vector<double>, var, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_843;
typedef std::tuple<
    std::vector<double>, var, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_844;
typedef std::tuple<
    std::vector<double>, var, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_845;
typedef std::tuple<std::vector<double>, var, std::vector<int>, double, double,
                   empty>
    type_vv_real_real_int_real_real_846;
typedef std::tuple<std::vector<double>, var, std::vector<int>, double,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_847;
typedef std::tuple<std::vector<double>, var, std::vector<int>, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_848;
typedef std::tuple<std::vector<double>, var, std::vector<int>, double, var,
                   empty>
    type_vv_real_real_int_real_real_849;
typedef std::tuple<std::vector<double>, var, std::vector<int>, double,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_850;
typedef std::tuple<
    std::vector<double>, var, std::vector<int>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_851;
typedef std::tuple<std::vector<double>, var, std::vector<int>,
                   std::vector<double>, double, empty>
    type_vv_real_real_int_real_real_852;
typedef std::tuple<std::vector<double>, var, std::vector<int>,
                   std::vector<double>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_853;
typedef std::tuple<std::vector<double>, var, std::vector<int>,
                   std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_854;
typedef std::tuple<std::vector<double>, var, std::vector<int>,
                   std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_855;
typedef std::tuple<std::vector<double>, var, std::vector<int>,
                   std::vector<double>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_856;
typedef std::tuple<
    std::vector<double>, var, std::vector<int>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_857;
typedef std::tuple<std::vector<double>, var, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_858;
typedef std::tuple<std::vector<double>, var, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_859;
typedef std::tuple<std::vector<double>, var, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_860;
typedef std::tuple<std::vector<double>, var, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_861;
typedef std::tuple<std::vector<double>, var, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_862;
typedef std::tuple<
    std::vector<double>, var, std::vector<int>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_863;
typedef std::tuple<std::vector<double>, var, std::vector<int>, var, double,
                   empty>
    type_vv_real_real_int_real_real_864;
typedef std::tuple<std::vector<double>, var, std::vector<int>, var,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_865;
typedef std::tuple<std::vector<double>, var, std::vector<int>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_866;
typedef std::tuple<std::vector<double>, var, std::vector<int>, var, var, empty>
    type_vv_real_real_int_real_real_867;
typedef std::tuple<std::vector<double>, var, std::vector<int>, var,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_868;
typedef std::tuple<
    std::vector<double>, var, std::vector<int>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_869;
typedef std::tuple<std::vector<double>, var, std::vector<int>, std::vector<var>,
                   double, empty>
    type_vv_real_real_int_real_real_870;
typedef std::tuple<std::vector<double>, var, std::vector<int>, std::vector<var>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_871;
typedef std::tuple<std::vector<double>, var, std::vector<int>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_872;
typedef std::tuple<std::vector<double>, var, std::vector<int>, std::vector<var>,
                   var, empty>
    type_vv_real_real_int_real_real_873;
typedef std::tuple<std::vector<double>, var, std::vector<int>, std::vector<var>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_874;
typedef std::tuple<
    std::vector<double>, var, std::vector<int>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_875;
typedef std::tuple<
    std::vector<double>, var, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_876;
typedef std::tuple<
    std::vector<double>, var, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_877;
typedef std::tuple<
    std::vector<double>, var, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_878;
typedef std::tuple<
    std::vector<double>, var, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_879;
typedef std::tuple<
    std::vector<double>, var, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_880;
typedef std::tuple<
    std::vector<double>, var, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_881;
typedef std::tuple<std::vector<double>, var,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double, double, empty>
    type_vv_real_real_int_real_real_882;
typedef std::tuple<std::vector<double>, var,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_883;
typedef std::tuple<std::vector<double>, var,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_884;
typedef std::tuple<std::vector<double>, var,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double, var, empty>
    type_vv_real_real_int_real_real_885;
typedef std::tuple<std::vector<double>, var,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_886;
typedef std::tuple<
    std::vector<double>, var, Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_887;
typedef std::tuple<std::vector<double>, var,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   double, empty>
    type_vv_real_real_int_real_real_888;
typedef std::tuple<std::vector<double>, var,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_889;
typedef std::tuple<std::vector<double>, var,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_890;
typedef std::tuple<std::vector<double>, var,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   var, empty>
    type_vv_real_real_int_real_real_891;
typedef std::tuple<std::vector<double>, var,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_892;
typedef std::tuple<
    std::vector<double>, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_893;
typedef std::tuple<std::vector<double>, var,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_894;
typedef std::tuple<
    std::vector<double>, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_895;
typedef std::tuple<std::vector<double>, var,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_896;
typedef std::tuple<std::vector<double>, var,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_897;
typedef std::tuple<
    std::vector<double>, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_898;
typedef std::tuple<
    std::vector<double>, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_899;
typedef std::tuple<std::vector<double>, var,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var, double, empty>
    type_vv_real_real_int_real_real_900;
typedef std::tuple<std::vector<double>, var,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_901;
typedef std::tuple<std::vector<double>, var,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_902;
typedef std::tuple<std::vector<double>, var,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var, var, empty>
    type_vv_real_real_int_real_real_903;
typedef std::tuple<std::vector<double>, var,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_904;
typedef std::tuple<
    std::vector<double>, var, Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_905;
typedef std::tuple<std::vector<double>, var,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   double, empty>
    type_vv_real_real_int_real_real_906;
typedef std::tuple<std::vector<double>, var,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_907;
typedef std::tuple<std::vector<double>, var,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_908;
typedef std::tuple<std::vector<double>, var,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>, var,
                   empty>
    type_vv_real_real_int_real_real_909;
typedef std::tuple<std::vector<double>, var,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_910;
typedef std::tuple<
    std::vector<double>, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_911;
typedef std::tuple<
    std::vector<double>, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_912;
typedef std::tuple<
    std::vector<double>, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_913;
typedef std::tuple<
    std::vector<double>, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_914;
typedef std::tuple<
    std::vector<double>, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_915;
typedef std::tuple<
    std::vector<double>, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_916;
typedef std::tuple<
    std::vector<double>, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_917;
typedef std::tuple<std::vector<double>, std::vector<var>, int, double, double,
                   empty>
    type_vv_real_real_int_real_real_918;
typedef std::tuple<std::vector<double>, std::vector<var>, int, double,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_919;
typedef std::tuple<std::vector<double>, std::vector<var>, int, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_920;
typedef std::tuple<std::vector<double>, std::vector<var>, int, double, var,
                   empty>
    type_vv_real_real_int_real_real_921;
typedef std::tuple<std::vector<double>, std::vector<var>, int, double,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_922;
typedef std::tuple<
    std::vector<double>, std::vector<var>, int, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_923;
typedef std::tuple<std::vector<double>, std::vector<var>, int,
                   std::vector<double>, double, empty>
    type_vv_real_real_int_real_real_924;
typedef std::tuple<std::vector<double>, std::vector<var>, int,
                   std::vector<double>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_925;
typedef std::tuple<std::vector<double>, std::vector<var>, int,
                   std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_926;
typedef std::tuple<std::vector<double>, std::vector<var>, int,
                   std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_927;
typedef std::tuple<std::vector<double>, std::vector<var>, int,
                   std::vector<double>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_928;
typedef std::tuple<
    std::vector<double>, std::vector<var>, int, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_929;
typedef std::tuple<std::vector<double>, std::vector<var>, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_930;
typedef std::tuple<std::vector<double>, std::vector<var>, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_931;
typedef std::tuple<std::vector<double>, std::vector<var>, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_932;
typedef std::tuple<std::vector<double>, std::vector<var>, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_933;
typedef std::tuple<std::vector<double>, std::vector<var>, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_934;
typedef std::tuple<
    std::vector<double>, std::vector<var>, int,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_935;
typedef std::tuple<std::vector<double>, std::vector<var>, int, var, double,
                   empty>
    type_vv_real_real_int_real_real_936;
typedef std::tuple<std::vector<double>, std::vector<var>, int, var,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_937;
typedef std::tuple<std::vector<double>, std::vector<var>, int, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_938;
typedef std::tuple<std::vector<double>, std::vector<var>, int, var, var, empty>
    type_vv_real_real_int_real_real_939;
typedef std::tuple<std::vector<double>, std::vector<var>, int, var,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_940;
typedef std::tuple<
    std::vector<double>, std::vector<var>, int, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_941;
typedef std::tuple<std::vector<double>, std::vector<var>, int, std::vector<var>,
                   double, empty>
    type_vv_real_real_int_real_real_942;
typedef std::tuple<std::vector<double>, std::vector<var>, int, std::vector<var>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_943;
typedef std::tuple<std::vector<double>, std::vector<var>, int, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_944;
typedef std::tuple<std::vector<double>, std::vector<var>, int, std::vector<var>,
                   var, empty>
    type_vv_real_real_int_real_real_945;
typedef std::tuple<std::vector<double>, std::vector<var>, int, std::vector<var>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_946;
typedef std::tuple<
    std::vector<double>, std::vector<var>, int, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_947;
typedef std::tuple<
    std::vector<double>, std::vector<var>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_948;
typedef std::tuple<
    std::vector<double>, std::vector<var>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_949;
typedef std::tuple<
    std::vector<double>, std::vector<var>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_950;
typedef std::tuple<
    std::vector<double>, std::vector<var>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_951;
typedef std::tuple<
    std::vector<double>, std::vector<var>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_952;
typedef std::tuple<
    std::vector<double>, std::vector<var>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_953;
typedef std::tuple<std::vector<double>, std::vector<var>, std::vector<int>,
                   double, double, empty>
    type_vv_real_real_int_real_real_954;
typedef std::tuple<std::vector<double>, std::vector<var>, std::vector<int>,
                   double, std::vector<double>, empty>
    type_vv_real_real_int_real_real_955;
typedef std::tuple<std::vector<double>, std::vector<var>, std::vector<int>,
                   double, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_956;
typedef std::tuple<std::vector<double>, std::vector<var>, std::vector<int>,
                   double, var, empty>
    type_vv_real_real_int_real_real_957;
typedef std::tuple<std::vector<double>, std::vector<var>, std::vector<int>,
                   double, std::vector<var>, empty>
    type_vv_real_real_int_real_real_958;
typedef std::tuple<
    std::vector<double>, std::vector<var>, std::vector<int>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_959;
typedef std::tuple<std::vector<double>, std::vector<var>, std::vector<int>,
                   std::vector<double>, double, empty>
    type_vv_real_real_int_real_real_960;
typedef std::tuple<std::vector<double>, std::vector<var>, std::vector<int>,
                   std::vector<double>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_961;
typedef std::tuple<std::vector<double>, std::vector<var>, std::vector<int>,
                   std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_962;
typedef std::tuple<std::vector<double>, std::vector<var>, std::vector<int>,
                   std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_963;
typedef std::tuple<std::vector<double>, std::vector<var>, std::vector<int>,
                   std::vector<double>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_964;
typedef std::tuple<
    std::vector<double>, std::vector<var>, std::vector<int>,
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_965;
typedef std::tuple<std::vector<double>, std::vector<var>, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_966;
typedef std::tuple<std::vector<double>, std::vector<var>, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_967;
typedef std::tuple<std::vector<double>, std::vector<var>, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_968;
typedef std::tuple<std::vector<double>, std::vector<var>, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_969;
typedef std::tuple<std::vector<double>, std::vector<var>, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_970;
typedef std::tuple<
    std::vector<double>, std::vector<var>, std::vector<int>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_971;
typedef std::tuple<std::vector<double>, std::vector<var>, std::vector<int>, var,
                   double, empty>
    type_vv_real_real_int_real_real_972;
typedef std::tuple<std::vector<double>, std::vector<var>, std::vector<int>, var,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_973;
typedef std::tuple<std::vector<double>, std::vector<var>, std::vector<int>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_974;
typedef std::tuple<std::vector<double>, std::vector<var>, std::vector<int>, var,
                   var, empty>
    type_vv_real_real_int_real_real_975;
typedef std::tuple<std::vector<double>, std::vector<var>, std::vector<int>, var,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_976;
typedef std::tuple<
    std::vector<double>, std::vector<var>, std::vector<int>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_977;
typedef std::tuple<std::vector<double>, std::vector<var>, std::vector<int>,
                   std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_978;
typedef std::tuple<std::vector<double>, std::vector<var>, std::vector<int>,
                   std::vector<var>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_979;
typedef std::tuple<std::vector<double>, std::vector<var>, std::vector<int>,
                   std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   empty>
    type_vv_real_real_int_real_real_980;
typedef std::tuple<std::vector<double>, std::vector<var>, std::vector<int>,
                   std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_981;
typedef std::tuple<std::vector<double>, std::vector<var>, std::vector<int>,
                   std::vector<var>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_982;
typedef std::tuple<
    std::vector<double>, std::vector<var>, std::vector<int>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_983;
typedef std::tuple<
    std::vector<double>, std::vector<var>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_984;
typedef std::tuple<
    std::vector<double>, std::vector<var>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_985;
typedef std::tuple<
    std::vector<double>, std::vector<var>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_986;
typedef std::tuple<
    std::vector<double>, std::vector<var>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_987;
typedef std::tuple<
    std::vector<double>, std::vector<var>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_988;
typedef std::tuple<
    std::vector<double>, std::vector<var>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_989;
typedef std::tuple<std::vector<double>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double, double, empty>
    type_vv_real_real_int_real_real_990;
typedef std::tuple<std::vector<double>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_991;
typedef std::tuple<std::vector<double>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_992;
typedef std::tuple<std::vector<double>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double, var, empty>
    type_vv_real_real_int_real_real_993;
typedef std::tuple<std::vector<double>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_994;
typedef std::tuple<
    std::vector<double>, std::vector<var>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_995;
typedef std::tuple<std::vector<double>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   double, empty>
    type_vv_real_real_int_real_real_996;
typedef std::tuple<std::vector<double>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_997;
typedef std::tuple<std::vector<double>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_998;
typedef std::tuple<std::vector<double>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   var, empty>
    type_vv_real_real_int_real_real_999;
typedef std::tuple<std::vector<double>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_1000;
typedef std::tuple<
    std::vector<double>, std::vector<var>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1001;
typedef std::tuple<std::vector<double>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_1002;
typedef std::tuple<std::vector<double>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_1003;
typedef std::tuple<std::vector<double>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1004;
typedef std::tuple<std::vector<double>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_1005;
typedef std::tuple<std::vector<double>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_1006;
typedef std::tuple<
    std::vector<double>, std::vector<var>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1007;
typedef std::tuple<std::vector<double>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var, double, empty>
    type_vv_real_real_int_real_real_1008;
typedef std::tuple<std::vector<double>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_1009;
typedef std::tuple<std::vector<double>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1010;
typedef std::tuple<std::vector<double>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var, var, empty>
    type_vv_real_real_int_real_real_1011;
typedef std::tuple<std::vector<double>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_1012;
typedef std::tuple<
    std::vector<double>, std::vector<var>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1013;
typedef std::tuple<std::vector<double>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   double, empty>
    type_vv_real_real_int_real_real_1014;
typedef std::tuple<std::vector<double>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_1015;
typedef std::tuple<std::vector<double>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1016;
typedef std::tuple<std::vector<double>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>, var,
                   empty>
    type_vv_real_real_int_real_real_1017;
typedef std::tuple<std::vector<double>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_1018;
typedef std::tuple<
    std::vector<double>, std::vector<var>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1019;
typedef std::tuple<
    std::vector<double>, std::vector<var>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_1020;
typedef std::tuple<
    std::vector<double>, std::vector<var>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_1021;
typedef std::tuple<
    std::vector<double>, std::vector<var>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1022;
typedef std::tuple<
    std::vector<double>, std::vector<var>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_1023;
typedef std::tuple<
    std::vector<double>, std::vector<var>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_1024;
typedef std::tuple<
    std::vector<double>, std::vector<var>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1025;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    double, double, empty>
    type_vv_real_real_int_real_real_1026;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    double, std::vector<double>, empty>
    type_vv_real_real_int_real_real_1027;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    double, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1028;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    double, var, empty>
    type_vv_real_real_int_real_real_1029;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    double, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1030;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    empty>
    type_vv_real_real_int_real_real_1031;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<double>, double, empty>
    type_vv_real_real_int_real_real_1032;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<double>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_1033;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1034;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_1035;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<double>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1036;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1037;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_1038;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_1039;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1040;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_1041;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1042;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1043;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int, var,
    double, empty>
    type_vv_real_real_int_real_real_1044;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int, var,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_1045;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int, var,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1046;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int, var,
    var, empty>
    type_vv_real_real_int_real_real_1047;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int, var,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_1048;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1049;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_1050;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<var>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_1051;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1052;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_1053;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<var>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1054;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1055;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_1056;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_1057;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1058;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_1059;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_1060;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1061;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, double, double, empty>
    type_vv_real_real_int_real_real_1062;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, double, std::vector<double>, empty>
    type_vv_real_real_int_real_real_1063;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, double, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1064;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, double, var, empty>
    type_vv_real_real_int_real_real_1065;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, double, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1066;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1067;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<double>, double, empty>
    type_vv_real_real_int_real_real_1068;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<double>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_1069;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<double>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1070;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_1071;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<double>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1072;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1073;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_1074;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_1075;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1076;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_1077;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_1078;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1079;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, var, double, empty>
    type_vv_real_real_int_real_real_1080;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, var, std::vector<double>, empty>
    type_vv_real_real_int_real_real_1081;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, var, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1082;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, var, var, empty>
    type_vv_real_real_int_real_real_1083;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, var, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1084;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1085;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_1086;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<var>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_1087;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<var>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1088;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_1089;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<var>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1090;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1091;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_1092;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_1093;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1094;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_1095;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_1096;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1097;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double, double, empty>
    type_vv_real_real_int_real_real_1098;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double, std::vector<double>, empty>
    type_vv_real_real_int_real_real_1099;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1100;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double, var, empty>
    type_vv_real_real_int_real_real_1101;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1102;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1103;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>, double, empty>
    type_vv_real_real_int_real_real_1104;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_1105;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1106;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_1107;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_1108;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1109;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_1110;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_1111;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1112;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_1113;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1114;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1115;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var, double, empty>
    type_vv_real_real_int_real_real_1116;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var, std::vector<double>, empty>
    type_vv_real_real_int_real_real_1117;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1118;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var, var, empty>
    type_vv_real_real_int_real_real_1119;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1120;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1121;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_1122;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_1123;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1124;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_1125;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>, std::vector<var>,
    empty>
    type_vv_real_real_int_real_real_1126;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1127;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_1128;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_1129;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1130;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_1131;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_1132;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1133;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double, int,
                   double, var, empty>
    type_vv_real_real_int_real_real_1134;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double, int,
                   double, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1135;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, int, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1136;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double, int,
                   std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_1137;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double, int,
                   std::vector<double>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1138;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, int, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1139;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_1140;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_1141;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, int,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1142;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double, int, var,
                   double, empty>
    type_vv_real_real_int_real_real_1143;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double, int, var,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_1144;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double, int, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1145;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double, int, var,
                   var, empty>
    type_vv_real_real_int_real_real_1146;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double, int, var,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_1147;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, int, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1148;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double, int,
                   std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_1149;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double, int,
                   std::vector<var>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_1150;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double, int,
                   std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   empty>
    type_vv_real_real_int_real_real_1151;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double, int,
                   std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_1152;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double, int,
                   std::vector<var>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1153;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, int, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1154;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_1155;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_1156;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1157;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_1158;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_1159;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1160;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
                   std::vector<int>, double, var, empty>
    type_vv_real_real_int_real_real_1161;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
                   std::vector<int>, double, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1162;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, std::vector<int>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1163;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
                   std::vector<int>, std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_1164;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
                   std::vector<int>, std::vector<double>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_1165;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, std::vector<int>,
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1166;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
                   std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   var, empty>
    type_vv_real_real_int_real_real_1167;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
                   std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_1168;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, std::vector<int>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1169;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
                   std::vector<int>, var, double, empty>
    type_vv_real_real_int_real_real_1170;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
                   std::vector<int>, var, std::vector<double>, empty>
    type_vv_real_real_int_real_real_1171;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
                   std::vector<int>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1172;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
                   std::vector<int>, var, var, empty>
    type_vv_real_real_int_real_real_1173;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
                   std::vector<int>, var, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1174;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, std::vector<int>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1175;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
                   std::vector<int>, std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_1176;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
                   std::vector<int>, std::vector<var>, std::vector<double>,
                   empty>
    type_vv_real_real_int_real_real_1177;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
                   std::vector<int>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1178;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
                   std::vector<int>, std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_1179;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
                   std::vector<int>, std::vector<var>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1180;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, std::vector<int>,
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1181;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_1182;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_1183;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1184;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_1185;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_1186;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1187;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double, var, empty>
    type_vv_real_real_int_real_real_1188;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_1189;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1190;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   var, empty>
    type_vv_real_real_int_real_real_1191;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_1192;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1193;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_1194;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_1195;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1196;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var, double, empty>
    type_vv_real_real_int_real_real_1197;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_1198;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1199;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var, var, empty>
    type_vv_real_real_int_real_real_1200;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_1201;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1202;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   double, empty>
    type_vv_real_real_int_real_real_1203;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_1204;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1205;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>, var,
                   empty>
    type_vv_real_real_int_real_real_1206;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_1207;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1208;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_1209;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_1210;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1211;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_1212;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_1213;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1214;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, int, double, var, empty>
    type_vv_real_real_int_real_real_1215;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, int, double, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1216;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, int, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1217;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, int, std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_1218;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, int, std::vector<double>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_1219;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, int,
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1220;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_1221;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, int,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1222;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, int,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1223;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, int, var, double, empty>
    type_vv_real_real_int_real_real_1224;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, int, var, std::vector<double>, empty>
    type_vv_real_real_int_real_real_1225;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, int, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1226;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, int, var, var, empty>
    type_vv_real_real_int_real_real_1227;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, int, var, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1228;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, int, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1229;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, int, std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_1230;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, int, std::vector<var>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_1231;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, int, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1232;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, int, std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_1233;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, int, std::vector<var>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_1234;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, int,
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1235;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_1236;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_1237;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1238;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_1239;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_1240;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1241;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, std::vector<int>, double, var, empty>
    type_vv_real_real_int_real_real_1242;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, std::vector<int>, double,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_1243;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>,
    std::vector<int>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1244;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, std::vector<int>, std::vector<double>,
                   var, empty>
    type_vv_real_real_int_real_real_1245;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, std::vector<int>, std::vector<double>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_1246;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>,
    std::vector<int>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1247;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_1248;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_1249;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>,
    std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1250;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, std::vector<int>, var, double, empty>
    type_vv_real_real_int_real_real_1251;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, std::vector<int>, var,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_1252;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, std::vector<int>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1253;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, std::vector<int>, var, var, empty>
    type_vv_real_real_int_real_real_1254;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, std::vector<int>, var, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_1255;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>,
    std::vector<int>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1256;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, std::vector<int>, std::vector<var>,
                   double, empty>
    type_vv_real_real_int_real_real_1257;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, std::vector<int>, std::vector<var>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_1258;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, std::vector<int>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1259;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, std::vector<int>, std::vector<var>, var,
                   empty>
    type_vv_real_real_int_real_real_1260;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, std::vector<int>, std::vector<var>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_1261;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>,
    std::vector<int>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1262;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_1263;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_1264;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1265;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_1266;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_1267;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1268;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   double, var, empty>
    type_vv_real_real_int_real_real_1269;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   double, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1270;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1271;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_1272;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<double>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1273;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1274;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_1275;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_1276;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1277;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   var, double, empty>
    type_vv_real_real_int_real_real_1278;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   var, std::vector<double>, empty>
    type_vv_real_real_int_real_real_1279;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   var, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1280;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   var, var, empty>
    type_vv_real_real_int_real_real_1281;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   var, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1282;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1283;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_1284;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<var>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_1285;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   empty>
    type_vv_real_real_int_real_real_1286;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_1287;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<var>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1288;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1289;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_1290;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_1291;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1292;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_1293;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_1294;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1295;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, int, double, var,
                   empty>
    type_vv_real_real_int_real_real_1296;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, int, double,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_1297;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, int, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1298;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
                   std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_1299;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
                   std::vector<double>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1300;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, int, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1301;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_1302;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_1303;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1304;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, int, var, double,
                   empty>
    type_vv_real_real_int_real_real_1305;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, int, var,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_1306;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, int, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1307;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, int, var, var,
                   empty>
    type_vv_real_real_int_real_real_1308;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, int, var,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_1309;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, int, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1310;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
                   std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_1311;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
                   std::vector<var>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_1312;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
                   std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   empty>
    type_vv_real_real_int_real_real_1313;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
                   std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_1314;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
                   std::vector<var>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1315;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, int, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1316;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_1317;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_1318;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1319;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_1320;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_1321;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1322;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
                   double, var, empty>
    type_vv_real_real_int_real_real_1323;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
                   double, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1324;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1325;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
                   std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_1326;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
                   std::vector<double>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1327;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1328;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_1329;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_1330;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1331;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
                   var, double, empty>
    type_vv_real_real_int_real_real_1332;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
                   var, std::vector<double>, empty>
    type_vv_real_real_int_real_real_1333;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
                   var, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1334;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
                   var, var, empty>
    type_vv_real_real_int_real_real_1335;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
                   var, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1336;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1337;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
                   std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_1338;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
                   std::vector<var>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_1339;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
                   std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   empty>
    type_vv_real_real_int_real_real_1340;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
                   std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_1341;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
                   std::vector<var>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1342;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1343;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_1344;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_1345;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1346;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_1347;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_1348;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1349;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double, var, empty>
    type_vv_real_real_int_real_real_1350;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_1351;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1352;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   var, empty>
    type_vv_real_real_int_real_real_1353;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_1354;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1355;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_1356;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_1357;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1358;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var, double, empty>
    type_vv_real_real_int_real_real_1359;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_1360;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1361;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var, var, empty>
    type_vv_real_real_int_real_real_1362;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_1363;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1364;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   double, empty>
    type_vv_real_real_int_real_real_1365;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_1366;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1367;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>, var,
                   empty>
    type_vv_real_real_int_real_real_1368;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_1369;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1370;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_1371;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_1372;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1373;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_1374;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_1375;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1376;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var, int, double,
                   double, empty>
    type_vv_real_real_int_real_real_1377;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var, int, double,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_1378;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var, int, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1379;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var, int, double,
                   var, empty>
    type_vv_real_real_int_real_real_1380;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var, int, double,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_1381;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, int, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1382;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var, int,
                   std::vector<double>, double, empty>
    type_vv_real_real_int_real_real_1383;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var, int,
                   std::vector<double>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_1384;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var, int,
                   std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1385;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var, int,
                   std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_1386;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var, int,
                   std::vector<double>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1387;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, int, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1388;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_1389;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_1390;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1391;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_1392;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_1393;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, int,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1394;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var, int, var,
                   double, empty>
    type_vv_real_real_int_real_real_1395;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var, int, var,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_1396;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var, int, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1397;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var, int, var, var,
                   empty>
    type_vv_real_real_int_real_real_1398;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var, int, var,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_1399;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, int, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1400;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var, int,
                   std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_1401;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var, int,
                   std::vector<var>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_1402;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var, int,
                   std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   empty>
    type_vv_real_real_int_real_real_1403;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var, int,
                   std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_1404;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var, int,
                   std::vector<var>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1405;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, int, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1406;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_1407;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_1408;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1409;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_1410;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_1411;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1412;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   std::vector<int>, double, double, empty>
    type_vv_real_real_int_real_real_1413;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   std::vector<int>, double, std::vector<double>, empty>
    type_vv_real_real_int_real_real_1414;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   std::vector<int>, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1415;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   std::vector<int>, double, var, empty>
    type_vv_real_real_int_real_real_1416;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   std::vector<int>, double, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1417;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, std::vector<int>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1418;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   std::vector<int>, std::vector<double>, double, empty>
    type_vv_real_real_int_real_real_1419;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   std::vector<int>, std::vector<double>, std::vector<double>,
                   empty>
    type_vv_real_real_int_real_real_1420;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   std::vector<int>, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1421;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   std::vector<int>, std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_1422;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   std::vector<int>, std::vector<double>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_1423;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, std::vector<int>,
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1424;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   double, empty>
    type_vv_real_real_int_real_real_1425;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_1426;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1427;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   var, empty>
    type_vv_real_real_int_real_real_1428;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_1429;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, std::vector<int>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1430;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   std::vector<int>, var, double, empty>
    type_vv_real_real_int_real_real_1431;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   std::vector<int>, var, std::vector<double>, empty>
    type_vv_real_real_int_real_real_1432;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   std::vector<int>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1433;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   std::vector<int>, var, var, empty>
    type_vv_real_real_int_real_real_1434;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   std::vector<int>, var, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1435;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, std::vector<int>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1436;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   std::vector<int>, std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_1437;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   std::vector<int>, std::vector<var>, std::vector<double>,
                   empty>
    type_vv_real_real_int_real_real_1438;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   std::vector<int>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1439;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   std::vector<int>, std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_1440;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   std::vector<int>, std::vector<var>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1441;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, std::vector<int>,
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1442;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_1443;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_1444;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1445;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_1446;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_1447;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1448;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double, double, empty>
    type_vv_real_real_int_real_real_1449;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_1450;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1451;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double, var, empty>
    type_vv_real_real_int_real_real_1452;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_1453;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1454;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   double, empty>
    type_vv_real_real_int_real_real_1455;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_1456;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1457;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   var, empty>
    type_vv_real_real_int_real_real_1458;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_1459;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1460;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_1461;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_1462;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1463;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_1464;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_1465;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1466;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var, double, empty>
    type_vv_real_real_int_real_real_1467;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_1468;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1469;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var, var, empty>
    type_vv_real_real_int_real_real_1470;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_1471;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1472;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   double, empty>
    type_vv_real_real_int_real_real_1473;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_1474;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1475;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>, var,
                   empty>
    type_vv_real_real_int_real_real_1476;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_1477;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1478;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_1479;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_1480;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1481;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_1482;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_1483;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1484;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   int, double, double, empty>
    type_vv_real_real_int_real_real_1485;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   int, double, std::vector<double>, empty>
    type_vv_real_real_int_real_real_1486;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   int, double, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1487;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   int, double, var, empty>
    type_vv_real_real_int_real_real_1488;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   int, double, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1489;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, int, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1490;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   int, std::vector<double>, double, empty>
    type_vv_real_real_int_real_real_1491;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   int, std::vector<double>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_1492;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   int, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1493;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   int, std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_1494;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   int, std::vector<double>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1495;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, int,
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1496;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   int, Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_1497;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   int, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_1498;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   int, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1499;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   int, Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_1500;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   int, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_1501;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, int,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1502;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   int, var, double, empty>
    type_vv_real_real_int_real_real_1503;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   int, var, std::vector<double>, empty>
    type_vv_real_real_int_real_real_1504;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   int, var, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1505;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   int, var, var, empty>
    type_vv_real_real_int_real_real_1506;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   int, var, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1507;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, int, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1508;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   int, std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_1509;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   int, std::vector<var>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_1510;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   int, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1511;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   int, std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_1512;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   int, std::vector<var>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1513;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, int,
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1514;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_1515;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_1516;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1517;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_1518;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_1519;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1520;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<int>, double, double, empty>
    type_vv_real_real_int_real_real_1521;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<int>, double, std::vector<double>, empty>
    type_vv_real_real_int_real_real_1522;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<int>, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1523;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<int>, double, var, empty>
    type_vv_real_real_int_real_real_1524;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<int>, double, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1525;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
    std::vector<int>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1526;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<int>, std::vector<double>, double, empty>
    type_vv_real_real_int_real_real_1527;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<int>, std::vector<double>, std::vector<double>,
                   empty>
    type_vv_real_real_int_real_real_1528;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<int>, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1529;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<int>, std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_1530;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<int>, std::vector<double>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_1531;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
    std::vector<int>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1532;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   double, empty>
    type_vv_real_real_int_real_real_1533;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_1534;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1535;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   var, empty>
    type_vv_real_real_int_real_real_1536;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_1537;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
    std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1538;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<int>, var, double, empty>
    type_vv_real_real_int_real_real_1539;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<int>, var, std::vector<double>, empty>
    type_vv_real_real_int_real_real_1540;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<int>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1541;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<int>, var, var, empty>
    type_vv_real_real_int_real_real_1542;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<int>, var, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1543;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
    std::vector<int>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1544;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<int>, std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_1545;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<int>, std::vector<var>, std::vector<double>,
                   empty>
    type_vv_real_real_int_real_real_1546;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<int>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1547;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<int>, std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_1548;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<int>, std::vector<var>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1549;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
    std::vector<int>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1550;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_1551;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_1552;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1553;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_1554;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_1555;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1556;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double, double, empty>
    type_vv_real_real_int_real_real_1557;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_1558;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1559;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double, var, empty>
    type_vv_real_real_int_real_real_1560;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_1561;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1562;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   double, empty>
    type_vv_real_real_int_real_real_1563;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_1564;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1565;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   var, empty>
    type_vv_real_real_int_real_real_1566;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_1567;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1568;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_1569;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_1570;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1571;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_1572;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_1573;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1574;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var, double, empty>
    type_vv_real_real_int_real_real_1575;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_1576;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1577;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var, var, empty>
    type_vv_real_real_int_real_real_1578;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_1579;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1580;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   double, empty>
    type_vv_real_real_int_real_real_1581;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_1582;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1583;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>, var,
                   empty>
    type_vv_real_real_int_real_real_1584;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_1585;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1586;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_1587;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_1588;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1589;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_1590;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_1591;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1592;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    double, double, empty>
    type_vv_real_real_int_real_real_1593;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    double, std::vector<double>, empty>
    type_vv_real_real_int_real_real_1594;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    double, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1595;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    double, var, empty>
    type_vv_real_real_int_real_real_1596;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    double, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1597;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    empty>
    type_vv_real_real_int_real_real_1598;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<double>, double, empty>
    type_vv_real_real_int_real_real_1599;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<double>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_1600;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1601;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_1602;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<double>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1603;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1604;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_1605;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_1606;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1607;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_1608;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1609;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1610;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int, var,
    double, empty>
    type_vv_real_real_int_real_real_1611;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int, var,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_1612;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int, var,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1613;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int, var,
    var, empty>
    type_vv_real_real_int_real_real_1614;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int, var,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_1615;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1616;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_1617;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<var>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_1618;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1619;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_1620;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<var>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1621;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1622;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_1623;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_1624;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1625;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_1626;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_1627;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1628;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, double, double, empty>
    type_vv_real_real_int_real_real_1629;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, double, std::vector<double>, empty>
    type_vv_real_real_int_real_real_1630;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, double, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1631;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, double, var, empty>
    type_vv_real_real_int_real_real_1632;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, double, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1633;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1634;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<double>, double, empty>
    type_vv_real_real_int_real_real_1635;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<double>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_1636;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<double>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1637;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_1638;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<double>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1639;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1640;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_1641;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_1642;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1643;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_1644;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_1645;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1646;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, var, double, empty>
    type_vv_real_real_int_real_real_1647;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, var, std::vector<double>, empty>
    type_vv_real_real_int_real_real_1648;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, var, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1649;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, var, var, empty>
    type_vv_real_real_int_real_real_1650;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, var, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1651;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1652;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_1653;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<var>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_1654;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<var>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1655;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_1656;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<var>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1657;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1658;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_1659;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_1660;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1661;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_1662;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_1663;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1664;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double, double, empty>
    type_vv_real_real_int_real_real_1665;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double, std::vector<double>, empty>
    type_vv_real_real_int_real_real_1666;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1667;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double, var, empty>
    type_vv_real_real_int_real_real_1668;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1669;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1670;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>, double, empty>
    type_vv_real_real_int_real_real_1671;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_1672;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1673;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_1674;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_1675;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1676;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_1677;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_1678;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1679;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_1680;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1681;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1682;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var, double, empty>
    type_vv_real_real_int_real_real_1683;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var, std::vector<double>, empty>
    type_vv_real_real_int_real_real_1684;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1685;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var, var, empty>
    type_vv_real_real_int_real_real_1686;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1687;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1688;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_1689;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_1690;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1691;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_1692;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>, std::vector<var>,
    empty>
    type_vv_real_real_int_real_real_1693;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1694;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_1695;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_1696;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1697;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_1698;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_1699;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1700;
typedef std::tuple<var, double, int, double, double, empty>
    type_vv_real_real_int_real_real_1701;
typedef std::tuple<var, double, int, double, std::vector<double>, empty>
    type_vv_real_real_int_real_real_1702;
typedef std::tuple<var, double, int, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1703;
typedef std::tuple<var, double, int, double, var, empty>
    type_vv_real_real_int_real_real_1704;
typedef std::tuple<var, double, int, double, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1705;
typedef std::tuple<
    var, double, int, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1706;
typedef std::tuple<var, double, int, std::vector<double>, double, empty>
    type_vv_real_real_int_real_real_1707;
typedef std::tuple<var, double, int, std::vector<double>, std::vector<double>,
                   empty>
    type_vv_real_real_int_real_real_1708;
typedef std::tuple<var, double, int, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1709;
typedef std::tuple<var, double, int, std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_1710;
typedef std::tuple<var, double, int, std::vector<double>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_1711;
typedef std::tuple<
    var, double, int, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1712;
typedef std::tuple<var, double, int, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   double, empty>
    type_vv_real_real_int_real_real_1713;
typedef std::tuple<var, double, int, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_1714;
typedef std::tuple<var, double, int, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1715;
typedef std::tuple<var, double, int, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   var, empty>
    type_vv_real_real_int_real_real_1716;
typedef std::tuple<var, double, int, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_1717;
typedef std::tuple<
    var, double, int, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1718;
typedef std::tuple<var, double, int, var, double, empty>
    type_vv_real_real_int_real_real_1719;
typedef std::tuple<var, double, int, var, std::vector<double>, empty>
    type_vv_real_real_int_real_real_1720;
typedef std::tuple<var, double, int, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1721;
typedef std::tuple<var, double, int, var, var, empty>
    type_vv_real_real_int_real_real_1722;
typedef std::tuple<var, double, int, var, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1723;
typedef std::tuple<
    var, double, int, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1724;
typedef std::tuple<var, double, int, std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_1725;
typedef std::tuple<var, double, int, std::vector<var>, std::vector<double>,
                   empty>
    type_vv_real_real_int_real_real_1726;
typedef std::tuple<var, double, int, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1727;
typedef std::tuple<var, double, int, std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_1728;
typedef std::tuple<var, double, int, std::vector<var>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1729;
typedef std::tuple<
    var, double, int, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1730;
typedef std::tuple<
    var, double, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_1731;
typedef std::tuple<
    var, double, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_1732;
typedef std::tuple<
    var, double, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1733;
typedef std::tuple<
    var, double, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_1734;
typedef std::tuple<
    var, double, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_1735;
typedef std::tuple<
    var, double, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1736;
typedef std::tuple<var, double, std::vector<int>, double, double, empty>
    type_vv_real_real_int_real_real_1737;
typedef std::tuple<var, double, std::vector<int>, double, std::vector<double>,
                   empty>
    type_vv_real_real_int_real_real_1738;
typedef std::tuple<var, double, std::vector<int>, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1739;
typedef std::tuple<var, double, std::vector<int>, double, var, empty>
    type_vv_real_real_int_real_real_1740;
typedef std::tuple<var, double, std::vector<int>, double, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_1741;
typedef std::tuple<
    var, double, std::vector<int>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1742;
typedef std::tuple<var, double, std::vector<int>, std::vector<double>, double,
                   empty>
    type_vv_real_real_int_real_real_1743;
typedef std::tuple<var, double, std::vector<int>, std::vector<double>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_1744;
typedef std::tuple<var, double, std::vector<int>, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1745;
typedef std::tuple<var, double, std::vector<int>, std::vector<double>, var,
                   empty>
    type_vv_real_real_int_real_real_1746;
typedef std::tuple<var, double, std::vector<int>, std::vector<double>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_1747;
typedef std::tuple<
    var, double, std::vector<int>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1748;
typedef std::tuple<var, double, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_1749;
typedef std::tuple<var, double, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_1750;
typedef std::tuple<var, double, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1751;
typedef std::tuple<var, double, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_1752;
typedef std::tuple<var, double, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_1753;
typedef std::tuple<
    var, double, std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1754;
typedef std::tuple<var, double, std::vector<int>, var, double, empty>
    type_vv_real_real_int_real_real_1755;
typedef std::tuple<var, double, std::vector<int>, var, std::vector<double>,
                   empty>
    type_vv_real_real_int_real_real_1756;
typedef std::tuple<var, double, std::vector<int>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1757;
typedef std::tuple<var, double, std::vector<int>, var, var, empty>
    type_vv_real_real_int_real_real_1758;
typedef std::tuple<var, double, std::vector<int>, var, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1759;
typedef std::tuple<
    var, double, std::vector<int>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1760;
typedef std::tuple<var, double, std::vector<int>, std::vector<var>, double,
                   empty>
    type_vv_real_real_int_real_real_1761;
typedef std::tuple<var, double, std::vector<int>, std::vector<var>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_1762;
typedef std::tuple<var, double, std::vector<int>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1763;
typedef std::tuple<var, double, std::vector<int>, std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_1764;
typedef std::tuple<var, double, std::vector<int>, std::vector<var>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_1765;
typedef std::tuple<
    var, double, std::vector<int>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1766;
typedef std::tuple<
    var, double, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_1767;
typedef std::tuple<
    var, double, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_1768;
typedef std::tuple<
    var, double, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1769;
typedef std::tuple<
    var, double, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_1770;
typedef std::tuple<
    var, double, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_1771;
typedef std::tuple<
    var, double, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1772;
typedef std::tuple<var, double, Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
                   double, empty>
    type_vv_real_real_int_real_real_1773;
typedef std::tuple<var, double, Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_1774;
typedef std::tuple<var, double, Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1775;
typedef std::tuple<var, double, Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
                   var, empty>
    type_vv_real_real_int_real_real_1776;
typedef std::tuple<var, double, Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_1777;
typedef std::tuple<
    var, double, Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1778;
typedef std::tuple<var, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<double>, double, empty>
    type_vv_real_real_int_real_real_1779;
typedef std::tuple<var, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<double>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_1780;
typedef std::tuple<var, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1781;
typedef std::tuple<var, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_1782;
typedef std::tuple<var, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<double>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1783;
typedef std::tuple<
    var, double, Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1784;
typedef std::tuple<var, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_1785;
typedef std::tuple<var, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_1786;
typedef std::tuple<var, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1787;
typedef std::tuple<var, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_1788;
typedef std::tuple<var, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_1789;
typedef std::tuple<
    var, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1790;
typedef std::tuple<var, double, Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
                   double, empty>
    type_vv_real_real_int_real_real_1791;
typedef std::tuple<var, double, Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_1792;
typedef std::tuple<var, double, Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1793;
typedef std::tuple<var, double, Eigen::Matrix<int, Eigen::Dynamic, 1>, var, var,
                   empty>
    type_vv_real_real_int_real_real_1794;
typedef std::tuple<var, double, Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_1795;
typedef std::tuple<
    var, double, Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1796;
typedef std::tuple<var, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_1797;
typedef std::tuple<var, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<var>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_1798;
typedef std::tuple<var, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   empty>
    type_vv_real_real_int_real_real_1799;
typedef std::tuple<var, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_1800;
typedef std::tuple<var, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<var>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1801;
typedef std::tuple<
    var, double, Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1802;
typedef std::tuple<
    var, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_1803;
typedef std::tuple<
    var, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_1804;
typedef std::tuple<
    var, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1805;
typedef std::tuple<
    var, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_1806;
typedef std::tuple<
    var, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_1807;
typedef std::tuple<
    var, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1808;
typedef std::tuple<var, std::vector<double>, int, double, double, empty>
    type_vv_real_real_int_real_real_1809;
typedef std::tuple<var, std::vector<double>, int, double, std::vector<double>,
                   empty>
    type_vv_real_real_int_real_real_1810;
typedef std::tuple<var, std::vector<double>, int, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1811;
typedef std::tuple<var, std::vector<double>, int, double, var, empty>
    type_vv_real_real_int_real_real_1812;
typedef std::tuple<var, std::vector<double>, int, double, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_1813;
typedef std::tuple<
    var, std::vector<double>, int, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1814;
typedef std::tuple<var, std::vector<double>, int, std::vector<double>, double,
                   empty>
    type_vv_real_real_int_real_real_1815;
typedef std::tuple<var, std::vector<double>, int, std::vector<double>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_1816;
typedef std::tuple<var, std::vector<double>, int, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1817;
typedef std::tuple<var, std::vector<double>, int, std::vector<double>, var,
                   empty>
    type_vv_real_real_int_real_real_1818;
typedef std::tuple<var, std::vector<double>, int, std::vector<double>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_1819;
typedef std::tuple<
    var, std::vector<double>, int, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1820;
typedef std::tuple<var, std::vector<double>, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_1821;
typedef std::tuple<var, std::vector<double>, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_1822;
typedef std::tuple<var, std::vector<double>, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1823;
typedef std::tuple<var, std::vector<double>, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_1824;
typedef std::tuple<var, std::vector<double>, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_1825;
typedef std::tuple<
    var, std::vector<double>, int, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1826;
typedef std::tuple<var, std::vector<double>, int, var, double, empty>
    type_vv_real_real_int_real_real_1827;
typedef std::tuple<var, std::vector<double>, int, var, std::vector<double>,
                   empty>
    type_vv_real_real_int_real_real_1828;
typedef std::tuple<var, std::vector<double>, int, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1829;
typedef std::tuple<var, std::vector<double>, int, var, var, empty>
    type_vv_real_real_int_real_real_1830;
typedef std::tuple<var, std::vector<double>, int, var, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1831;
typedef std::tuple<
    var, std::vector<double>, int, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1832;
typedef std::tuple<var, std::vector<double>, int, std::vector<var>, double,
                   empty>
    type_vv_real_real_int_real_real_1833;
typedef std::tuple<var, std::vector<double>, int, std::vector<var>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_1834;
typedef std::tuple<var, std::vector<double>, int, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1835;
typedef std::tuple<var, std::vector<double>, int, std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_1836;
typedef std::tuple<var, std::vector<double>, int, std::vector<var>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_1837;
typedef std::tuple<
    var, std::vector<double>, int, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1838;
typedef std::tuple<
    var, std::vector<double>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_1839;
typedef std::tuple<
    var, std::vector<double>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_1840;
typedef std::tuple<
    var, std::vector<double>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1841;
typedef std::tuple<
    var, std::vector<double>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_1842;
typedef std::tuple<
    var, std::vector<double>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_1843;
typedef std::tuple<
    var, std::vector<double>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1844;
typedef std::tuple<var, std::vector<double>, std::vector<int>, double, double,
                   empty>
    type_vv_real_real_int_real_real_1845;
typedef std::tuple<var, std::vector<double>, std::vector<int>, double,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_1846;
typedef std::tuple<var, std::vector<double>, std::vector<int>, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1847;
typedef std::tuple<var, std::vector<double>, std::vector<int>, double, var,
                   empty>
    type_vv_real_real_int_real_real_1848;
typedef std::tuple<var, std::vector<double>, std::vector<int>, double,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_1849;
typedef std::tuple<
    var, std::vector<double>, std::vector<int>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1850;
typedef std::tuple<var, std::vector<double>, std::vector<int>,
                   std::vector<double>, double, empty>
    type_vv_real_real_int_real_real_1851;
typedef std::tuple<var, std::vector<double>, std::vector<int>,
                   std::vector<double>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_1852;
typedef std::tuple<var, std::vector<double>, std::vector<int>,
                   std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1853;
typedef std::tuple<var, std::vector<double>, std::vector<int>,
                   std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_1854;
typedef std::tuple<var, std::vector<double>, std::vector<int>,
                   std::vector<double>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1855;
typedef std::tuple<
    var, std::vector<double>, std::vector<int>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1856;
typedef std::tuple<var, std::vector<double>, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_1857;
typedef std::tuple<var, std::vector<double>, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_1858;
typedef std::tuple<var, std::vector<double>, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1859;
typedef std::tuple<var, std::vector<double>, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_1860;
typedef std::tuple<var, std::vector<double>, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_1861;
typedef std::tuple<
    var, std::vector<double>, std::vector<int>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1862;
typedef std::tuple<var, std::vector<double>, std::vector<int>, var, double,
                   empty>
    type_vv_real_real_int_real_real_1863;
typedef std::tuple<var, std::vector<double>, std::vector<int>, var,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_1864;
typedef std::tuple<var, std::vector<double>, std::vector<int>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1865;
typedef std::tuple<var, std::vector<double>, std::vector<int>, var, var, empty>
    type_vv_real_real_int_real_real_1866;
typedef std::tuple<var, std::vector<double>, std::vector<int>, var,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_1867;
typedef std::tuple<
    var, std::vector<double>, std::vector<int>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1868;
typedef std::tuple<var, std::vector<double>, std::vector<int>, std::vector<var>,
                   double, empty>
    type_vv_real_real_int_real_real_1869;
typedef std::tuple<var, std::vector<double>, std::vector<int>, std::vector<var>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_1870;
typedef std::tuple<var, std::vector<double>, std::vector<int>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1871;
typedef std::tuple<var, std::vector<double>, std::vector<int>, std::vector<var>,
                   var, empty>
    type_vv_real_real_int_real_real_1872;
typedef std::tuple<var, std::vector<double>, std::vector<int>, std::vector<var>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_1873;
typedef std::tuple<
    var, std::vector<double>, std::vector<int>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1874;
typedef std::tuple<
    var, std::vector<double>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_1875;
typedef std::tuple<
    var, std::vector<double>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_1876;
typedef std::tuple<
    var, std::vector<double>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1877;
typedef std::tuple<
    var, std::vector<double>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_1878;
typedef std::tuple<
    var, std::vector<double>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_1879;
typedef std::tuple<
    var, std::vector<double>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1880;
typedef std::tuple<var, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double, double, empty>
    type_vv_real_real_int_real_real_1881;
typedef std::tuple<var, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_1882;
typedef std::tuple<var, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1883;
typedef std::tuple<var, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double, var, empty>
    type_vv_real_real_int_real_real_1884;
typedef std::tuple<var, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_1885;
typedef std::tuple<
    var, std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1886;
typedef std::tuple<var, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   double, empty>
    type_vv_real_real_int_real_real_1887;
typedef std::tuple<var, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_1888;
typedef std::tuple<var, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1889;
typedef std::tuple<var, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   var, empty>
    type_vv_real_real_int_real_real_1890;
typedef std::tuple<var, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_1891;
typedef std::tuple<
    var, std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1892;
typedef std::tuple<var, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_1893;
typedef std::tuple<
    var, std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_1894;
typedef std::tuple<var, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1895;
typedef std::tuple<var, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_1896;
typedef std::tuple<
    var, std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1897;
typedef std::tuple<
    var, std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1898;
typedef std::tuple<var, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var, double, empty>
    type_vv_real_real_int_real_real_1899;
typedef std::tuple<var, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_1900;
typedef std::tuple<var, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1901;
typedef std::tuple<var, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var, var, empty>
    type_vv_real_real_int_real_real_1902;
typedef std::tuple<var, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_1903;
typedef std::tuple<
    var, std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1904;
typedef std::tuple<var, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   double, empty>
    type_vv_real_real_int_real_real_1905;
typedef std::tuple<var, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_1906;
typedef std::tuple<var, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1907;
typedef std::tuple<var, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>, var,
                   empty>
    type_vv_real_real_int_real_real_1908;
typedef std::tuple<var, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_1909;
typedef std::tuple<
    var, std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1910;
typedef std::tuple<
    var, std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_1911;
typedef std::tuple<
    var, std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_1912;
typedef std::tuple<
    var, std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1913;
typedef std::tuple<
    var, std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_1914;
typedef std::tuple<
    var, std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_1915;
typedef std::tuple<
    var, std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1916;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>, int, double,
                   double, empty>
    type_vv_real_real_int_real_real_1917;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>, int, double,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_1918;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>, int, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1919;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>, int, double,
                   var, empty>
    type_vv_real_real_int_real_real_1920;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>, int, double,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_1921;
typedef std::tuple<
    var, Eigen::Matrix<double, Eigen::Dynamic, 1>, int, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1922;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
                   std::vector<double>, double, empty>
    type_vv_real_real_int_real_real_1923;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
                   std::vector<double>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_1924;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
                   std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1925;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
                   std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_1926;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
                   std::vector<double>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1927;
typedef std::tuple<
    var, Eigen::Matrix<double, Eigen::Dynamic, 1>, int, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1928;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_1929;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_1930;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1931;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_1932;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_1933;
typedef std::tuple<
    var, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1934;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>, int, var,
                   double, empty>
    type_vv_real_real_int_real_real_1935;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>, int, var,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_1936;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>, int, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1937;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>, int, var, var,
                   empty>
    type_vv_real_real_int_real_real_1938;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>, int, var,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_1939;
typedef std::tuple<
    var, Eigen::Matrix<double, Eigen::Dynamic, 1>, int, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1940;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
                   std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_1941;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
                   std::vector<var>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_1942;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
                   std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   empty>
    type_vv_real_real_int_real_real_1943;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
                   std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_1944;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
                   std::vector<var>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1945;
typedef std::tuple<
    var, Eigen::Matrix<double, Eigen::Dynamic, 1>, int, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1946;
typedef std::tuple<
    var, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_1947;
typedef std::tuple<
    var, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_1948;
typedef std::tuple<
    var, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1949;
typedef std::tuple<
    var, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_1950;
typedef std::tuple<
    var, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_1951;
typedef std::tuple<
    var, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1952;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, double, double, empty>
    type_vv_real_real_int_real_real_1953;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, double, std::vector<double>, empty>
    type_vv_real_real_int_real_real_1954;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1955;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, double, var, empty>
    type_vv_real_real_int_real_real_1956;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, double, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1957;
typedef std::tuple<
    var, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1958;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, std::vector<double>, double, empty>
    type_vv_real_real_int_real_real_1959;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, std::vector<double>, std::vector<double>,
                   empty>
    type_vv_real_real_int_real_real_1960;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1961;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_1962;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, std::vector<double>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_1963;
typedef std::tuple<
    var, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1964;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   double, empty>
    type_vv_real_real_int_real_real_1965;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_1966;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1967;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   var, empty>
    type_vv_real_real_int_real_real_1968;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_1969;
typedef std::tuple<
    var, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1970;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, var, double, empty>
    type_vv_real_real_int_real_real_1971;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, var, std::vector<double>, empty>
    type_vv_real_real_int_real_real_1972;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1973;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, var, var, empty>
    type_vv_real_real_int_real_real_1974;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, var, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1975;
typedef std::tuple<
    var, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1976;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_1977;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, std::vector<var>, std::vector<double>,
                   empty>
    type_vv_real_real_int_real_real_1978;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1979;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_1980;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, std::vector<var>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_1981;
typedef std::tuple<
    var, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1982;
typedef std::tuple<
    var, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_1983;
typedef std::tuple<
    var, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_1984;
typedef std::tuple<
    var, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1985;
typedef std::tuple<
    var, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_1986;
typedef std::tuple<
    var, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_1987;
typedef std::tuple<
    var, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1988;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double, double, empty>
    type_vv_real_real_int_real_real_1989;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_1990;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1991;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double, var, empty>
    type_vv_real_real_int_real_real_1992;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_1993;
typedef std::tuple<
    var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_1994;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   double, empty>
    type_vv_real_real_int_real_real_1995;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_1996;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_1997;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   var, empty>
    type_vv_real_real_int_real_real_1998;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_1999;
typedef std::tuple<
    var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2000;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_2001;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2002;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2003;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_2004;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_2005;
typedef std::tuple<
    var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2006;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var, double, empty>
    type_vv_real_real_int_real_real_2007;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2008;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2009;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var, var, empty>
    type_vv_real_real_int_real_real_2010;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_2011;
typedef std::tuple<
    var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2012;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   double, empty>
    type_vv_real_real_int_real_real_2013;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2014;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2015;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>, var,
                   empty>
    type_vv_real_real_int_real_real_2016;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_2017;
typedef std::tuple<
    var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2018;
typedef std::tuple<
    var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_2019;
typedef std::tuple<
    var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_2020;
typedef std::tuple<
    var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2021;
typedef std::tuple<
    var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_2022;
typedef std::tuple<
    var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_2023;
typedef std::tuple<
    var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2024;
typedef std::tuple<var, var, int, double, double, empty>
    type_vv_real_real_int_real_real_2025;
typedef std::tuple<var, var, int, double, std::vector<double>, empty>
    type_vv_real_real_int_real_real_2026;
typedef std::tuple<var, var, int, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2027;
typedef std::tuple<var, var, int, double, var, empty>
    type_vv_real_real_int_real_real_2028;
typedef std::tuple<var, var, int, double, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2029;
typedef std::tuple<
    var, var, int, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2030;
typedef std::tuple<var, var, int, std::vector<double>, double, empty>
    type_vv_real_real_int_real_real_2031;
typedef std::tuple<var, var, int, std::vector<double>, std::vector<double>,
                   empty>
    type_vv_real_real_int_real_real_2032;
typedef std::tuple<var, var, int, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2033;
typedef std::tuple<var, var, int, std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_2034;
typedef std::tuple<var, var, int, std::vector<double>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2035;
typedef std::tuple<
    var, var, int, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2036;
typedef std::tuple<var, var, int, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   double, empty>
    type_vv_real_real_int_real_real_2037;
typedef std::tuple<var, var, int, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2038;
typedef std::tuple<var, var, int, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2039;
typedef std::tuple<var, var, int, Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   empty>
    type_vv_real_real_int_real_real_2040;
typedef std::tuple<var, var, int, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_2041;
typedef std::tuple<
    var, var, int, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2042;
typedef std::tuple<var, var, int, var, double, empty>
    type_vv_real_real_int_real_real_2043;
typedef std::tuple<var, var, int, var, std::vector<double>, empty>
    type_vv_real_real_int_real_real_2044;
typedef std::tuple<var, var, int, var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   empty>
    type_vv_real_real_int_real_real_2045;
typedef std::tuple<var, var, int, var, var, empty>
    type_vv_real_real_int_real_real_2046;
typedef std::tuple<var, var, int, var, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2047;
typedef std::tuple<
    var, var, int, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2048;
typedef std::tuple<var, var, int, std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_2049;
typedef std::tuple<var, var, int, std::vector<var>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_2050;
typedef std::tuple<var, var, int, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2051;
typedef std::tuple<var, var, int, std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_2052;
typedef std::tuple<var, var, int, std::vector<var>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2053;
typedef std::tuple<
    var, var, int, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2054;
typedef std::tuple<
    var, var, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_2055;
typedef std::tuple<
    var, var, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_2056;
typedef std::tuple<
    var, var, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2057;
typedef std::tuple<
    var, var, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_2058;
typedef std::tuple<
    var, var, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_2059;
typedef std::tuple<
    var, var, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2060;
typedef std::tuple<var, var, std::vector<int>, double, double, empty>
    type_vv_real_real_int_real_real_2061;
typedef std::tuple<var, var, std::vector<int>, double, std::vector<double>,
                   empty>
    type_vv_real_real_int_real_real_2062;
typedef std::tuple<var, var, std::vector<int>, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2063;
typedef std::tuple<var, var, std::vector<int>, double, var, empty>
    type_vv_real_real_int_real_real_2064;
typedef std::tuple<var, var, std::vector<int>, double, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2065;
typedef std::tuple<
    var, var, std::vector<int>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2066;
typedef std::tuple<var, var, std::vector<int>, std::vector<double>, double,
                   empty>
    type_vv_real_real_int_real_real_2067;
typedef std::tuple<var, var, std::vector<int>, std::vector<double>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2068;
typedef std::tuple<var, var, std::vector<int>, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2069;
typedef std::tuple<var, var, std::vector<int>, std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_2070;
typedef std::tuple<var, var, std::vector<int>, std::vector<double>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_2071;
typedef std::tuple<
    var, var, std::vector<int>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2072;
typedef std::tuple<var, var, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_2073;
typedef std::tuple<var, var, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2074;
typedef std::tuple<var, var, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2075;
typedef std::tuple<var, var, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_2076;
typedef std::tuple<var, var, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_2077;
typedef std::tuple<
    var, var, std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2078;
typedef std::tuple<var, var, std::vector<int>, var, double, empty>
    type_vv_real_real_int_real_real_2079;
typedef std::tuple<var, var, std::vector<int>, var, std::vector<double>, empty>
    type_vv_real_real_int_real_real_2080;
typedef std::tuple<var, var, std::vector<int>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2081;
typedef std::tuple<var, var, std::vector<int>, var, var, empty>
    type_vv_real_real_int_real_real_2082;
typedef std::tuple<var, var, std::vector<int>, var, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2083;
typedef std::tuple<
    var, var, std::vector<int>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2084;
typedef std::tuple<var, var, std::vector<int>, std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_2085;
typedef std::tuple<var, var, std::vector<int>, std::vector<var>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2086;
typedef std::tuple<var, var, std::vector<int>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2087;
typedef std::tuple<var, var, std::vector<int>, std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_2088;
typedef std::tuple<var, var, std::vector<int>, std::vector<var>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_2089;
typedef std::tuple<
    var, var, std::vector<int>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2090;
typedef std::tuple<
    var, var, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_2091;
typedef std::tuple<
    var, var, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_2092;
typedef std::tuple<
    var, var, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2093;
typedef std::tuple<
    var, var, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_2094;
typedef std::tuple<
    var, var, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_2095;
typedef std::tuple<
    var, var, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2096;
typedef std::tuple<var, var, Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
                   double, empty>
    type_vv_real_real_int_real_real_2097;
typedef std::tuple<var, var, Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2098;
typedef std::tuple<var, var, Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2099;
typedef std::tuple<var, var, Eigen::Matrix<int, Eigen::Dynamic, 1>, double, var,
                   empty>
    type_vv_real_real_int_real_real_2100;
typedef std::tuple<var, var, Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_2101;
typedef std::tuple<
    var, var, Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2102;
typedef std::tuple<var, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<double>, double, empty>
    type_vv_real_real_int_real_real_2103;
typedef std::tuple<var, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<double>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_2104;
typedef std::tuple<var, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2105;
typedef std::tuple<var, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_2106;
typedef std::tuple<var, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<double>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2107;
typedef std::tuple<
    var, var, Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2108;
typedef std::tuple<var, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_2109;
typedef std::tuple<var, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2110;
typedef std::tuple<var, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2111;
typedef std::tuple<var, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_2112;
typedef std::tuple<var, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_2113;
typedef std::tuple<
    var, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2114;
typedef std::tuple<var, var, Eigen::Matrix<int, Eigen::Dynamic, 1>, var, double,
                   empty>
    type_vv_real_real_int_real_real_2115;
typedef std::tuple<var, var, Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2116;
typedef std::tuple<var, var, Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2117;
typedef std::tuple<var, var, Eigen::Matrix<int, Eigen::Dynamic, 1>, var, var,
                   empty>
    type_vv_real_real_int_real_real_2118;
typedef std::tuple<var, var, Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_2119;
typedef std::tuple<
    var, var, Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2120;
typedef std::tuple<var, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_2121;
typedef std::tuple<var, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<var>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_2122;
typedef std::tuple<var, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   empty>
    type_vv_real_real_int_real_real_2123;
typedef std::tuple<var, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_2124;
typedef std::tuple<var, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<var>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2125;
typedef std::tuple<
    var, var, Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2126;
typedef std::tuple<
    var, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_2127;
typedef std::tuple<
    var, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_2128;
typedef std::tuple<
    var, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2129;
typedef std::tuple<
    var, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_2130;
typedef std::tuple<
    var, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_2131;
typedef std::tuple<
    var, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2132;
typedef std::tuple<var, std::vector<var>, int, double, double, empty>
    type_vv_real_real_int_real_real_2133;
typedef std::tuple<var, std::vector<var>, int, double, std::vector<double>,
                   empty>
    type_vv_real_real_int_real_real_2134;
typedef std::tuple<var, std::vector<var>, int, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2135;
typedef std::tuple<var, std::vector<var>, int, double, var, empty>
    type_vv_real_real_int_real_real_2136;
typedef std::tuple<var, std::vector<var>, int, double, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2137;
typedef std::tuple<
    var, std::vector<var>, int, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2138;
typedef std::tuple<var, std::vector<var>, int, std::vector<double>, double,
                   empty>
    type_vv_real_real_int_real_real_2139;
typedef std::tuple<var, std::vector<var>, int, std::vector<double>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2140;
typedef std::tuple<var, std::vector<var>, int, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2141;
typedef std::tuple<var, std::vector<var>, int, std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_2142;
typedef std::tuple<var, std::vector<var>, int, std::vector<double>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_2143;
typedef std::tuple<
    var, std::vector<var>, int, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2144;
typedef std::tuple<var, std::vector<var>, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_2145;
typedef std::tuple<var, std::vector<var>, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2146;
typedef std::tuple<var, std::vector<var>, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2147;
typedef std::tuple<var, std::vector<var>, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_2148;
typedef std::tuple<var, std::vector<var>, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_2149;
typedef std::tuple<
    var, std::vector<var>, int, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2150;
typedef std::tuple<var, std::vector<var>, int, var, double, empty>
    type_vv_real_real_int_real_real_2151;
typedef std::tuple<var, std::vector<var>, int, var, std::vector<double>, empty>
    type_vv_real_real_int_real_real_2152;
typedef std::tuple<var, std::vector<var>, int, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2153;
typedef std::tuple<var, std::vector<var>, int, var, var, empty>
    type_vv_real_real_int_real_real_2154;
typedef std::tuple<var, std::vector<var>, int, var, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2155;
typedef std::tuple<
    var, std::vector<var>, int, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2156;
typedef std::tuple<var, std::vector<var>, int, std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_2157;
typedef std::tuple<var, std::vector<var>, int, std::vector<var>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2158;
typedef std::tuple<var, std::vector<var>, int, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2159;
typedef std::tuple<var, std::vector<var>, int, std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_2160;
typedef std::tuple<var, std::vector<var>, int, std::vector<var>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_2161;
typedef std::tuple<
    var, std::vector<var>, int, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2162;
typedef std::tuple<
    var, std::vector<var>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_2163;
typedef std::tuple<
    var, std::vector<var>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_2164;
typedef std::tuple<
    var, std::vector<var>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2165;
typedef std::tuple<
    var, std::vector<var>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_2166;
typedef std::tuple<
    var, std::vector<var>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_2167;
typedef std::tuple<
    var, std::vector<var>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2168;
typedef std::tuple<var, std::vector<var>, std::vector<int>, double, double,
                   empty>
    type_vv_real_real_int_real_real_2169;
typedef std::tuple<var, std::vector<var>, std::vector<int>, double,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2170;
typedef std::tuple<var, std::vector<var>, std::vector<int>, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2171;
typedef std::tuple<var, std::vector<var>, std::vector<int>, double, var, empty>
    type_vv_real_real_int_real_real_2172;
typedef std::tuple<var, std::vector<var>, std::vector<int>, double,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_2173;
typedef std::tuple<
    var, std::vector<var>, std::vector<int>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2174;
typedef std::tuple<var, std::vector<var>, std::vector<int>, std::vector<double>,
                   double, empty>
    type_vv_real_real_int_real_real_2175;
typedef std::tuple<var, std::vector<var>, std::vector<int>, std::vector<double>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2176;
typedef std::tuple<var, std::vector<var>, std::vector<int>, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2177;
typedef std::tuple<var, std::vector<var>, std::vector<int>, std::vector<double>,
                   var, empty>
    type_vv_real_real_int_real_real_2178;
typedef std::tuple<var, std::vector<var>, std::vector<int>, std::vector<double>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_2179;
typedef std::tuple<
    var, std::vector<var>, std::vector<int>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2180;
typedef std::tuple<var, std::vector<var>, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_2181;
typedef std::tuple<var, std::vector<var>, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2182;
typedef std::tuple<var, std::vector<var>, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2183;
typedef std::tuple<var, std::vector<var>, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_2184;
typedef std::tuple<var, std::vector<var>, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_2185;
typedef std::tuple<
    var, std::vector<var>, std::vector<int>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2186;
typedef std::tuple<var, std::vector<var>, std::vector<int>, var, double, empty>
    type_vv_real_real_int_real_real_2187;
typedef std::tuple<var, std::vector<var>, std::vector<int>, var,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2188;
typedef std::tuple<var, std::vector<var>, std::vector<int>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2189;
typedef std::tuple<var, std::vector<var>, std::vector<int>, var, var, empty>
    type_vv_real_real_int_real_real_2190;
typedef std::tuple<var, std::vector<var>, std::vector<int>, var,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_2191;
typedef std::tuple<
    var, std::vector<var>, std::vector<int>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2192;
typedef std::tuple<var, std::vector<var>, std::vector<int>, std::vector<var>,
                   double, empty>
    type_vv_real_real_int_real_real_2193;
typedef std::tuple<var, std::vector<var>, std::vector<int>, std::vector<var>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2194;
typedef std::tuple<var, std::vector<var>, std::vector<int>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2195;
typedef std::tuple<var, std::vector<var>, std::vector<int>, std::vector<var>,
                   var, empty>
    type_vv_real_real_int_real_real_2196;
typedef std::tuple<var, std::vector<var>, std::vector<int>, std::vector<var>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_2197;
typedef std::tuple<
    var, std::vector<var>, std::vector<int>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2198;
typedef std::tuple<
    var, std::vector<var>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_2199;
typedef std::tuple<
    var, std::vector<var>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_2200;
typedef std::tuple<
    var, std::vector<var>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2201;
typedef std::tuple<
    var, std::vector<var>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_2202;
typedef std::tuple<
    var, std::vector<var>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_2203;
typedef std::tuple<
    var, std::vector<var>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2204;
typedef std::tuple<var, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   double, double, empty>
    type_vv_real_real_int_real_real_2205;
typedef std::tuple<var, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   double, std::vector<double>, empty>
    type_vv_real_real_int_real_real_2206;
typedef std::tuple<var, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   double, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2207;
typedef std::tuple<var, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   double, var, empty>
    type_vv_real_real_int_real_real_2208;
typedef std::tuple<var, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   double, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2209;
typedef std::tuple<
    var, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2210;
typedef std::tuple<var, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<double>, double, empty>
    type_vv_real_real_int_real_real_2211;
typedef std::tuple<var, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<double>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_2212;
typedef std::tuple<var, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2213;
typedef std::tuple<var, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_2214;
typedef std::tuple<var, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<double>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2215;
typedef std::tuple<
    var, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2216;
typedef std::tuple<var, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_2217;
typedef std::tuple<var, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2218;
typedef std::tuple<var, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2219;
typedef std::tuple<var, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_2220;
typedef std::tuple<var, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_2221;
typedef std::tuple<
    var, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2222;
typedef std::tuple<var, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   var, double, empty>
    type_vv_real_real_int_real_real_2223;
typedef std::tuple<var, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   var, std::vector<double>, empty>
    type_vv_real_real_int_real_real_2224;
typedef std::tuple<var, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   var, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2225;
typedef std::tuple<var, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   var, var, empty>
    type_vv_real_real_int_real_real_2226;
typedef std::tuple<var, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   var, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2227;
typedef std::tuple<
    var, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2228;
typedef std::tuple<var, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_2229;
typedef std::tuple<var, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<var>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_2230;
typedef std::tuple<var, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   empty>
    type_vv_real_real_int_real_real_2231;
typedef std::tuple<var, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_2232;
typedef std::tuple<var, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<var>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2233;
typedef std::tuple<
    var, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2234;
typedef std::tuple<
    var, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_2235;
typedef std::tuple<
    var, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_2236;
typedef std::tuple<
    var, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2237;
typedef std::tuple<
    var, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_2238;
typedef std::tuple<
    var, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_2239;
typedef std::tuple<
    var, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2240;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    double, double, empty>
    type_vv_real_real_int_real_real_2241;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    double, std::vector<double>, empty>
    type_vv_real_real_int_real_real_2242;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    double, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2243;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    double, var, empty>
    type_vv_real_real_int_real_real_2244;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    double, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2245;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    empty>
    type_vv_real_real_int_real_real_2246;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<double>, double, empty>
    type_vv_real_real_int_real_real_2247;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<double>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_2248;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2249;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_2250;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<double>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2251;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2252;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_2253;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_2254;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2255;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_2256;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2257;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2258;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    var, double, empty>
    type_vv_real_real_int_real_real_2259;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    var, std::vector<double>, empty>
    type_vv_real_real_int_real_real_2260;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    var, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2261;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    var, var, empty>
    type_vv_real_real_int_real_real_2262;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    var, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2263;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2264;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_2265;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<var>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_2266;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2267;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_2268;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<var>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2269;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2270;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_2271;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_2272;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2273;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_2274;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_2275;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2276;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, double, double, empty>
    type_vv_real_real_int_real_real_2277;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, double, std::vector<double>, empty>
    type_vv_real_real_int_real_real_2278;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, double, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2279;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, double, var, empty>
    type_vv_real_real_int_real_real_2280;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, double, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2281;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2282;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<double>, double, empty>
    type_vv_real_real_int_real_real_2283;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<double>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_2284;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<double>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2285;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_2286;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<double>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2287;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2288;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_2289;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_2290;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2291;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_2292;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_2293;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2294;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, var, double, empty>
    type_vv_real_real_int_real_real_2295;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, var, std::vector<double>, empty>
    type_vv_real_real_int_real_real_2296;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, var, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2297;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, var, var, empty>
    type_vv_real_real_int_real_real_2298;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, var, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2299;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2300;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_2301;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<var>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_2302;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<var>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2303;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_2304;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<var>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2305;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2306;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_2307;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_2308;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2309;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_2310;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_2311;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2312;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double, double, empty>
    type_vv_real_real_int_real_real_2313;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double, std::vector<double>, empty>
    type_vv_real_real_int_real_real_2314;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2315;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double, var, empty>
    type_vv_real_real_int_real_real_2316;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2317;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2318;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>, double, empty>
    type_vv_real_real_int_real_real_2319;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_2320;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2321;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_2322;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_2323;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2324;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_2325;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_2326;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2327;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_2328;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2329;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2330;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var, double, empty>
    type_vv_real_real_int_real_real_2331;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var, std::vector<double>, empty>
    type_vv_real_real_int_real_real_2332;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2333;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var, var, empty>
    type_vv_real_real_int_real_real_2334;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2335;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2336;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_2337;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_2338;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2339;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_2340;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>, std::vector<var>,
    empty>
    type_vv_real_real_int_real_real_2341;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2342;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_2343;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_2344;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2345;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_2346;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_2347;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2348;
typedef std::tuple<std::vector<var>, double, int, double, double, empty>
    type_vv_real_real_int_real_real_2349;
typedef std::tuple<std::vector<var>, double, int, double, std::vector<double>,
                   empty>
    type_vv_real_real_int_real_real_2350;
typedef std::tuple<std::vector<var>, double, int, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2351;
typedef std::tuple<std::vector<var>, double, int, double, var, empty>
    type_vv_real_real_int_real_real_2352;
typedef std::tuple<std::vector<var>, double, int, double, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_2353;
typedef std::tuple<
    std::vector<var>, double, int, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2354;
typedef std::tuple<std::vector<var>, double, int, std::vector<double>, double,
                   empty>
    type_vv_real_real_int_real_real_2355;
typedef std::tuple<std::vector<var>, double, int, std::vector<double>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2356;
typedef std::tuple<std::vector<var>, double, int, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2357;
typedef std::tuple<std::vector<var>, double, int, std::vector<double>, var,
                   empty>
    type_vv_real_real_int_real_real_2358;
typedef std::tuple<std::vector<var>, double, int, std::vector<double>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_2359;
typedef std::tuple<
    std::vector<var>, double, int, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2360;
typedef std::tuple<std::vector<var>, double, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_2361;
typedef std::tuple<std::vector<var>, double, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2362;
typedef std::tuple<std::vector<var>, double, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2363;
typedef std::tuple<std::vector<var>, double, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_2364;
typedef std::tuple<std::vector<var>, double, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_2365;
typedef std::tuple<
    std::vector<var>, double, int, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2366;
typedef std::tuple<std::vector<var>, double, int, var, double, empty>
    type_vv_real_real_int_real_real_2367;
typedef std::tuple<std::vector<var>, double, int, var, std::vector<double>,
                   empty>
    type_vv_real_real_int_real_real_2368;
typedef std::tuple<std::vector<var>, double, int, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2369;
typedef std::tuple<std::vector<var>, double, int, var, var, empty>
    type_vv_real_real_int_real_real_2370;
typedef std::tuple<std::vector<var>, double, int, var, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2371;
typedef std::tuple<
    std::vector<var>, double, int, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2372;
typedef std::tuple<std::vector<var>, double, int, std::vector<var>, double,
                   empty>
    type_vv_real_real_int_real_real_2373;
typedef std::tuple<std::vector<var>, double, int, std::vector<var>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2374;
typedef std::tuple<std::vector<var>, double, int, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2375;
typedef std::tuple<std::vector<var>, double, int, std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_2376;
typedef std::tuple<std::vector<var>, double, int, std::vector<var>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_2377;
typedef std::tuple<
    std::vector<var>, double, int, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2378;
typedef std::tuple<
    std::vector<var>, double, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_2379;
typedef std::tuple<
    std::vector<var>, double, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_2380;
typedef std::tuple<
    std::vector<var>, double, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2381;
typedef std::tuple<
    std::vector<var>, double, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_2382;
typedef std::tuple<
    std::vector<var>, double, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_2383;
typedef std::tuple<
    std::vector<var>, double, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2384;
typedef std::tuple<std::vector<var>, double, std::vector<int>, double, double,
                   empty>
    type_vv_real_real_int_real_real_2385;
typedef std::tuple<std::vector<var>, double, std::vector<int>, double,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2386;
typedef std::tuple<std::vector<var>, double, std::vector<int>, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2387;
typedef std::tuple<std::vector<var>, double, std::vector<int>, double, var,
                   empty>
    type_vv_real_real_int_real_real_2388;
typedef std::tuple<std::vector<var>, double, std::vector<int>, double,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_2389;
typedef std::tuple<
    std::vector<var>, double, std::vector<int>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2390;
typedef std::tuple<std::vector<var>, double, std::vector<int>,
                   std::vector<double>, double, empty>
    type_vv_real_real_int_real_real_2391;
typedef std::tuple<std::vector<var>, double, std::vector<int>,
                   std::vector<double>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_2392;
typedef std::tuple<std::vector<var>, double, std::vector<int>,
                   std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2393;
typedef std::tuple<std::vector<var>, double, std::vector<int>,
                   std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_2394;
typedef std::tuple<std::vector<var>, double, std::vector<int>,
                   std::vector<double>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2395;
typedef std::tuple<
    std::vector<var>, double, std::vector<int>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2396;
typedef std::tuple<std::vector<var>, double, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_2397;
typedef std::tuple<std::vector<var>, double, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2398;
typedef std::tuple<std::vector<var>, double, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2399;
typedef std::tuple<std::vector<var>, double, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_2400;
typedef std::tuple<std::vector<var>, double, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_2401;
typedef std::tuple<
    std::vector<var>, double, std::vector<int>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2402;
typedef std::tuple<std::vector<var>, double, std::vector<int>, var, double,
                   empty>
    type_vv_real_real_int_real_real_2403;
typedef std::tuple<std::vector<var>, double, std::vector<int>, var,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2404;
typedef std::tuple<std::vector<var>, double, std::vector<int>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2405;
typedef std::tuple<std::vector<var>, double, std::vector<int>, var, var, empty>
    type_vv_real_real_int_real_real_2406;
typedef std::tuple<std::vector<var>, double, std::vector<int>, var,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_2407;
typedef std::tuple<
    std::vector<var>, double, std::vector<int>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2408;
typedef std::tuple<std::vector<var>, double, std::vector<int>, std::vector<var>,
                   double, empty>
    type_vv_real_real_int_real_real_2409;
typedef std::tuple<std::vector<var>, double, std::vector<int>, std::vector<var>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2410;
typedef std::tuple<std::vector<var>, double, std::vector<int>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2411;
typedef std::tuple<std::vector<var>, double, std::vector<int>, std::vector<var>,
                   var, empty>
    type_vv_real_real_int_real_real_2412;
typedef std::tuple<std::vector<var>, double, std::vector<int>, std::vector<var>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_2413;
typedef std::tuple<
    std::vector<var>, double, std::vector<int>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2414;
typedef std::tuple<
    std::vector<var>, double, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_2415;
typedef std::tuple<
    std::vector<var>, double, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_2416;
typedef std::tuple<
    std::vector<var>, double, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2417;
typedef std::tuple<
    std::vector<var>, double, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_2418;
typedef std::tuple<
    std::vector<var>, double, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_2419;
typedef std::tuple<
    std::vector<var>, double, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2420;
typedef std::tuple<std::vector<var>, double,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double, double, empty>
    type_vv_real_real_int_real_real_2421;
typedef std::tuple<std::vector<var>, double,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2422;
typedef std::tuple<std::vector<var>, double,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2423;
typedef std::tuple<std::vector<var>, double,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double, var, empty>
    type_vv_real_real_int_real_real_2424;
typedef std::tuple<std::vector<var>, double,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_2425;
typedef std::tuple<
    std::vector<var>, double, Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2426;
typedef std::tuple<std::vector<var>, double,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   double, empty>
    type_vv_real_real_int_real_real_2427;
typedef std::tuple<std::vector<var>, double,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2428;
typedef std::tuple<std::vector<var>, double,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2429;
typedef std::tuple<std::vector<var>, double,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   var, empty>
    type_vv_real_real_int_real_real_2430;
typedef std::tuple<std::vector<var>, double,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_2431;
typedef std::tuple<
    std::vector<var>, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2432;
typedef std::tuple<std::vector<var>, double,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_2433;
typedef std::tuple<
    std::vector<var>, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_2434;
typedef std::tuple<std::vector<var>, double,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2435;
typedef std::tuple<std::vector<var>, double,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_2436;
typedef std::tuple<
    std::vector<var>, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2437;
typedef std::tuple<
    std::vector<var>, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2438;
typedef std::tuple<std::vector<var>, double,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var, double, empty>
    type_vv_real_real_int_real_real_2439;
typedef std::tuple<std::vector<var>, double,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2440;
typedef std::tuple<std::vector<var>, double,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2441;
typedef std::tuple<std::vector<var>, double,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var, var, empty>
    type_vv_real_real_int_real_real_2442;
typedef std::tuple<std::vector<var>, double,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_2443;
typedef std::tuple<
    std::vector<var>, double, Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2444;
typedef std::tuple<std::vector<var>, double,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   double, empty>
    type_vv_real_real_int_real_real_2445;
typedef std::tuple<std::vector<var>, double,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2446;
typedef std::tuple<std::vector<var>, double,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2447;
typedef std::tuple<std::vector<var>, double,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>, var,
                   empty>
    type_vv_real_real_int_real_real_2448;
typedef std::tuple<std::vector<var>, double,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_2449;
typedef std::tuple<
    std::vector<var>, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2450;
typedef std::tuple<
    std::vector<var>, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_2451;
typedef std::tuple<
    std::vector<var>, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_2452;
typedef std::tuple<
    std::vector<var>, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2453;
typedef std::tuple<
    std::vector<var>, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_2454;
typedef std::tuple<
    std::vector<var>, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_2455;
typedef std::tuple<
    std::vector<var>, double, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2456;
typedef std::tuple<std::vector<var>, std::vector<double>, int, double, double,
                   empty>
    type_vv_real_real_int_real_real_2457;
typedef std::tuple<std::vector<var>, std::vector<double>, int, double,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2458;
typedef std::tuple<std::vector<var>, std::vector<double>, int, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2459;
typedef std::tuple<std::vector<var>, std::vector<double>, int, double, var,
                   empty>
    type_vv_real_real_int_real_real_2460;
typedef std::tuple<std::vector<var>, std::vector<double>, int, double,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_2461;
typedef std::tuple<
    std::vector<var>, std::vector<double>, int, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2462;
typedef std::tuple<std::vector<var>, std::vector<double>, int,
                   std::vector<double>, double, empty>
    type_vv_real_real_int_real_real_2463;
typedef std::tuple<std::vector<var>, std::vector<double>, int,
                   std::vector<double>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_2464;
typedef std::tuple<std::vector<var>, std::vector<double>, int,
                   std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2465;
typedef std::tuple<std::vector<var>, std::vector<double>, int,
                   std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_2466;
typedef std::tuple<std::vector<var>, std::vector<double>, int,
                   std::vector<double>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2467;
typedef std::tuple<
    std::vector<var>, std::vector<double>, int, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2468;
typedef std::tuple<std::vector<var>, std::vector<double>, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_2469;
typedef std::tuple<std::vector<var>, std::vector<double>, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2470;
typedef std::tuple<std::vector<var>, std::vector<double>, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2471;
typedef std::tuple<std::vector<var>, std::vector<double>, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_2472;
typedef std::tuple<std::vector<var>, std::vector<double>, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_2473;
typedef std::tuple<
    std::vector<var>, std::vector<double>, int,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2474;
typedef std::tuple<std::vector<var>, std::vector<double>, int, var, double,
                   empty>
    type_vv_real_real_int_real_real_2475;
typedef std::tuple<std::vector<var>, std::vector<double>, int, var,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2476;
typedef std::tuple<std::vector<var>, std::vector<double>, int, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2477;
typedef std::tuple<std::vector<var>, std::vector<double>, int, var, var, empty>
    type_vv_real_real_int_real_real_2478;
typedef std::tuple<std::vector<var>, std::vector<double>, int, var,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_2479;
typedef std::tuple<
    std::vector<var>, std::vector<double>, int, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2480;
typedef std::tuple<std::vector<var>, std::vector<double>, int, std::vector<var>,
                   double, empty>
    type_vv_real_real_int_real_real_2481;
typedef std::tuple<std::vector<var>, std::vector<double>, int, std::vector<var>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2482;
typedef std::tuple<std::vector<var>, std::vector<double>, int, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2483;
typedef std::tuple<std::vector<var>, std::vector<double>, int, std::vector<var>,
                   var, empty>
    type_vv_real_real_int_real_real_2484;
typedef std::tuple<std::vector<var>, std::vector<double>, int, std::vector<var>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_2485;
typedef std::tuple<
    std::vector<var>, std::vector<double>, int, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2486;
typedef std::tuple<
    std::vector<var>, std::vector<double>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_2487;
typedef std::tuple<
    std::vector<var>, std::vector<double>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_2488;
typedef std::tuple<
    std::vector<var>, std::vector<double>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2489;
typedef std::tuple<
    std::vector<var>, std::vector<double>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_2490;
typedef std::tuple<
    std::vector<var>, std::vector<double>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_2491;
typedef std::tuple<
    std::vector<var>, std::vector<double>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2492;
typedef std::tuple<std::vector<var>, std::vector<double>, std::vector<int>,
                   double, double, empty>
    type_vv_real_real_int_real_real_2493;
typedef std::tuple<std::vector<var>, std::vector<double>, std::vector<int>,
                   double, std::vector<double>, empty>
    type_vv_real_real_int_real_real_2494;
typedef std::tuple<std::vector<var>, std::vector<double>, std::vector<int>,
                   double, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2495;
typedef std::tuple<std::vector<var>, std::vector<double>, std::vector<int>,
                   double, var, empty>
    type_vv_real_real_int_real_real_2496;
typedef std::tuple<std::vector<var>, std::vector<double>, std::vector<int>,
                   double, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2497;
typedef std::tuple<
    std::vector<var>, std::vector<double>, std::vector<int>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2498;
typedef std::tuple<std::vector<var>, std::vector<double>, std::vector<int>,
                   std::vector<double>, double, empty>
    type_vv_real_real_int_real_real_2499;
typedef std::tuple<std::vector<var>, std::vector<double>, std::vector<int>,
                   std::vector<double>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_2500;
typedef std::tuple<std::vector<var>, std::vector<double>, std::vector<int>,
                   std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2501;
typedef std::tuple<std::vector<var>, std::vector<double>, std::vector<int>,
                   std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_2502;
typedef std::tuple<std::vector<var>, std::vector<double>, std::vector<int>,
                   std::vector<double>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2503;
typedef std::tuple<
    std::vector<var>, std::vector<double>, std::vector<int>,
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2504;
typedef std::tuple<std::vector<var>, std::vector<double>, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_2505;
typedef std::tuple<std::vector<var>, std::vector<double>, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2506;
typedef std::tuple<std::vector<var>, std::vector<double>, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2507;
typedef std::tuple<std::vector<var>, std::vector<double>, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_2508;
typedef std::tuple<std::vector<var>, std::vector<double>, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_2509;
typedef std::tuple<
    std::vector<var>, std::vector<double>, std::vector<int>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2510;
typedef std::tuple<std::vector<var>, std::vector<double>, std::vector<int>, var,
                   double, empty>
    type_vv_real_real_int_real_real_2511;
typedef std::tuple<std::vector<var>, std::vector<double>, std::vector<int>, var,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2512;
typedef std::tuple<std::vector<var>, std::vector<double>, std::vector<int>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2513;
typedef std::tuple<std::vector<var>, std::vector<double>, std::vector<int>, var,
                   var, empty>
    type_vv_real_real_int_real_real_2514;
typedef std::tuple<std::vector<var>, std::vector<double>, std::vector<int>, var,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_2515;
typedef std::tuple<
    std::vector<var>, std::vector<double>, std::vector<int>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2516;
typedef std::tuple<std::vector<var>, std::vector<double>, std::vector<int>,
                   std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_2517;
typedef std::tuple<std::vector<var>, std::vector<double>, std::vector<int>,
                   std::vector<var>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_2518;
typedef std::tuple<std::vector<var>, std::vector<double>, std::vector<int>,
                   std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   empty>
    type_vv_real_real_int_real_real_2519;
typedef std::tuple<std::vector<var>, std::vector<double>, std::vector<int>,
                   std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_2520;
typedef std::tuple<std::vector<var>, std::vector<double>, std::vector<int>,
                   std::vector<var>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2521;
typedef std::tuple<
    std::vector<var>, std::vector<double>, std::vector<int>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2522;
typedef std::tuple<
    std::vector<var>, std::vector<double>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_2523;
typedef std::tuple<
    std::vector<var>, std::vector<double>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_2524;
typedef std::tuple<
    std::vector<var>, std::vector<double>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2525;
typedef std::tuple<
    std::vector<var>, std::vector<double>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_2526;
typedef std::tuple<
    std::vector<var>, std::vector<double>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_2527;
typedef std::tuple<
    std::vector<var>, std::vector<double>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2528;
typedef std::tuple<std::vector<var>, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double, double, empty>
    type_vv_real_real_int_real_real_2529;
typedef std::tuple<std::vector<var>, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2530;
typedef std::tuple<std::vector<var>, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2531;
typedef std::tuple<std::vector<var>, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double, var, empty>
    type_vv_real_real_int_real_real_2532;
typedef std::tuple<std::vector<var>, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_2533;
typedef std::tuple<
    std::vector<var>, std::vector<double>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2534;
typedef std::tuple<std::vector<var>, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   double, empty>
    type_vv_real_real_int_real_real_2535;
typedef std::tuple<std::vector<var>, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2536;
typedef std::tuple<std::vector<var>, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2537;
typedef std::tuple<std::vector<var>, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   var, empty>
    type_vv_real_real_int_real_real_2538;
typedef std::tuple<std::vector<var>, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_2539;
typedef std::tuple<
    std::vector<var>, std::vector<double>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2540;
typedef std::tuple<std::vector<var>, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_2541;
typedef std::tuple<std::vector<var>, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2542;
typedef std::tuple<std::vector<var>, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2543;
typedef std::tuple<std::vector<var>, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_2544;
typedef std::tuple<std::vector<var>, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_2545;
typedef std::tuple<
    std::vector<var>, std::vector<double>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2546;
typedef std::tuple<std::vector<var>, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var, double, empty>
    type_vv_real_real_int_real_real_2547;
typedef std::tuple<std::vector<var>, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2548;
typedef std::tuple<std::vector<var>, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2549;
typedef std::tuple<std::vector<var>, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var, var, empty>
    type_vv_real_real_int_real_real_2550;
typedef std::tuple<std::vector<var>, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_2551;
typedef std::tuple<
    std::vector<var>, std::vector<double>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2552;
typedef std::tuple<std::vector<var>, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   double, empty>
    type_vv_real_real_int_real_real_2553;
typedef std::tuple<std::vector<var>, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2554;
typedef std::tuple<std::vector<var>, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2555;
typedef std::tuple<std::vector<var>, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>, var,
                   empty>
    type_vv_real_real_int_real_real_2556;
typedef std::tuple<std::vector<var>, std::vector<double>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_2557;
typedef std::tuple<
    std::vector<var>, std::vector<double>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2558;
typedef std::tuple<
    std::vector<var>, std::vector<double>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_2559;
typedef std::tuple<
    std::vector<var>, std::vector<double>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_2560;
typedef std::tuple<
    std::vector<var>, std::vector<double>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2561;
typedef std::tuple<
    std::vector<var>, std::vector<double>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_2562;
typedef std::tuple<
    std::vector<var>, std::vector<double>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_2563;
typedef std::tuple<
    std::vector<var>, std::vector<double>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2564;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   int, double, double, empty>
    type_vv_real_real_int_real_real_2565;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   int, double, std::vector<double>, empty>
    type_vv_real_real_int_real_real_2566;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   int, double, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2567;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   int, double, var, empty>
    type_vv_real_real_int_real_real_2568;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   int, double, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2569;
typedef std::tuple<
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>, int, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2570;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   int, std::vector<double>, double, empty>
    type_vv_real_real_int_real_real_2571;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   int, std::vector<double>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_2572;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   int, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2573;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   int, std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_2574;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   int, std::vector<double>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2575;
typedef std::tuple<
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2576;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   int, Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_2577;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   int, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2578;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   int, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2579;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   int, Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_2580;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   int, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_2581;
typedef std::tuple<
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2582;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   int, var, double, empty>
    type_vv_real_real_int_real_real_2583;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   int, var, std::vector<double>, empty>
    type_vv_real_real_int_real_real_2584;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   int, var, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2585;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   int, var, var, empty>
    type_vv_real_real_int_real_real_2586;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   int, var, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2587;
typedef std::tuple<
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>, int, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2588;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   int, std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_2589;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   int, std::vector<var>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_2590;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   int, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2591;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   int, std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_2592;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   int, std::vector<var>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2593;
typedef std::tuple<
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2594;
typedef std::tuple<
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_2595;
typedef std::tuple<
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_2596;
typedef std::tuple<
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2597;
typedef std::tuple<
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_2598;
typedef std::tuple<
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_2599;
typedef std::tuple<
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2600;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, double, double, empty>
    type_vv_real_real_int_real_real_2601;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, double, std::vector<double>, empty>
    type_vv_real_real_int_real_real_2602;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2603;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, double, var, empty>
    type_vv_real_real_int_real_real_2604;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, double, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2605;
typedef std::tuple<
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    std::vector<int>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2606;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, std::vector<double>, double, empty>
    type_vv_real_real_int_real_real_2607;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, std::vector<double>, std::vector<double>,
                   empty>
    type_vv_real_real_int_real_real_2608;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2609;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_2610;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, std::vector<double>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_2611;
typedef std::tuple<
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    std::vector<int>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2612;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   double, empty>
    type_vv_real_real_int_real_real_2613;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2614;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2615;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   var, empty>
    type_vv_real_real_int_real_real_2616;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_2617;
typedef std::tuple<
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2618;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, var, double, empty>
    type_vv_real_real_int_real_real_2619;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, var, std::vector<double>, empty>
    type_vv_real_real_int_real_real_2620;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2621;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, var, var, empty>
    type_vv_real_real_int_real_real_2622;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, var, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2623;
typedef std::tuple<
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    std::vector<int>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2624;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_2625;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, std::vector<var>, std::vector<double>,
                   empty>
    type_vv_real_real_int_real_real_2626;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2627;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_2628;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<int>, std::vector<var>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2629;
typedef std::tuple<
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    std::vector<int>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2630;
typedef std::tuple<
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_2631;
typedef std::tuple<
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_2632;
typedef std::tuple<
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2633;
typedef std::tuple<
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_2634;
typedef std::tuple<
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_2635;
typedef std::tuple<
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2636;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double, double, empty>
    type_vv_real_real_int_real_real_2637;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2638;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2639;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double, var, empty>
    type_vv_real_real_int_real_real_2640;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_2641;
typedef std::tuple<
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2642;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   double, empty>
    type_vv_real_real_int_real_real_2643;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2644;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2645;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   var, empty>
    type_vv_real_real_int_real_real_2646;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_2647;
typedef std::tuple<
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2648;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_2649;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2650;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2651;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_2652;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_2653;
typedef std::tuple<
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2654;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var, double, empty>
    type_vv_real_real_int_real_real_2655;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2656;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2657;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var, var, empty>
    type_vv_real_real_int_real_real_2658;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_2659;
typedef std::tuple<
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2660;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   double, empty>
    type_vv_real_real_int_real_real_2661;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2662;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2663;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>, var,
                   empty>
    type_vv_real_real_int_real_real_2664;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_2665;
typedef std::tuple<
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2666;
typedef std::tuple<
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_2667;
typedef std::tuple<
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_2668;
typedef std::tuple<
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2669;
typedef std::tuple<
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_2670;
typedef std::tuple<
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_2671;
typedef std::tuple<
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2672;
typedef std::tuple<std::vector<var>, var, int, double, double, empty>
    type_vv_real_real_int_real_real_2673;
typedef std::tuple<std::vector<var>, var, int, double, std::vector<double>,
                   empty>
    type_vv_real_real_int_real_real_2674;
typedef std::tuple<std::vector<var>, var, int, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2675;
typedef std::tuple<std::vector<var>, var, int, double, var, empty>
    type_vv_real_real_int_real_real_2676;
typedef std::tuple<std::vector<var>, var, int, double, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2677;
typedef std::tuple<
    std::vector<var>, var, int, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2678;
typedef std::tuple<std::vector<var>, var, int, std::vector<double>, double,
                   empty>
    type_vv_real_real_int_real_real_2679;
typedef std::tuple<std::vector<var>, var, int, std::vector<double>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2680;
typedef std::tuple<std::vector<var>, var, int, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2681;
typedef std::tuple<std::vector<var>, var, int, std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_2682;
typedef std::tuple<std::vector<var>, var, int, std::vector<double>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_2683;
typedef std::tuple<
    std::vector<var>, var, int, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2684;
typedef std::tuple<std::vector<var>, var, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_2685;
typedef std::tuple<std::vector<var>, var, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2686;
typedef std::tuple<std::vector<var>, var, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2687;
typedef std::tuple<std::vector<var>, var, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_2688;
typedef std::tuple<std::vector<var>, var, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_2689;
typedef std::tuple<
    std::vector<var>, var, int, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2690;
typedef std::tuple<std::vector<var>, var, int, var, double, empty>
    type_vv_real_real_int_real_real_2691;
typedef std::tuple<std::vector<var>, var, int, var, std::vector<double>, empty>
    type_vv_real_real_int_real_real_2692;
typedef std::tuple<std::vector<var>, var, int, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2693;
typedef std::tuple<std::vector<var>, var, int, var, var, empty>
    type_vv_real_real_int_real_real_2694;
typedef std::tuple<std::vector<var>, var, int, var, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2695;
typedef std::tuple<
    std::vector<var>, var, int, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2696;
typedef std::tuple<std::vector<var>, var, int, std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_2697;
typedef std::tuple<std::vector<var>, var, int, std::vector<var>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2698;
typedef std::tuple<std::vector<var>, var, int, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2699;
typedef std::tuple<std::vector<var>, var, int, std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_2700;
typedef std::tuple<std::vector<var>, var, int, std::vector<var>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_2701;
typedef std::tuple<
    std::vector<var>, var, int, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2702;
typedef std::tuple<
    std::vector<var>, var, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_2703;
typedef std::tuple<
    std::vector<var>, var, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_2704;
typedef std::tuple<
    std::vector<var>, var, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2705;
typedef std::tuple<
    std::vector<var>, var, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_2706;
typedef std::tuple<
    std::vector<var>, var, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_2707;
typedef std::tuple<
    std::vector<var>, var, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2708;
typedef std::tuple<std::vector<var>, var, std::vector<int>, double, double,
                   empty>
    type_vv_real_real_int_real_real_2709;
typedef std::tuple<std::vector<var>, var, std::vector<int>, double,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2710;
typedef std::tuple<std::vector<var>, var, std::vector<int>, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2711;
typedef std::tuple<std::vector<var>, var, std::vector<int>, double, var, empty>
    type_vv_real_real_int_real_real_2712;
typedef std::tuple<std::vector<var>, var, std::vector<int>, double,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_2713;
typedef std::tuple<
    std::vector<var>, var, std::vector<int>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2714;
typedef std::tuple<std::vector<var>, var, std::vector<int>, std::vector<double>,
                   double, empty>
    type_vv_real_real_int_real_real_2715;
typedef std::tuple<std::vector<var>, var, std::vector<int>, std::vector<double>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2716;
typedef std::tuple<std::vector<var>, var, std::vector<int>, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2717;
typedef std::tuple<std::vector<var>, var, std::vector<int>, std::vector<double>,
                   var, empty>
    type_vv_real_real_int_real_real_2718;
typedef std::tuple<std::vector<var>, var, std::vector<int>, std::vector<double>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_2719;
typedef std::tuple<
    std::vector<var>, var, std::vector<int>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2720;
typedef std::tuple<std::vector<var>, var, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_2721;
typedef std::tuple<std::vector<var>, var, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2722;
typedef std::tuple<std::vector<var>, var, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2723;
typedef std::tuple<std::vector<var>, var, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_2724;
typedef std::tuple<std::vector<var>, var, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_2725;
typedef std::tuple<
    std::vector<var>, var, std::vector<int>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2726;
typedef std::tuple<std::vector<var>, var, std::vector<int>, var, double, empty>
    type_vv_real_real_int_real_real_2727;
typedef std::tuple<std::vector<var>, var, std::vector<int>, var,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2728;
typedef std::tuple<std::vector<var>, var, std::vector<int>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2729;
typedef std::tuple<std::vector<var>, var, std::vector<int>, var, var, empty>
    type_vv_real_real_int_real_real_2730;
typedef std::tuple<std::vector<var>, var, std::vector<int>, var,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_2731;
typedef std::tuple<
    std::vector<var>, var, std::vector<int>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2732;
typedef std::tuple<std::vector<var>, var, std::vector<int>, std::vector<var>,
                   double, empty>
    type_vv_real_real_int_real_real_2733;
typedef std::tuple<std::vector<var>, var, std::vector<int>, std::vector<var>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2734;
typedef std::tuple<std::vector<var>, var, std::vector<int>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2735;
typedef std::tuple<std::vector<var>, var, std::vector<int>, std::vector<var>,
                   var, empty>
    type_vv_real_real_int_real_real_2736;
typedef std::tuple<std::vector<var>, var, std::vector<int>, std::vector<var>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_2737;
typedef std::tuple<
    std::vector<var>, var, std::vector<int>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2738;
typedef std::tuple<
    std::vector<var>, var, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_2739;
typedef std::tuple<
    std::vector<var>, var, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_2740;
typedef std::tuple<
    std::vector<var>, var, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2741;
typedef std::tuple<
    std::vector<var>, var, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_2742;
typedef std::tuple<
    std::vector<var>, var, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_2743;
typedef std::tuple<
    std::vector<var>, var, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2744;
typedef std::tuple<std::vector<var>, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   double, double, empty>
    type_vv_real_real_int_real_real_2745;
typedef std::tuple<std::vector<var>, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   double, std::vector<double>, empty>
    type_vv_real_real_int_real_real_2746;
typedef std::tuple<std::vector<var>, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   double, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2747;
typedef std::tuple<std::vector<var>, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   double, var, empty>
    type_vv_real_real_int_real_real_2748;
typedef std::tuple<std::vector<var>, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   double, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2749;
typedef std::tuple<
    std::vector<var>, var, Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2750;
typedef std::tuple<std::vector<var>, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<double>, double, empty>
    type_vv_real_real_int_real_real_2751;
typedef std::tuple<std::vector<var>, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<double>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_2752;
typedef std::tuple<std::vector<var>, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2753;
typedef std::tuple<std::vector<var>, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_2754;
typedef std::tuple<std::vector<var>, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<double>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2755;
typedef std::tuple<
    std::vector<var>, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2756;
typedef std::tuple<std::vector<var>, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_2757;
typedef std::tuple<std::vector<var>, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2758;
typedef std::tuple<std::vector<var>, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2759;
typedef std::tuple<std::vector<var>, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_2760;
typedef std::tuple<std::vector<var>, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_2761;
typedef std::tuple<
    std::vector<var>, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2762;
typedef std::tuple<std::vector<var>, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   var, double, empty>
    type_vv_real_real_int_real_real_2763;
typedef std::tuple<std::vector<var>, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   var, std::vector<double>, empty>
    type_vv_real_real_int_real_real_2764;
typedef std::tuple<std::vector<var>, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   var, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2765;
typedef std::tuple<std::vector<var>, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   var, var, empty>
    type_vv_real_real_int_real_real_2766;
typedef std::tuple<std::vector<var>, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   var, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2767;
typedef std::tuple<
    std::vector<var>, var, Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2768;
typedef std::tuple<std::vector<var>, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_2769;
typedef std::tuple<std::vector<var>, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<var>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_2770;
typedef std::tuple<std::vector<var>, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   empty>
    type_vv_real_real_int_real_real_2771;
typedef std::tuple<std::vector<var>, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_2772;
typedef std::tuple<std::vector<var>, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   std::vector<var>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2773;
typedef std::tuple<
    std::vector<var>, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2774;
typedef std::tuple<
    std::vector<var>, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_2775;
typedef std::tuple<
    std::vector<var>, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_2776;
typedef std::tuple<
    std::vector<var>, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2777;
typedef std::tuple<
    std::vector<var>, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_2778;
typedef std::tuple<
    std::vector<var>, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_2779;
typedef std::tuple<
    std::vector<var>, var, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2780;
typedef std::tuple<std::vector<var>, std::vector<var>, int, double, double,
                   empty>
    type_vv_real_real_int_real_real_2781;
typedef std::tuple<std::vector<var>, std::vector<var>, int, double,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2782;
typedef std::tuple<std::vector<var>, std::vector<var>, int, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2783;
typedef std::tuple<std::vector<var>, std::vector<var>, int, double, var, empty>
    type_vv_real_real_int_real_real_2784;
typedef std::tuple<std::vector<var>, std::vector<var>, int, double,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_2785;
typedef std::tuple<
    std::vector<var>, std::vector<var>, int, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2786;
typedef std::tuple<std::vector<var>, std::vector<var>, int, std::vector<double>,
                   double, empty>
    type_vv_real_real_int_real_real_2787;
typedef std::tuple<std::vector<var>, std::vector<var>, int, std::vector<double>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2788;
typedef std::tuple<std::vector<var>, std::vector<var>, int, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2789;
typedef std::tuple<std::vector<var>, std::vector<var>, int, std::vector<double>,
                   var, empty>
    type_vv_real_real_int_real_real_2790;
typedef std::tuple<std::vector<var>, std::vector<var>, int, std::vector<double>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_2791;
typedef std::tuple<
    std::vector<var>, std::vector<var>, int, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2792;
typedef std::tuple<std::vector<var>, std::vector<var>, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_2793;
typedef std::tuple<std::vector<var>, std::vector<var>, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2794;
typedef std::tuple<std::vector<var>, std::vector<var>, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2795;
typedef std::tuple<std::vector<var>, std::vector<var>, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_2796;
typedef std::tuple<std::vector<var>, std::vector<var>, int,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_2797;
typedef std::tuple<
    std::vector<var>, std::vector<var>, int,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2798;
typedef std::tuple<std::vector<var>, std::vector<var>, int, var, double, empty>
    type_vv_real_real_int_real_real_2799;
typedef std::tuple<std::vector<var>, std::vector<var>, int, var,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2800;
typedef std::tuple<std::vector<var>, std::vector<var>, int, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2801;
typedef std::tuple<std::vector<var>, std::vector<var>, int, var, var, empty>
    type_vv_real_real_int_real_real_2802;
typedef std::tuple<std::vector<var>, std::vector<var>, int, var,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_2803;
typedef std::tuple<
    std::vector<var>, std::vector<var>, int, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2804;
typedef std::tuple<std::vector<var>, std::vector<var>, int, std::vector<var>,
                   double, empty>
    type_vv_real_real_int_real_real_2805;
typedef std::tuple<std::vector<var>, std::vector<var>, int, std::vector<var>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2806;
typedef std::tuple<std::vector<var>, std::vector<var>, int, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2807;
typedef std::tuple<std::vector<var>, std::vector<var>, int, std::vector<var>,
                   var, empty>
    type_vv_real_real_int_real_real_2808;
typedef std::tuple<std::vector<var>, std::vector<var>, int, std::vector<var>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_2809;
typedef std::tuple<
    std::vector<var>, std::vector<var>, int, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2810;
typedef std::tuple<
    std::vector<var>, std::vector<var>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_2811;
typedef std::tuple<
    std::vector<var>, std::vector<var>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_2812;
typedef std::tuple<
    std::vector<var>, std::vector<var>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2813;
typedef std::tuple<
    std::vector<var>, std::vector<var>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_2814;
typedef std::tuple<
    std::vector<var>, std::vector<var>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_2815;
typedef std::tuple<
    std::vector<var>, std::vector<var>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2816;
typedef std::tuple<std::vector<var>, std::vector<var>, std::vector<int>, double,
                   double, empty>
    type_vv_real_real_int_real_real_2817;
typedef std::tuple<std::vector<var>, std::vector<var>, std::vector<int>, double,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2818;
typedef std::tuple<std::vector<var>, std::vector<var>, std::vector<int>, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2819;
typedef std::tuple<std::vector<var>, std::vector<var>, std::vector<int>, double,
                   var, empty>
    type_vv_real_real_int_real_real_2820;
typedef std::tuple<std::vector<var>, std::vector<var>, std::vector<int>, double,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_2821;
typedef std::tuple<
    std::vector<var>, std::vector<var>, std::vector<int>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2822;
typedef std::tuple<std::vector<var>, std::vector<var>, std::vector<int>,
                   std::vector<double>, double, empty>
    type_vv_real_real_int_real_real_2823;
typedef std::tuple<std::vector<var>, std::vector<var>, std::vector<int>,
                   std::vector<double>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_2824;
typedef std::tuple<std::vector<var>, std::vector<var>, std::vector<int>,
                   std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2825;
typedef std::tuple<std::vector<var>, std::vector<var>, std::vector<int>,
                   std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_2826;
typedef std::tuple<std::vector<var>, std::vector<var>, std::vector<int>,
                   std::vector<double>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2827;
typedef std::tuple<
    std::vector<var>, std::vector<var>, std::vector<int>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2828;
typedef std::tuple<std::vector<var>, std::vector<var>, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_2829;
typedef std::tuple<std::vector<var>, std::vector<var>, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2830;
typedef std::tuple<std::vector<var>, std::vector<var>, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2831;
typedef std::tuple<std::vector<var>, std::vector<var>, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_2832;
typedef std::tuple<std::vector<var>, std::vector<var>, std::vector<int>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_2833;
typedef std::tuple<
    std::vector<var>, std::vector<var>, std::vector<int>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2834;
typedef std::tuple<std::vector<var>, std::vector<var>, std::vector<int>, var,
                   double, empty>
    type_vv_real_real_int_real_real_2835;
typedef std::tuple<std::vector<var>, std::vector<var>, std::vector<int>, var,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2836;
typedef std::tuple<std::vector<var>, std::vector<var>, std::vector<int>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2837;
typedef std::tuple<std::vector<var>, std::vector<var>, std::vector<int>, var,
                   var, empty>
    type_vv_real_real_int_real_real_2838;
typedef std::tuple<std::vector<var>, std::vector<var>, std::vector<int>, var,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_2839;
typedef std::tuple<
    std::vector<var>, std::vector<var>, std::vector<int>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2840;
typedef std::tuple<std::vector<var>, std::vector<var>, std::vector<int>,
                   std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_2841;
typedef std::tuple<std::vector<var>, std::vector<var>, std::vector<int>,
                   std::vector<var>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_2842;
typedef std::tuple<std::vector<var>, std::vector<var>, std::vector<int>,
                   std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   empty>
    type_vv_real_real_int_real_real_2843;
typedef std::tuple<std::vector<var>, std::vector<var>, std::vector<int>,
                   std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_2844;
typedef std::tuple<std::vector<var>, std::vector<var>, std::vector<int>,
                   std::vector<var>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2845;
typedef std::tuple<
    std::vector<var>, std::vector<var>, std::vector<int>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2846;
typedef std::tuple<
    std::vector<var>, std::vector<var>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_2847;
typedef std::tuple<
    std::vector<var>, std::vector<var>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_2848;
typedef std::tuple<
    std::vector<var>, std::vector<var>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2849;
typedef std::tuple<
    std::vector<var>, std::vector<var>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_2850;
typedef std::tuple<
    std::vector<var>, std::vector<var>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_2851;
typedef std::tuple<
    std::vector<var>, std::vector<var>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2852;
typedef std::tuple<std::vector<var>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double, double, empty>
    type_vv_real_real_int_real_real_2853;
typedef std::tuple<std::vector<var>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2854;
typedef std::tuple<std::vector<var>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2855;
typedef std::tuple<std::vector<var>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double, var, empty>
    type_vv_real_real_int_real_real_2856;
typedef std::tuple<std::vector<var>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_2857;
typedef std::tuple<
    std::vector<var>, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    empty>
    type_vv_real_real_int_real_real_2858;
typedef std::tuple<std::vector<var>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   double, empty>
    type_vv_real_real_int_real_real_2859;
typedef std::tuple<std::vector<var>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2860;
typedef std::tuple<std::vector<var>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2861;
typedef std::tuple<std::vector<var>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   var, empty>
    type_vv_real_real_int_real_real_2862;
typedef std::tuple<std::vector<var>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_2863;
typedef std::tuple<
    std::vector<var>, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2864;
typedef std::tuple<std::vector<var>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_2865;
typedef std::tuple<
    std::vector<var>, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_2866;
typedef std::tuple<std::vector<var>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2867;
typedef std::tuple<std::vector<var>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_2868;
typedef std::tuple<
    std::vector<var>, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2869;
typedef std::tuple<
    std::vector<var>, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2870;
typedef std::tuple<std::vector<var>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var, double, empty>
    type_vv_real_real_int_real_real_2871;
typedef std::tuple<std::vector<var>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2872;
typedef std::tuple<std::vector<var>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2873;
typedef std::tuple<std::vector<var>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var, var, empty>
    type_vv_real_real_int_real_real_2874;
typedef std::tuple<std::vector<var>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, var, std::vector<var>,
                   empty>
    type_vv_real_real_int_real_real_2875;
typedef std::tuple<
    std::vector<var>, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2876;
typedef std::tuple<std::vector<var>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   double, empty>
    type_vv_real_real_int_real_real_2877;
typedef std::tuple<std::vector<var>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<double>, empty>
    type_vv_real_real_int_real_real_2878;
typedef std::tuple<std::vector<var>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2879;
typedef std::tuple<std::vector<var>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>, var,
                   empty>
    type_vv_real_real_int_real_real_2880;
typedef std::tuple<std::vector<var>, std::vector<var>,
                   Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<var>, empty>
    type_vv_real_real_int_real_real_2881;
typedef std::tuple<
    std::vector<var>, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2882;
typedef std::tuple<
    std::vector<var>, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_2883;
typedef std::tuple<
    std::vector<var>, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_2884;
typedef std::tuple<
    std::vector<var>, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2885;
typedef std::tuple<
    std::vector<var>, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_2886;
typedef std::tuple<
    std::vector<var>, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_2887;
typedef std::tuple<
    std::vector<var>, std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2888;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    double, double, empty>
    type_vv_real_real_int_real_real_2889;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    double, std::vector<double>, empty>
    type_vv_real_real_int_real_real_2890;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    double, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2891;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    double, var, empty>
    type_vv_real_real_int_real_real_2892;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    double, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2893;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    empty>
    type_vv_real_real_int_real_real_2894;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<double>, double, empty>
    type_vv_real_real_int_real_real_2895;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<double>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_2896;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2897;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_2898;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<double>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2899;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2900;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_2901;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_2902;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2903;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_2904;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2905;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2906;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int, var,
    double, empty>
    type_vv_real_real_int_real_real_2907;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int, var,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_2908;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int, var,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2909;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int, var,
    var, empty>
    type_vv_real_real_int_real_real_2910;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int, var,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_2911;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2912;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_2913;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<var>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_2914;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2915;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_2916;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<var>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2917;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2918;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_2919;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_2920;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2921;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_2922;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_2923;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2924;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, double, double, empty>
    type_vv_real_real_int_real_real_2925;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, double, std::vector<double>, empty>
    type_vv_real_real_int_real_real_2926;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, double, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2927;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, double, var, empty>
    type_vv_real_real_int_real_real_2928;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, double, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2929;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2930;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<double>, double, empty>
    type_vv_real_real_int_real_real_2931;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<double>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_2932;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<double>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2933;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_2934;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<double>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2935;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2936;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_2937;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_2938;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2939;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_2940;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_2941;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2942;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, var, double, empty>
    type_vv_real_real_int_real_real_2943;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, var, std::vector<double>, empty>
    type_vv_real_real_int_real_real_2944;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, var, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2945;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, var, var, empty>
    type_vv_real_real_int_real_real_2946;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, var, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2947;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2948;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_2949;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<var>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_2950;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<var>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2951;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_2952;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<var>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2953;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2954;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_2955;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_2956;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2957;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_2958;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_2959;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2960;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double, double, empty>
    type_vv_real_real_int_real_real_2961;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double, std::vector<double>, empty>
    type_vv_real_real_int_real_real_2962;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2963;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double, var, empty>
    type_vv_real_real_int_real_real_2964;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2965;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2966;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>, double, empty>
    type_vv_real_real_int_real_real_2967;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_2968;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2969;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_2970;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_2971;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2972;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_2973;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_2974;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2975;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_2976;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2977;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2978;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var, double, empty>
    type_vv_real_real_int_real_real_2979;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var, std::vector<double>, empty>
    type_vv_real_real_int_real_real_2980;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2981;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var, var, empty>
    type_vv_real_real_int_real_real_2982;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var, std::vector<var>, empty>
    type_vv_real_real_int_real_real_2983;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2984;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_2985;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_2986;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2987;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_2988;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>, std::vector<var>,
    empty>
    type_vv_real_real_int_real_real_2989;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2990;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_2991;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_2992;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2993;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_2994;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_2995;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_2996;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    int, double, double, empty>
    type_vv_real_real_int_real_real_2997;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    int, double, std::vector<double>, empty>
    type_vv_real_real_int_real_real_2998;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    int, double, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_2999;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    int, double, var, empty>
    type_vv_real_real_int_real_real_3000;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    int, double, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3001;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    int, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3002;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    int, std::vector<double>, double, empty>
    type_vv_real_real_int_real_real_3003;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    int, std::vector<double>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3004;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    int, std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3005;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    int, std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_3006;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    int, std::vector<double>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3007;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    int, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3008;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    int, Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_3009;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    int, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3010;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    int, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3011;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    int, Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_3012;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    int, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3013;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    int, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3014;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    int, var, double, empty>
    type_vv_real_real_int_real_real_3015;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    int, var, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3016;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    int, var, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3017;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    int, var, var, empty>
    type_vv_real_real_int_real_real_3018;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    int, var, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3019;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    int, var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    empty>
    type_vv_real_real_int_real_real_3020;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    int, std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_3021;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    int, std::vector<var>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3022;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    int, std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3023;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    int, std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_3024;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    int, std::vector<var>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3025;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    int, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3026;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    int, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    double, empty>
    type_vv_real_real_int_real_real_3027;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    int, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_3028;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    int, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3029;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    int, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    empty>
    type_vv_real_real_int_real_real_3030;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    int, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_3031;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    int, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3032;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    std::vector<int>, double, double, empty>
    type_vv_real_real_int_real_real_3033;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    std::vector<int>, double, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3034;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    std::vector<int>, double, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3035;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    std::vector<int>, double, var, empty>
    type_vv_real_real_int_real_real_3036;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    std::vector<int>, double, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3037;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    std::vector<int>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3038;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    std::vector<int>, std::vector<double>, double, empty>
    type_vv_real_real_int_real_real_3039;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    std::vector<int>, std::vector<double>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3040;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    std::vector<int>, std::vector<double>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3041;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    std::vector<int>, std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_3042;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    std::vector<int>, std::vector<double>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3043;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    std::vector<int>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3044;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_3045;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_3046;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3047;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_3048;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_3049;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3050;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    std::vector<int>, var, double, empty>
    type_vv_real_real_int_real_real_3051;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    std::vector<int>, var, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3052;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    std::vector<int>, var, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3053;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    std::vector<int>, var, var, empty>
    type_vv_real_real_int_real_real_3054;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    std::vector<int>, var, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3055;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    std::vector<int>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3056;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    std::vector<int>, std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_3057;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    std::vector<int>, std::vector<var>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3058;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    std::vector<int>, std::vector<var>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3059;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    std::vector<int>, std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_3060;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    std::vector<int>, std::vector<var>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3061;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    std::vector<int>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3062;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_3063;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_3064;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3065;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_3066;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_3067;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3068;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double, double, empty>
    type_vv_real_real_int_real_real_3069;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3070;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3071;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double, var, empty>
    type_vv_real_real_int_real_real_3072;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3073;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3074;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>, double, empty>
    type_vv_real_real_int_real_real_3075;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_3076;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3077;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_3078;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_3079;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3080;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_3081;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3082;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3083;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_3084;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3085;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3086;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var, double, empty>
    type_vv_real_real_int_real_real_3087;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3088;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3089;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var, var, empty>
    type_vv_real_real_int_real_real_3090;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3091;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3092;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_3093;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_3094;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3095;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_3096;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>, std::vector<var>,
    empty>
    type_vv_real_real_int_real_real_3097;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3098;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_3099;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_3100;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3101;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_3102;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_3103;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3104;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, int, double, double, empty>
    type_vv_real_real_int_real_real_3105;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, int, double, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3106;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, int, double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    empty>
    type_vv_real_real_int_real_real_3107;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, int, double, var, empty>
    type_vv_real_real_int_real_real_3108;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, int, double, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3109;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, int, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3110;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, int, std::vector<double>, double, empty>
    type_vv_real_real_int_real_real_3111;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, int, std::vector<double>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3112;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, int, std::vector<double>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3113;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, int, std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_3114;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, int, std::vector<double>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3115;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, int, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3116;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, int, Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
    empty>
    type_vv_real_real_int_real_real_3117;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, int, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_3118;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, int, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3119;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, int, Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
    empty>
    type_vv_real_real_int_real_real_3120;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, int, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_3121;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, int, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3122;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, int, var, double, empty>
    type_vv_real_real_int_real_real_3123;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, int, var, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3124;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, int, var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    empty>
    type_vv_real_real_int_real_real_3125;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, int, var, var, empty>
    type_vv_real_real_int_real_real_3126;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, int, var, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3127;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, int, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3128;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, int, std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_3129;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, int, std::vector<var>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3130;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, int, std::vector<var>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3131;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, int, std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_3132;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, int, std::vector<var>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3133;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, int, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3134;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_3135;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_3136;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3137;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_3138;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_3139;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3140;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<int>, double, double, empty>
    type_vv_real_real_int_real_real_3141;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<int>, double, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3142;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<int>, double,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3143;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<int>, double, var, empty>
    type_vv_real_real_int_real_real_3144;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<int>, double, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3145;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<int>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3146;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<int>, std::vector<double>, double, empty>
    type_vv_real_real_int_real_real_3147;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<int>, std::vector<double>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_3148;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<int>, std::vector<double>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3149;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<int>, std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_3150;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<int>, std::vector<double>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_3151;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<int>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3152;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<int>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_3153;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<int>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3154;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<int>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3155;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<int>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_3156;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<int>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3157;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<int>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3158;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<int>, var, double, empty>
    type_vv_real_real_int_real_real_3159;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<int>, var, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3160;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<int>, var,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3161;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<int>, var, var, empty>
    type_vv_real_real_int_real_real_3162;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<int>, var, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3163;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<int>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3164;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<int>, std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_3165;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<int>, std::vector<var>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_3166;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<int>, std::vector<var>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3167;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<int>, std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_3168;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<int>, std::vector<var>, std::vector<var>,
    empty>
    type_vv_real_real_int_real_real_3169;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<int>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3170;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_3171;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_3172;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3173;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_3174;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_3175;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3176;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>, double, double,
    empty>
    type_vv_real_real_int_real_real_3177;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_3178;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3179;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>, double, var,
    empty>
    type_vv_real_real_int_real_real_3180;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_3181;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3182;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    std::vector<double>, double, empty>
    type_vv_real_real_int_real_real_3183;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    std::vector<double>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3184;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3185;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_3186;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    std::vector<double>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3187;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3188;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_3189;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3190;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3191;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_3192;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3193;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3194;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>, var, double,
    empty>
    type_vv_real_real_int_real_real_3195;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_3196;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3197;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>, var, var, empty>
    type_vv_real_real_int_real_real_3198;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_3199;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3200;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_3201;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    std::vector<var>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3202;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3203;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_3204;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    std::vector<var>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3205;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3206;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_3207;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_3208;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3209;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_3210;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_3211;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3212;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, int, double, double, empty>
    type_vv_real_real_int_real_real_3213;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, int, double, std::vector<double>,
    empty>
    type_vv_real_real_int_real_real_3214;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, int, double,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3215;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, int, double, var, empty>
    type_vv_real_real_int_real_real_3216;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, int, double, std::vector<var>,
    empty>
    type_vv_real_real_int_real_real_3217;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, int, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3218;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, int, std::vector<double>, double,
    empty>
    type_vv_real_real_int_real_real_3219;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, int, std::vector<double>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_3220;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, int, std::vector<double>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3221;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, int, std::vector<double>, var,
    empty>
    type_vv_real_real_int_real_real_3222;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, int, std::vector<double>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_3223;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, int, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3224;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_3225;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3226;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3227;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_3228;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3229;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3230;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, int, var, double, empty>
    type_vv_real_real_int_real_real_3231;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, int, var, std::vector<double>,
    empty>
    type_vv_real_real_int_real_real_3232;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, int, var,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3233;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, int, var, var, empty>
    type_vv_real_real_int_real_real_3234;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, int, var, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3235;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, int, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3236;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, int, std::vector<var>, double,
    empty>
    type_vv_real_real_int_real_real_3237;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, int, std::vector<var>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_3238;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, int, std::vector<var>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3239;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, int, std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_3240;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, int, std::vector<var>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_3241;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, int, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3242;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_3243;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_3244;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3245;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_3246;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_3247;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3248;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>, double, double,
    empty>
    type_vv_real_real_int_real_real_3249;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>, double,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_3250;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>, double,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3251;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>, double, var,
    empty>
    type_vv_real_real_int_real_real_3252;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>, double,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_3253;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3254;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
    std::vector<double>, double, empty>
    type_vv_real_real_int_real_real_3255;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
    std::vector<double>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3256;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3257;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
    std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_3258;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
    std::vector<double>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3259;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3260;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_3261;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3262;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3263;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_3264;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3265;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3266;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>, var, double,
    empty>
    type_vv_real_real_int_real_real_3267;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>, var,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_3268;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>, var,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3269;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>, var, var, empty>
    type_vv_real_real_int_real_real_3270;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>, var,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_3271;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3272;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
    std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_3273;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
    std::vector<var>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3274;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3275;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
    std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_3276;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
    std::vector<var>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3277;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3278;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_3279;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_3280;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3281;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_3282;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_3283;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3284;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double, double, empty>
    type_vv_real_real_int_real_real_3285;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3286;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3287;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double, var, empty>
    type_vv_real_real_int_real_real_3288;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3289;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3290;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>, double, empty>
    type_vv_real_real_int_real_real_3291;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_3292;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3293;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_3294;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_3295;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3296;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_3297;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3298;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3299;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_3300;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3301;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3302;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var, double, empty>
    type_vv_real_real_int_real_real_3303;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3304;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3305;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var, var, empty>
    type_vv_real_real_int_real_real_3306;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3307;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3308;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_3309;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_3310;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3311;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_3312;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>, std::vector<var>,
    empty>
    type_vv_real_real_int_real_real_3313;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3314;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_3315;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_3316;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3317;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_3318;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_3319;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3320;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, int,
    double, double, empty>
    type_vv_real_real_int_real_real_3321;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, int,
    double, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3322;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, int,
    double, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3323;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, int,
    double, var, empty>
    type_vv_real_real_int_real_real_3324;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, int,
    double, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3325;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, int,
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    empty>
    type_vv_real_real_int_real_real_3326;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, int,
    std::vector<double>, double, empty>
    type_vv_real_real_int_real_real_3327;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, int,
    std::vector<double>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3328;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, int,
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3329;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, int,
    std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_3330;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, int,
    std::vector<double>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3331;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, int,
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3332;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, int,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_3333;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, int,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3334;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, int,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3335;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, int,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_3336;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, int,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3337;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, int,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3338;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, int,
    var, double, empty>
    type_vv_real_real_int_real_real_3339;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, int,
    var, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3340;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, int,
    var, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3341;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, int,
    var, var, empty>
    type_vv_real_real_int_real_real_3342;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, int,
    var, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3343;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, int,
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3344;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, int,
    std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_3345;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, int,
    std::vector<var>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3346;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, int,
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3347;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, int,
    std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_3348;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, int,
    std::vector<var>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3349;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, int,
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3350;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_3351;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_3352;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3353;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_3354;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_3355;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3356;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<int>, double, double, empty>
    type_vv_real_real_int_real_real_3357;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<int>, double, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3358;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<int>, double, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3359;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<int>, double, var, empty>
    type_vv_real_real_int_real_real_3360;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<int>, double, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3361;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<int>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3362;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<int>, std::vector<double>, double, empty>
    type_vv_real_real_int_real_real_3363;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<int>, std::vector<double>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3364;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<int>, std::vector<double>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3365;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<int>, std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_3366;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<int>, std::vector<double>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3367;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<int>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3368;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_3369;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_3370;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3371;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_3372;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_3373;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3374;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<int>, var, double, empty>
    type_vv_real_real_int_real_real_3375;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<int>, var, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3376;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<int>, var, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3377;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<int>, var, var, empty>
    type_vv_real_real_int_real_real_3378;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<int>, var, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3379;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<int>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3380;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<int>, std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_3381;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<int>, std::vector<var>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3382;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<int>, std::vector<var>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3383;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<int>, std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_3384;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<int>, std::vector<var>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3385;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<int>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3386;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_3387;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_3388;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3389;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_3390;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_3391;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3392;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double, double, empty>
    type_vv_real_real_int_real_real_3393;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3394;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3395;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double, var, empty>
    type_vv_real_real_int_real_real_3396;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3397;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3398;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>, double, empty>
    type_vv_real_real_int_real_real_3399;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_3400;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3401;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_3402;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_3403;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3404;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_3405;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3406;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3407;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_3408;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3409;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3410;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var, double, empty>
    type_vv_real_real_int_real_real_3411;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3412;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3413;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var, var, empty>
    type_vv_real_real_int_real_real_3414;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3415;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3416;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_3417;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_3418;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3419;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_3420;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>, std::vector<var>,
    empty>
    type_vv_real_real_int_real_real_3421;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3422;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_3423;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_3424;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3425;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_3426;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_3427;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3428;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, int, double, double, empty>
    type_vv_real_real_int_real_real_3429;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, int, double, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3430;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, int, double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    empty>
    type_vv_real_real_int_real_real_3431;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, int, double, var, empty>
    type_vv_real_real_int_real_real_3432;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, int, double, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3433;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, int, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3434;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, int, std::vector<double>, double, empty>
    type_vv_real_real_int_real_real_3435;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, int, std::vector<double>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3436;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, int, std::vector<double>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3437;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, int, std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_3438;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, int, std::vector<double>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3439;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, int, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3440;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, int, Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
    empty>
    type_vv_real_real_int_real_real_3441;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, int, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_3442;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, int, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3443;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, int, Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_3444;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, int, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_3445;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, int, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3446;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, int, var, double, empty>
    type_vv_real_real_int_real_real_3447;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, int, var, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3448;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, int, var, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3449;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, int, var, var, empty>
    type_vv_real_real_int_real_real_3450;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, int, var, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3451;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, int, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3452;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, int, std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_3453;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, int, std::vector<var>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3454;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, int, std::vector<var>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3455;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, int, std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_3456;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, int, std::vector<var>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3457;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, int, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3458;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_3459;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_3460;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3461;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_3462;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_3463;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3464;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<int>, double, double, empty>
    type_vv_real_real_int_real_real_3465;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<int>, double, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3466;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<int>, double,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3467;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<int>, double, var, empty>
    type_vv_real_real_int_real_real_3468;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<int>, double, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3469;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<int>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3470;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<int>, std::vector<double>, double, empty>
    type_vv_real_real_int_real_real_3471;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<int>, std::vector<double>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_3472;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<int>, std::vector<double>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3473;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<int>, std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_3474;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<int>, std::vector<double>, std::vector<var>,
    empty>
    type_vv_real_real_int_real_real_3475;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<int>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3476;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<int>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_3477;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<int>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3478;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<int>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3479;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<int>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_3480;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<int>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3481;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<int>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3482;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<int>, var, double, empty>
    type_vv_real_real_int_real_real_3483;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<int>, var, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3484;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<int>, var,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3485;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<int>, var, var, empty>
    type_vv_real_real_int_real_real_3486;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<int>, var, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3487;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<int>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3488;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<int>, std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_3489;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<int>, std::vector<var>, std::vector<double>,
    empty>
    type_vv_real_real_int_real_real_3490;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<int>, std::vector<var>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3491;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<int>, std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_3492;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<int>, std::vector<var>, std::vector<var>,
    empty>
    type_vv_real_real_int_real_real_3493;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<int>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3494;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_3495;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_3496;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3497;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_3498;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_3499;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3500;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>, double, double,
    empty>
    type_vv_real_real_int_real_real_3501;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_3502;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3503;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>, double, var, empty>
    type_vv_real_real_int_real_real_3504;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_3505;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3506;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    std::vector<double>, double, empty>
    type_vv_real_real_int_real_real_3507;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    std::vector<double>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3508;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3509;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_3510;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    std::vector<double>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3511;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3512;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_3513;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3514;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3515;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_3516;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3517;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3518;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>, var, double, empty>
    type_vv_real_real_int_real_real_3519;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_3520;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3521;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>, var, var, empty>
    type_vv_real_real_int_real_real_3522;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_3523;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3524;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
    double, empty>
    type_vv_real_real_int_real_real_3525;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_3526;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3527;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
    var, empty>
    type_vv_real_real_int_real_real_3528;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_3529;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3530;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_3531;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_3532;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3533;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_3534;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_3535;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3536;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    double, double, empty>
    type_vv_real_real_int_real_real_3537;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    double, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3538;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    double, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3539;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    double, var, empty>
    type_vv_real_real_int_real_real_3540;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    double, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3541;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    empty>
    type_vv_real_real_int_real_real_3542;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<double>, double, empty>
    type_vv_real_real_int_real_real_3543;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<double>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3544;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3545;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_3546;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<double>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3547;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3548;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_3549;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3550;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3551;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_3552;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3553;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3554;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int, var,
    double, empty>
    type_vv_real_real_int_real_real_3555;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int, var,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_3556;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int, var,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3557;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int, var,
    var, empty>
    type_vv_real_real_int_real_real_3558;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int, var,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_3559;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3560;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_3561;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<var>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3562;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3563;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_3564;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<var>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3565;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3566;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_3567;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_3568;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3569;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_3570;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_3571;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, int,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3572;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, double, double, empty>
    type_vv_real_real_int_real_real_3573;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, double, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3574;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, double, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3575;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, double, var, empty>
    type_vv_real_real_int_real_real_3576;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, double, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3577;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3578;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<double>, double, empty>
    type_vv_real_real_int_real_real_3579;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<double>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3580;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<double>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3581;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_3582;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<double>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3583;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3584;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_3585;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_3586;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3587;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_3588;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_3589;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3590;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, var, double, empty>
    type_vv_real_real_int_real_real_3591;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, var, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3592;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, var, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3593;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, var, var, empty>
    type_vv_real_real_int_real_real_3594;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, var, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3595;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3596;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_3597;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<var>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3598;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<var>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3599;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_3600;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<var>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3601;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3602;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_3603;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_3604;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3605;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_3606;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_3607;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<int>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3608;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double, double, empty>
    type_vv_real_real_int_real_real_3609;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3610;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3611;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double, var, empty>
    type_vv_real_real_int_real_real_3612;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3613;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3614;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>, double, empty>
    type_vv_real_real_int_real_real_3615;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_3616;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3617;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>, var, empty>
    type_vv_real_real_int_real_real_3618;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_3619;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3620;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty>
    type_vv_real_real_int_real_real_3621;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3622;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3623;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty>
    type_vv_real_real_int_real_real_3624;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3625;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3626;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var, double, empty>
    type_vv_real_real_int_real_real_3627;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var, std::vector<double>, empty>
    type_vv_real_real_int_real_real_3628;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3629;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var, var, empty>
    type_vv_real_real_int_real_real_3630;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var, std::vector<var>, empty>
    type_vv_real_real_int_real_real_3631;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3632;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>, double, empty>
    type_vv_real_real_int_real_real_3633;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_3634;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3635;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>, var, empty>
    type_vv_real_real_int_real_real_3636;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>, std::vector<var>,
    empty>
    type_vv_real_real_int_real_real_3637;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3638;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty>
    type_vv_real_real_int_real_real_3639;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty>
    type_vv_real_real_int_real_real_3640;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty>
    type_vv_real_real_int_real_real_3641;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty>
    type_vv_real_real_int_real_real_3642;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty>
    type_vv_real_real_int_real_real_3643;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<int, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty>
    type_vv_real_real_int_real_real_3644;
