#include <gtest/gtest.h>
#include <tuple>
#include <test/prob/test_fixture_distr.hpp>
#include <test/prob/test_fixture_cdf.hpp>
#include <test/prob/test_fixture_cdf_log.hpp>
#include <test/prob/test_fixture_ccdf_log.hpp>

typedef std::tuple<double, double, double, var, empty, empty>
    type_vv_real_real_real_real_0;
typedef std::tuple<double, double, double, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_1;
typedef std::tuple<
    double, double, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_2;
typedef std::tuple<double, double, std::vector<double>, var, empty, empty>
    type_vv_real_real_real_real_3;
typedef std::tuple<double, double, std::vector<double>, std::vector<var>, empty,
                   empty>
    type_vv_real_real_real_real_4;
typedef std::tuple<
    double, double, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_5;
typedef std::tuple<double, double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   var, empty, empty>
    type_vv_real_real_real_real_6;
typedef std::tuple<double, double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<var>, empty, empty>
    type_vv_real_real_real_real_7;
typedef std::tuple<
    double, double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_8;
typedef std::tuple<double, double, var, double, empty, empty>
    type_vv_real_real_real_real_9;
typedef std::tuple<double, double, var, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_10;
typedef std::tuple<double, double, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_11;
typedef std::tuple<double, double, var, var, empty, empty>
    type_vv_real_real_real_real_12;
typedef std::tuple<double, double, var, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_13;
typedef std::tuple<
    double, double, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_14;
typedef std::tuple<double, double, std::vector<var>, double, empty, empty>
    type_vv_real_real_real_real_15;
typedef std::tuple<double, double, std::vector<var>, std::vector<double>, empty,
                   empty>
    type_vv_real_real_real_real_16;
typedef std::tuple<double, double, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_17;
typedef std::tuple<double, double, std::vector<var>, var, empty, empty>
    type_vv_real_real_real_real_18;
typedef std::tuple<double, double, std::vector<var>, std::vector<var>, empty,
                   empty>
    type_vv_real_real_real_real_19;
typedef std::tuple<
    double, double, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_20;
typedef std::tuple<
    double, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty, empty>
    type_vv_real_real_real_real_21;
typedef std::tuple<
    double, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty, empty>
    type_vv_real_real_real_real_22;
typedef std::tuple<
    double, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_23;
typedef std::tuple<
    double, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty,
    empty>
    type_vv_real_real_real_real_24;
typedef std::tuple<
    double, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty, empty>
    type_vv_real_real_real_real_25;
typedef std::tuple<
    double, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_26;
typedef std::tuple<double, std::vector<double>, double, var, empty, empty>
    type_vv_real_real_real_real_27;
typedef std::tuple<double, std::vector<double>, double, std::vector<var>, empty,
                   empty>
    type_vv_real_real_real_real_28;
typedef std::tuple<
    double, std::vector<double>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_29;
typedef std::tuple<double, std::vector<double>, std::vector<double>, var, empty,
                   empty>
    type_vv_real_real_real_real_30;
typedef std::tuple<double, std::vector<double>, std::vector<double>,
                   std::vector<var>, empty, empty>
    type_vv_real_real_real_real_31;
typedef std::tuple<
    double, std::vector<double>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_32;
typedef std::tuple<double, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty, empty>
    type_vv_real_real_real_real_33;
typedef std::tuple<double, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty, empty>
    type_vv_real_real_real_real_34;
typedef std::tuple<
    double, std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_35;
typedef std::tuple<double, std::vector<double>, var, double, empty, empty>
    type_vv_real_real_real_real_36;
typedef std::tuple<double, std::vector<double>, var, std::vector<double>, empty,
                   empty>
    type_vv_real_real_real_real_37;
typedef std::tuple<double, std::vector<double>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_38;
typedef std::tuple<double, std::vector<double>, var, var, empty, empty>
    type_vv_real_real_real_real_39;
typedef std::tuple<double, std::vector<double>, var, std::vector<var>, empty,
                   empty>
    type_vv_real_real_real_real_40;
typedef std::tuple<
    double, std::vector<double>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_41;
typedef std::tuple<double, std::vector<double>, std::vector<var>, double, empty,
                   empty>
    type_vv_real_real_real_real_42;
typedef std::tuple<double, std::vector<double>, std::vector<var>,
                   std::vector<double>, empty, empty>
    type_vv_real_real_real_real_43;
typedef std::tuple<double, std::vector<double>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_44;
typedef std::tuple<double, std::vector<double>, std::vector<var>, var, empty,
                   empty>
    type_vv_real_real_real_real_45;
typedef std::tuple<double, std::vector<double>, std::vector<var>,
                   std::vector<var>, empty, empty>
    type_vv_real_real_real_real_46;
typedef std::tuple<
    double, std::vector<double>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_47;
typedef std::tuple<
    double, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty, empty>
    type_vv_real_real_real_real_48;
typedef std::tuple<
    double, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty, empty>
    type_vv_real_real_real_real_49;
typedef std::tuple<
    double, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_50;
typedef std::tuple<
    double, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty,
    empty>
    type_vv_real_real_real_real_51;
typedef std::tuple<
    double, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty, empty>
    type_vv_real_real_real_real_52;
typedef std::tuple<
    double, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_53;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
                   var, empty, empty>
    type_vv_real_real_real_real_54;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
                   std::vector<var>, empty, empty>
    type_vv_real_real_real_real_55;
typedef std::tuple<
    double, Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_56;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, var, empty, empty>
    type_vv_real_real_real_real_57;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_58;
typedef std::tuple<
    double, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_59;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty, empty>
    type_vv_real_real_real_real_60;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty, empty>
    type_vv_real_real_real_real_61;
typedef std::tuple<
    double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_62;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   double, empty, empty>
    type_vv_real_real_real_real_63;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   std::vector<double>, empty, empty>
    type_vv_real_real_real_real_64;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_65;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>, var, var,
                   empty, empty>
    type_vv_real_real_real_real_66;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   std::vector<var>, empty, empty>
    type_vv_real_real_real_real_67;
typedef std::tuple<
    double, Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_68;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<var>, double, empty, empty>
    type_vv_real_real_real_real_69;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<var>, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_70;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   empty, empty>
    type_vv_real_real_real_real_71;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<var>, var, empty, empty>
    type_vv_real_real_real_real_72;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<var>, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_73;
typedef std::tuple<
    double, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_74;
typedef std::tuple<
    double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty, empty>
    type_vv_real_real_real_real_75;
typedef std::tuple<
    double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty, empty>
    type_vv_real_real_real_real_76;
typedef std::tuple<
    double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_77;
typedef std::tuple<
    double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty,
    empty>
    type_vv_real_real_real_real_78;
typedef std::tuple<
    double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty, empty>
    type_vv_real_real_real_real_79;
typedef std::tuple<
    double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_80;
typedef std::tuple<double, var, double, double, empty, empty>
    type_vv_real_real_real_real_81;
typedef std::tuple<double, var, double, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_82;
typedef std::tuple<double, var, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_83;
typedef std::tuple<double, var, double, var, empty, empty>
    type_vv_real_real_real_real_84;
typedef std::tuple<double, var, double, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_85;
typedef std::tuple<
    double, var, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_86;
typedef std::tuple<double, var, std::vector<double>, double, empty, empty>
    type_vv_real_real_real_real_87;
typedef std::tuple<double, var, std::vector<double>, std::vector<double>, empty,
                   empty>
    type_vv_real_real_real_real_88;
typedef std::tuple<double, var, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_89;
typedef std::tuple<double, var, std::vector<double>, var, empty, empty>
    type_vv_real_real_real_real_90;
typedef std::tuple<double, var, std::vector<double>, std::vector<var>, empty,
                   empty>
    type_vv_real_real_real_real_91;
typedef std::tuple<
    double, var, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_92;
typedef std::tuple<double, var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   double, empty, empty>
    type_vv_real_real_real_real_93;
typedef std::tuple<double, var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty, empty>
    type_vv_real_real_real_real_94;
typedef std::tuple<double, var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_95;
typedef std::tuple<double, var, Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   empty, empty>
    type_vv_real_real_real_real_96;
typedef std::tuple<double, var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<var>, empty, empty>
    type_vv_real_real_real_real_97;
typedef std::tuple<
    double, var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_98;
typedef std::tuple<double, var, var, double, empty, empty>
    type_vv_real_real_real_real_99;
typedef std::tuple<double, var, var, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_100;
typedef std::tuple<double, var, var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   empty, empty>
    type_vv_real_real_real_real_101;
typedef std::tuple<double, var, var, var, empty, empty>
    type_vv_real_real_real_real_102;
typedef std::tuple<double, var, var, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_103;
typedef std::tuple<
    double, var, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_104;
typedef std::tuple<double, var, std::vector<var>, double, empty, empty>
    type_vv_real_real_real_real_105;
typedef std::tuple<double, var, std::vector<var>, std::vector<double>, empty,
                   empty>
    type_vv_real_real_real_real_106;
typedef std::tuple<double, var, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_107;
typedef std::tuple<double, var, std::vector<var>, var, empty, empty>
    type_vv_real_real_real_real_108;
typedef std::tuple<double, var, std::vector<var>, std::vector<var>, empty,
                   empty>
    type_vv_real_real_real_real_109;
typedef std::tuple<
    double, var, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_110;
typedef std::tuple<
    double, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty, empty>
    type_vv_real_real_real_real_111;
typedef std::tuple<
    double, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty, empty>
    type_vv_real_real_real_real_112;
typedef std::tuple<
    double, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_113;
typedef std::tuple<
    double, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty,
    empty>
    type_vv_real_real_real_real_114;
typedef std::tuple<
    double, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty, empty>
    type_vv_real_real_real_real_115;
typedef std::tuple<
    double, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_116;
typedef std::tuple<double, std::vector<var>, double, double, empty, empty>
    type_vv_real_real_real_real_117;
typedef std::tuple<double, std::vector<var>, double, std::vector<double>, empty,
                   empty>
    type_vv_real_real_real_real_118;
typedef std::tuple<double, std::vector<var>, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_119;
typedef std::tuple<double, std::vector<var>, double, var, empty, empty>
    type_vv_real_real_real_real_120;
typedef std::tuple<double, std::vector<var>, double, std::vector<var>, empty,
                   empty>
    type_vv_real_real_real_real_121;
typedef std::tuple<
    double, std::vector<var>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_122;
typedef std::tuple<double, std::vector<var>, std::vector<double>, double, empty,
                   empty>
    type_vv_real_real_real_real_123;
typedef std::tuple<double, std::vector<var>, std::vector<double>,
                   std::vector<double>, empty, empty>
    type_vv_real_real_real_real_124;
typedef std::tuple<double, std::vector<var>, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_125;
typedef std::tuple<double, std::vector<var>, std::vector<double>, var, empty,
                   empty>
    type_vv_real_real_real_real_126;
typedef std::tuple<double, std::vector<var>, std::vector<double>,
                   std::vector<var>, empty, empty>
    type_vv_real_real_real_real_127;
typedef std::tuple<
    double, std::vector<var>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_128;
typedef std::tuple<double, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty,
                   empty>
    type_vv_real_real_real_real_129;
typedef std::tuple<double, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty, empty>
    type_vv_real_real_real_real_130;
typedef std::tuple<double, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_131;
typedef std::tuple<double, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty, empty>
    type_vv_real_real_real_real_132;
typedef std::tuple<double, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty, empty>
    type_vv_real_real_real_real_133;
typedef std::tuple<
    double, std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_134;
typedef std::tuple<double, std::vector<var>, var, double, empty, empty>
    type_vv_real_real_real_real_135;
typedef std::tuple<double, std::vector<var>, var, std::vector<double>, empty,
                   empty>
    type_vv_real_real_real_real_136;
typedef std::tuple<double, std::vector<var>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_137;
typedef std::tuple<double, std::vector<var>, var, var, empty, empty>
    type_vv_real_real_real_real_138;
typedef std::tuple<double, std::vector<var>, var, std::vector<var>, empty,
                   empty>
    type_vv_real_real_real_real_139;
typedef std::tuple<
    double, std::vector<var>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_140;
typedef std::tuple<double, std::vector<var>, std::vector<var>, double, empty,
                   empty>
    type_vv_real_real_real_real_141;
typedef std::tuple<double, std::vector<var>, std::vector<var>,
                   std::vector<double>, empty, empty>
    type_vv_real_real_real_real_142;
typedef std::tuple<double, std::vector<var>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_143;
typedef std::tuple<double, std::vector<var>, std::vector<var>, var, empty,
                   empty>
    type_vv_real_real_real_real_144;
typedef std::tuple<double, std::vector<var>, std::vector<var>, std::vector<var>,
                   empty, empty>
    type_vv_real_real_real_real_145;
typedef std::tuple<
    double, std::vector<var>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_146;
typedef std::tuple<
    double, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty, empty>
    type_vv_real_real_real_real_147;
typedef std::tuple<
    double, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty, empty>
    type_vv_real_real_real_real_148;
typedef std::tuple<
    double, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_149;
typedef std::tuple<
    double, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty,
    empty>
    type_vv_real_real_real_real_150;
typedef std::tuple<
    double, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty, empty>
    type_vv_real_real_real_real_151;
typedef std::tuple<
    double, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_152;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    double, double, empty, empty>
    type_vv_real_real_real_real_153;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    double, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_154;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    double, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_155;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    double, var, empty, empty>
    type_vv_real_real_real_real_156;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    double, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_157;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    empty, empty>
    type_vv_real_real_real_real_158;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, double, empty, empty>
    type_vv_real_real_real_real_159;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_160;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_161;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, var, empty, empty>
    type_vv_real_real_real_real_162;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_163;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_164;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty, empty>
    type_vv_real_real_real_real_165;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_166;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_167;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty, empty>
    type_vv_real_real_real_real_168;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_169;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_170;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    var, double, empty, empty>
    type_vv_real_real_real_real_171;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    var, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_172;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    var, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_173;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    var, var, empty, empty>
    type_vv_real_real_real_real_174;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    var, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_175;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_176;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, double, empty, empty>
    type_vv_real_real_real_real_177;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_178;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_179;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, var, empty, empty>
    type_vv_real_real_real_real_180;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_181;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_182;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty, empty>
    type_vv_real_real_real_real_183;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty, empty>
    type_vv_real_real_real_real_184;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_185;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty,
    empty>
    type_vv_real_real_real_real_186;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty, empty>
    type_vv_real_real_real_real_187;
typedef std::tuple<
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_188;
typedef std::tuple<std::vector<double>, double, double, var, empty, empty>
    type_vv_real_real_real_real_189;
typedef std::tuple<std::vector<double>, double, double, std::vector<var>, empty,
                   empty>
    type_vv_real_real_real_real_190;
typedef std::tuple<
    std::vector<double>, double, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_191;
typedef std::tuple<std::vector<double>, double, std::vector<double>, var, empty,
                   empty>
    type_vv_real_real_real_real_192;
typedef std::tuple<std::vector<double>, double, std::vector<double>,
                   std::vector<var>, empty, empty>
    type_vv_real_real_real_real_193;
typedef std::tuple<
    std::vector<double>, double, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_194;
typedef std::tuple<std::vector<double>, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty, empty>
    type_vv_real_real_real_real_195;
typedef std::tuple<std::vector<double>, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty, empty>
    type_vv_real_real_real_real_196;
typedef std::tuple<
    std::vector<double>, double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_197;
typedef std::tuple<std::vector<double>, double, var, double, empty, empty>
    type_vv_real_real_real_real_198;
typedef std::tuple<std::vector<double>, double, var, std::vector<double>, empty,
                   empty>
    type_vv_real_real_real_real_199;
typedef std::tuple<std::vector<double>, double, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_200;
typedef std::tuple<std::vector<double>, double, var, var, empty, empty>
    type_vv_real_real_real_real_201;
typedef std::tuple<std::vector<double>, double, var, std::vector<var>, empty,
                   empty>
    type_vv_real_real_real_real_202;
typedef std::tuple<
    std::vector<double>, double, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_203;
typedef std::tuple<std::vector<double>, double, std::vector<var>, double, empty,
                   empty>
    type_vv_real_real_real_real_204;
typedef std::tuple<std::vector<double>, double, std::vector<var>,
                   std::vector<double>, empty, empty>
    type_vv_real_real_real_real_205;
typedef std::tuple<std::vector<double>, double, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_206;
typedef std::tuple<std::vector<double>, double, std::vector<var>, var, empty,
                   empty>
    type_vv_real_real_real_real_207;
typedef std::tuple<std::vector<double>, double, std::vector<var>,
                   std::vector<var>, empty, empty>
    type_vv_real_real_real_real_208;
typedef std::tuple<
    std::vector<double>, double, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_209;
typedef std::tuple<
    std::vector<double>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty, empty>
    type_vv_real_real_real_real_210;
typedef std::tuple<
    std::vector<double>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty, empty>
    type_vv_real_real_real_real_211;
typedef std::tuple<
    std::vector<double>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_212;
typedef std::tuple<
    std::vector<double>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty,
    empty>
    type_vv_real_real_real_real_213;
typedef std::tuple<
    std::vector<double>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty, empty>
    type_vv_real_real_real_real_214;
typedef std::tuple<
    std::vector<double>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_215;
typedef std::tuple<std::vector<double>, std::vector<double>, double, var, empty,
                   empty>
    type_vv_real_real_real_real_216;
typedef std::tuple<std::vector<double>, std::vector<double>, double,
                   std::vector<var>, empty, empty>
    type_vv_real_real_real_real_217;
typedef std::tuple<
    std::vector<double>, std::vector<double>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_218;
typedef std::tuple<std::vector<double>, std::vector<double>,
                   std::vector<double>, var, empty, empty>
    type_vv_real_real_real_real_219;
typedef std::tuple<std::vector<double>, std::vector<double>,
                   std::vector<double>, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_220;
typedef std::tuple<
    std::vector<double>, std::vector<double>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_221;
typedef std::tuple<std::vector<double>, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty, empty>
    type_vv_real_real_real_real_222;
typedef std::tuple<std::vector<double>, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty, empty>
    type_vv_real_real_real_real_223;
typedef std::tuple<
    std::vector<double>, std::vector<double>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_224;
typedef std::tuple<std::vector<double>, std::vector<double>, var, double, empty,
                   empty>
    type_vv_real_real_real_real_225;
typedef std::tuple<std::vector<double>, std::vector<double>, var,
                   std::vector<double>, empty, empty>
    type_vv_real_real_real_real_226;
typedef std::tuple<std::vector<double>, std::vector<double>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_227;
typedef std::tuple<std::vector<double>, std::vector<double>, var, var, empty,
                   empty>
    type_vv_real_real_real_real_228;
typedef std::tuple<std::vector<double>, std::vector<double>, var,
                   std::vector<var>, empty, empty>
    type_vv_real_real_real_real_229;
typedef std::tuple<
    std::vector<double>, std::vector<double>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_230;
typedef std::tuple<std::vector<double>, std::vector<double>, std::vector<var>,
                   double, empty, empty>
    type_vv_real_real_real_real_231;
typedef std::tuple<std::vector<double>, std::vector<double>, std::vector<var>,
                   std::vector<double>, empty, empty>
    type_vv_real_real_real_real_232;
typedef std::tuple<std::vector<double>, std::vector<double>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_233;
typedef std::tuple<std::vector<double>, std::vector<double>, std::vector<var>,
                   var, empty, empty>
    type_vv_real_real_real_real_234;
typedef std::tuple<std::vector<double>, std::vector<double>, std::vector<var>,
                   std::vector<var>, empty, empty>
    type_vv_real_real_real_real_235;
typedef std::tuple<
    std::vector<double>, std::vector<double>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_236;
typedef std::tuple<
    std::vector<double>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty, empty>
    type_vv_real_real_real_real_237;
typedef std::tuple<
    std::vector<double>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty, empty>
    type_vv_real_real_real_real_238;
typedef std::tuple<
    std::vector<double>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_239;
typedef std::tuple<
    std::vector<double>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty,
    empty>
    type_vv_real_real_real_real_240;
typedef std::tuple<
    std::vector<double>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty, empty>
    type_vv_real_real_real_real_241;
typedef std::tuple<
    std::vector<double>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_242;
typedef std::tuple<std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, double, var, empty,
                   empty>
    type_vv_real_real_real_real_243;
typedef std::tuple<std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
                   std::vector<var>, empty, empty>
    type_vv_real_real_real_real_244;
typedef std::tuple<
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_245;
typedef std::tuple<std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, var, empty, empty>
    type_vv_real_real_real_real_246;
typedef std::tuple<std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_247;
typedef std::tuple<
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_248;
typedef std::tuple<std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty, empty>
    type_vv_real_real_real_real_249;
typedef std::tuple<
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_250;
typedef std::tuple<
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_251;
typedef std::tuple<std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, double, empty,
                   empty>
    type_vv_real_real_real_real_252;
typedef std::tuple<std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   std::vector<double>, empty, empty>
    type_vv_real_real_real_real_253;
typedef std::tuple<std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_254;
typedef std::tuple<std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, var, empty,
                   empty>
    type_vv_real_real_real_real_255;
typedef std::tuple<std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   std::vector<var>, empty, empty>
    type_vv_real_real_real_real_256;
typedef std::tuple<
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_257;
typedef std::tuple<std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   double, empty, empty>
    type_vv_real_real_real_real_258;
typedef std::tuple<std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<double>, empty, empty>
    type_vv_real_real_real_real_259;
typedef std::tuple<std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_260;
typedef std::tuple<std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   var, empty, empty>
    type_vv_real_real_real_real_261;
typedef std::tuple<std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<var>, empty, empty>
    type_vv_real_real_real_real_262;
typedef std::tuple<
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_263;
typedef std::tuple<
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty, empty>
    type_vv_real_real_real_real_264;
typedef std::tuple<
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty, empty>
    type_vv_real_real_real_real_265;
typedef std::tuple<
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_266;
typedef std::tuple<
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty,
    empty>
    type_vv_real_real_real_real_267;
typedef std::tuple<
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty, empty>
    type_vv_real_real_real_real_268;
typedef std::tuple<
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_269;
typedef std::tuple<std::vector<double>, var, double, double, empty, empty>
    type_vv_real_real_real_real_270;
typedef std::tuple<std::vector<double>, var, double, std::vector<double>, empty,
                   empty>
    type_vv_real_real_real_real_271;
typedef std::tuple<std::vector<double>, var, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_272;
typedef std::tuple<std::vector<double>, var, double, var, empty, empty>
    type_vv_real_real_real_real_273;
typedef std::tuple<std::vector<double>, var, double, std::vector<var>, empty,
                   empty>
    type_vv_real_real_real_real_274;
typedef std::tuple<
    std::vector<double>, var, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_275;
typedef std::tuple<std::vector<double>, var, std::vector<double>, double, empty,
                   empty>
    type_vv_real_real_real_real_276;
typedef std::tuple<std::vector<double>, var, std::vector<double>,
                   std::vector<double>, empty, empty>
    type_vv_real_real_real_real_277;
typedef std::tuple<std::vector<double>, var, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_278;
typedef std::tuple<std::vector<double>, var, std::vector<double>, var, empty,
                   empty>
    type_vv_real_real_real_real_279;
typedef std::tuple<std::vector<double>, var, std::vector<double>,
                   std::vector<var>, empty, empty>
    type_vv_real_real_real_real_280;
typedef std::tuple<
    std::vector<double>, var, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_281;
typedef std::tuple<std::vector<double>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty,
                   empty>
    type_vv_real_real_real_real_282;
typedef std::tuple<std::vector<double>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty, empty>
    type_vv_real_real_real_real_283;
typedef std::tuple<std::vector<double>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_284;
typedef std::tuple<std::vector<double>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty, empty>
    type_vv_real_real_real_real_285;
typedef std::tuple<std::vector<double>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty, empty>
    type_vv_real_real_real_real_286;
typedef std::tuple<
    std::vector<double>, var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_287;
typedef std::tuple<std::vector<double>, var, var, double, empty, empty>
    type_vv_real_real_real_real_288;
typedef std::tuple<std::vector<double>, var, var, std::vector<double>, empty,
                   empty>
    type_vv_real_real_real_real_289;
typedef std::tuple<std::vector<double>, var, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_290;
typedef std::tuple<std::vector<double>, var, var, var, empty, empty>
    type_vv_real_real_real_real_291;
typedef std::tuple<std::vector<double>, var, var, std::vector<var>, empty,
                   empty>
    type_vv_real_real_real_real_292;
typedef std::tuple<
    std::vector<double>, var, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_293;
typedef std::tuple<std::vector<double>, var, std::vector<var>, double, empty,
                   empty>
    type_vv_real_real_real_real_294;
typedef std::tuple<std::vector<double>, var, std::vector<var>,
                   std::vector<double>, empty, empty>
    type_vv_real_real_real_real_295;
typedef std::tuple<std::vector<double>, var, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_296;
typedef std::tuple<std::vector<double>, var, std::vector<var>, var, empty,
                   empty>
    type_vv_real_real_real_real_297;
typedef std::tuple<std::vector<double>, var, std::vector<var>, std::vector<var>,
                   empty, empty>
    type_vv_real_real_real_real_298;
typedef std::tuple<
    std::vector<double>, var, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_299;
typedef std::tuple<
    std::vector<double>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty, empty>
    type_vv_real_real_real_real_300;
typedef std::tuple<
    std::vector<double>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty, empty>
    type_vv_real_real_real_real_301;
typedef std::tuple<
    std::vector<double>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_302;
typedef std::tuple<
    std::vector<double>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty,
    empty>
    type_vv_real_real_real_real_303;
typedef std::tuple<
    std::vector<double>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty, empty>
    type_vv_real_real_real_real_304;
typedef std::tuple<
    std::vector<double>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_305;
typedef std::tuple<std::vector<double>, std::vector<var>, double, double, empty,
                   empty>
    type_vv_real_real_real_real_306;
typedef std::tuple<std::vector<double>, std::vector<var>, double,
                   std::vector<double>, empty, empty>
    type_vv_real_real_real_real_307;
typedef std::tuple<std::vector<double>, std::vector<var>, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_308;
typedef std::tuple<std::vector<double>, std::vector<var>, double, var, empty,
                   empty>
    type_vv_real_real_real_real_309;
typedef std::tuple<std::vector<double>, std::vector<var>, double,
                   std::vector<var>, empty, empty>
    type_vv_real_real_real_real_310;
typedef std::tuple<
    std::vector<double>, std::vector<var>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_311;
typedef std::tuple<std::vector<double>, std::vector<var>, std::vector<double>,
                   double, empty, empty>
    type_vv_real_real_real_real_312;
typedef std::tuple<std::vector<double>, std::vector<var>, std::vector<double>,
                   std::vector<double>, empty, empty>
    type_vv_real_real_real_real_313;
typedef std::tuple<std::vector<double>, std::vector<var>, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_314;
typedef std::tuple<std::vector<double>, std::vector<var>, std::vector<double>,
                   var, empty, empty>
    type_vv_real_real_real_real_315;
typedef std::tuple<std::vector<double>, std::vector<var>, std::vector<double>,
                   std::vector<var>, empty, empty>
    type_vv_real_real_real_real_316;
typedef std::tuple<
    std::vector<double>, std::vector<var>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_317;
typedef std::tuple<std::vector<double>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty,
                   empty>
    type_vv_real_real_real_real_318;
typedef std::tuple<std::vector<double>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty, empty>
    type_vv_real_real_real_real_319;
typedef std::tuple<std::vector<double>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_320;
typedef std::tuple<std::vector<double>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty, empty>
    type_vv_real_real_real_real_321;
typedef std::tuple<std::vector<double>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty, empty>
    type_vv_real_real_real_real_322;
typedef std::tuple<
    std::vector<double>, std::vector<var>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_323;
typedef std::tuple<std::vector<double>, std::vector<var>, var, double, empty,
                   empty>
    type_vv_real_real_real_real_324;
typedef std::tuple<std::vector<double>, std::vector<var>, var,
                   std::vector<double>, empty, empty>
    type_vv_real_real_real_real_325;
typedef std::tuple<std::vector<double>, std::vector<var>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_326;
typedef std::tuple<std::vector<double>, std::vector<var>, var, var, empty,
                   empty>
    type_vv_real_real_real_real_327;
typedef std::tuple<std::vector<double>, std::vector<var>, var, std::vector<var>,
                   empty, empty>
    type_vv_real_real_real_real_328;
typedef std::tuple<
    std::vector<double>, std::vector<var>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_329;
typedef std::tuple<std::vector<double>, std::vector<var>, std::vector<var>,
                   double, empty, empty>
    type_vv_real_real_real_real_330;
typedef std::tuple<std::vector<double>, std::vector<var>, std::vector<var>,
                   std::vector<double>, empty, empty>
    type_vv_real_real_real_real_331;
typedef std::tuple<std::vector<double>, std::vector<var>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_332;
typedef std::tuple<std::vector<double>, std::vector<var>, std::vector<var>, var,
                   empty, empty>
    type_vv_real_real_real_real_333;
typedef std::tuple<std::vector<double>, std::vector<var>, std::vector<var>,
                   std::vector<var>, empty, empty>
    type_vv_real_real_real_real_334;
typedef std::tuple<
    std::vector<double>, std::vector<var>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_335;
typedef std::tuple<
    std::vector<double>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty, empty>
    type_vv_real_real_real_real_336;
typedef std::tuple<
    std::vector<double>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty, empty>
    type_vv_real_real_real_real_337;
typedef std::tuple<
    std::vector<double>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_338;
typedef std::tuple<
    std::vector<double>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty,
    empty>
    type_vv_real_real_real_real_339;
typedef std::tuple<
    std::vector<double>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty, empty>
    type_vv_real_real_real_real_340;
typedef std::tuple<
    std::vector<double>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_341;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    double, empty, empty>
    type_vv_real_real_real_real_342;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    std::vector<double>, empty, empty>
    type_vv_real_real_real_real_343;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_344;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    var, empty, empty>
    type_vv_real_real_real_real_345;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    std::vector<var>, empty, empty>
    type_vv_real_real_real_real_346;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_347;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, double, empty, empty>
    type_vv_real_real_real_real_348;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_349;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_350;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, var, empty, empty>
    type_vv_real_real_real_real_351;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_352;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_353;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty, empty>
    type_vv_real_real_real_real_354;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_355;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_356;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty, empty>
    type_vv_real_real_real_real_357;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_358;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_359;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    double, empty, empty>
    type_vv_real_real_real_real_360;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<double>, empty, empty>
    type_vv_real_real_real_real_361;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_362;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, var,
    empty, empty>
    type_vv_real_real_real_real_363;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<var>, empty, empty>
    type_vv_real_real_real_real_364;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_365;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, double, empty, empty>
    type_vv_real_real_real_real_366;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_367;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_368;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, var, empty, empty>
    type_vv_real_real_real_real_369;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_370;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_371;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty, empty>
    type_vv_real_real_real_real_372;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty, empty>
    type_vv_real_real_real_real_373;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_374;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty,
    empty>
    type_vv_real_real_real_real_375;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty, empty>
    type_vv_real_real_real_real_376;
typedef std::tuple<
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_377;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double, double,
                   var, empty, empty>
    type_vv_real_real_real_real_378;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double, double,
                   std::vector<var>, empty, empty>
    type_vv_real_real_real_real_379;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_380;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
                   std::vector<double>, var, empty, empty>
    type_vv_real_real_real_real_381;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
                   std::vector<double>, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_382;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_383;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty, empty>
    type_vv_real_real_real_real_384;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty, empty>
    type_vv_real_real_real_real_385;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_386;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double, var,
                   double, empty, empty>
    type_vv_real_real_real_real_387;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double, var,
                   std::vector<double>, empty, empty>
    type_vv_real_real_real_real_388;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_389;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double, var, var,
                   empty, empty>
    type_vv_real_real_real_real_390;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double, var,
                   std::vector<var>, empty, empty>
    type_vv_real_real_real_real_391;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_392;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
                   std::vector<var>, double, empty, empty>
    type_vv_real_real_real_real_393;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
                   std::vector<var>, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_394;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
                   std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   empty, empty>
    type_vv_real_real_real_real_395;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
                   std::vector<var>, var, empty, empty>
    type_vv_real_real_real_real_396;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
                   std::vector<var>, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_397;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_398;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty, empty>
    type_vv_real_real_real_real_399;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty, empty>
    type_vv_real_real_real_real_400;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_401;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty,
    empty>
    type_vv_real_real_real_real_402;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty, empty>
    type_vv_real_real_real_real_403;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_404;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, double, var, empty, empty>
    type_vv_real_real_real_real_405;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, double, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_406;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_407;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, std::vector<double>, var, empty, empty>
    type_vv_real_real_real_real_408;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, std::vector<double>, std::vector<var>,
                   empty, empty>
    type_vv_real_real_real_real_409;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>,
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_410;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty, empty>
    type_vv_real_real_real_real_411;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_412;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_413;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, var, double, empty, empty>
    type_vv_real_real_real_real_414;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, var, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_415;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_416;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, var, var, empty, empty>
    type_vv_real_real_real_real_417;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, var, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_418;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_419;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, std::vector<var>, double, empty, empty>
    type_vv_real_real_real_real_420;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, std::vector<var>, std::vector<double>,
                   empty, empty>
    type_vv_real_real_real_real_421;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_422;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, std::vector<var>, var, empty, empty>
    type_vv_real_real_real_real_423;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, std::vector<var>, std::vector<var>,
                   empty, empty>
    type_vv_real_real_real_real_424;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>,
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_425;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty, empty>
    type_vv_real_real_real_real_426;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty, empty>
    type_vv_real_real_real_real_427;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_428;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty,
    empty>
    type_vv_real_real_real_real_429;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty, empty>
    type_vv_real_real_real_real_430;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_431;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, double, var, empty,
                   empty>
    type_vv_real_real_real_real_432;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
                   std::vector<var>, empty, empty>
    type_vv_real_real_real_real_433;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_434;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, var, empty, empty>
    type_vv_real_real_real_real_435;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_436;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_437;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty, empty>
    type_vv_real_real_real_real_438;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty, empty>
    type_vv_real_real_real_real_439;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_440;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, double, empty,
                   empty>
    type_vv_real_real_real_real_441;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   std::vector<double>, empty, empty>
    type_vv_real_real_real_real_442;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_443;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, var, empty,
                   empty>
    type_vv_real_real_real_real_444;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   std::vector<var>, empty, empty>
    type_vv_real_real_real_real_445;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_446;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   double, empty, empty>
    type_vv_real_real_real_real_447;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<double>, empty, empty>
    type_vv_real_real_real_real_448;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_449;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   var, empty, empty>
    type_vv_real_real_real_real_450;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<var>, empty, empty>
    type_vv_real_real_real_real_451;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_452;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty, empty>
    type_vv_real_real_real_real_453;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty, empty>
    type_vv_real_real_real_real_454;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_455;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty,
    empty>
    type_vv_real_real_real_real_456;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty, empty>
    type_vv_real_real_real_real_457;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_458;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var, double,
                   double, empty, empty>
    type_vv_real_real_real_real_459;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var, double,
                   std::vector<double>, empty, empty>
    type_vv_real_real_real_real_460;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_461;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var, double, var,
                   empty, empty>
    type_vv_real_real_real_real_462;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var, double,
                   std::vector<var>, empty, empty>
    type_vv_real_real_real_real_463;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_464;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   std::vector<double>, double, empty, empty>
    type_vv_real_real_real_real_465;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   std::vector<double>, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_466;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_467;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   std::vector<double>, var, empty, empty>
    type_vv_real_real_real_real_468;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   std::vector<double>, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_469;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_470;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty,
                   empty>
    type_vv_real_real_real_real_471;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty, empty>
    type_vv_real_real_real_real_472;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_473;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty, empty>
    type_vv_real_real_real_real_474;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty, empty>
    type_vv_real_real_real_real_475;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_476;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var, var, double,
                   empty, empty>
    type_vv_real_real_real_real_477;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var, var,
                   std::vector<double>, empty, empty>
    type_vv_real_real_real_real_478;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_479;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var, var, var,
                   empty, empty>
    type_vv_real_real_real_real_480;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var, var,
                   std::vector<var>, empty, empty>
    type_vv_real_real_real_real_481;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_482;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   std::vector<var>, double, empty, empty>
    type_vv_real_real_real_real_483;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   std::vector<var>, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_484;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   empty, empty>
    type_vv_real_real_real_real_485;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   std::vector<var>, var, empty, empty>
    type_vv_real_real_real_real_486;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   std::vector<var>, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_487;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_488;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty, empty>
    type_vv_real_real_real_real_489;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty, empty>
    type_vv_real_real_real_real_490;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_491;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty,
    empty>
    type_vv_real_real_real_real_492;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty, empty>
    type_vv_real_real_real_real_493;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_494;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   double, double, empty, empty>
    type_vv_real_real_real_real_495;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   double, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_496;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   double, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty,
                   empty>
    type_vv_real_real_real_real_497;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   double, var, empty, empty>
    type_vv_real_real_real_real_498;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   double, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_499;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_500;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<double>, double, empty, empty>
    type_vv_real_real_real_real_501;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<double>, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_502;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_503;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<double>, var, empty, empty>
    type_vv_real_real_real_real_504;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<double>, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_505;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_506;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty,
                   empty>
    type_vv_real_real_real_real_507;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty, empty>
    type_vv_real_real_real_real_508;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_509;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty, empty>
    type_vv_real_real_real_real_510;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty, empty>
    type_vv_real_real_real_real_511;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_512;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   var, double, empty, empty>
    type_vv_real_real_real_real_513;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   var, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_514;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   var, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_515;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   var, var, empty, empty>
    type_vv_real_real_real_real_516;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   var, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_517;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_518;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<var>, double, empty, empty>
    type_vv_real_real_real_real_519;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<var>, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_520;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   empty, empty>
    type_vv_real_real_real_real_521;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<var>, var, empty, empty>
    type_vv_real_real_real_real_522;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   std::vector<var>, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_523;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_524;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty, empty>
    type_vv_real_real_real_real_525;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty, empty>
    type_vv_real_real_real_real_526;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_527;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty,
    empty>
    type_vv_real_real_real_real_528;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty, empty>
    type_vv_real_real_real_real_529;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_530;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    double, empty, empty>
    type_vv_real_real_real_real_531;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    std::vector<double>, empty, empty>
    type_vv_real_real_real_real_532;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_533;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    var, empty, empty>
    type_vv_real_real_real_real_534;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    std::vector<var>, empty, empty>
    type_vv_real_real_real_real_535;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_536;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, double, empty, empty>
    type_vv_real_real_real_real_537;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_538;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_539;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, var, empty, empty>
    type_vv_real_real_real_real_540;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_541;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_542;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty, empty>
    type_vv_real_real_real_real_543;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_544;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_545;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty, empty>
    type_vv_real_real_real_real_546;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_547;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_548;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    double, empty, empty>
    type_vv_real_real_real_real_549;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<double>, empty, empty>
    type_vv_real_real_real_real_550;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_551;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, var,
    empty, empty>
    type_vv_real_real_real_real_552;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<var>, empty, empty>
    type_vv_real_real_real_real_553;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_554;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, double, empty, empty>
    type_vv_real_real_real_real_555;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_556;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_557;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, var, empty, empty>
    type_vv_real_real_real_real_558;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_559;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_560;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty, empty>
    type_vv_real_real_real_real_561;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty, empty>
    type_vv_real_real_real_real_562;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_563;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty,
    empty>
    type_vv_real_real_real_real_564;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty, empty>
    type_vv_real_real_real_real_565;
typedef std::tuple<
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_566;
typedef std::tuple<var, double, double, double, empty, empty>
    type_vv_real_real_real_real_567;
typedef std::tuple<var, double, double, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_568;
typedef std::tuple<var, double, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_569;
typedef std::tuple<var, double, double, var, empty, empty>
    type_vv_real_real_real_real_570;
typedef std::tuple<var, double, double, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_571;
typedef std::tuple<
    var, double, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_572;
typedef std::tuple<var, double, std::vector<double>, double, empty, empty>
    type_vv_real_real_real_real_573;
typedef std::tuple<var, double, std::vector<double>, std::vector<double>, empty,
                   empty>
    type_vv_real_real_real_real_574;
typedef std::tuple<var, double, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_575;
typedef std::tuple<var, double, std::vector<double>, var, empty, empty>
    type_vv_real_real_real_real_576;
typedef std::tuple<var, double, std::vector<double>, std::vector<var>, empty,
                   empty>
    type_vv_real_real_real_real_577;
typedef std::tuple<
    var, double, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_578;
typedef std::tuple<var, double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   double, empty, empty>
    type_vv_real_real_real_real_579;
typedef std::tuple<var, double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty, empty>
    type_vv_real_real_real_real_580;
typedef std::tuple<var, double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_581;
typedef std::tuple<var, double, Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   empty, empty>
    type_vv_real_real_real_real_582;
typedef std::tuple<var, double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<var>, empty, empty>
    type_vv_real_real_real_real_583;
typedef std::tuple<
    var, double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_584;
typedef std::tuple<var, double, var, double, empty, empty>
    type_vv_real_real_real_real_585;
typedef std::tuple<var, double, var, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_586;
typedef std::tuple<var, double, var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   empty, empty>
    type_vv_real_real_real_real_587;
typedef std::tuple<var, double, var, var, empty, empty>
    type_vv_real_real_real_real_588;
typedef std::tuple<var, double, var, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_589;
typedef std::tuple<
    var, double, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_590;
typedef std::tuple<var, double, std::vector<var>, double, empty, empty>
    type_vv_real_real_real_real_591;
typedef std::tuple<var, double, std::vector<var>, std::vector<double>, empty,
                   empty>
    type_vv_real_real_real_real_592;
typedef std::tuple<var, double, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_593;
typedef std::tuple<var, double, std::vector<var>, var, empty, empty>
    type_vv_real_real_real_real_594;
typedef std::tuple<var, double, std::vector<var>, std::vector<var>, empty,
                   empty>
    type_vv_real_real_real_real_595;
typedef std::tuple<
    var, double, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_596;
typedef std::tuple<
    var, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty, empty>
    type_vv_real_real_real_real_597;
typedef std::tuple<
    var, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty, empty>
    type_vv_real_real_real_real_598;
typedef std::tuple<
    var, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_599;
typedef std::tuple<
    var, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty,
    empty>
    type_vv_real_real_real_real_600;
typedef std::tuple<
    var, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty, empty>
    type_vv_real_real_real_real_601;
typedef std::tuple<
    var, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_602;
typedef std::tuple<var, std::vector<double>, double, double, empty, empty>
    type_vv_real_real_real_real_603;
typedef std::tuple<var, std::vector<double>, double, std::vector<double>, empty,
                   empty>
    type_vv_real_real_real_real_604;
typedef std::tuple<var, std::vector<double>, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_605;
typedef std::tuple<var, std::vector<double>, double, var, empty, empty>
    type_vv_real_real_real_real_606;
typedef std::tuple<var, std::vector<double>, double, std::vector<var>, empty,
                   empty>
    type_vv_real_real_real_real_607;
typedef std::tuple<
    var, std::vector<double>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_608;
typedef std::tuple<var, std::vector<double>, std::vector<double>, double, empty,
                   empty>
    type_vv_real_real_real_real_609;
typedef std::tuple<var, std::vector<double>, std::vector<double>,
                   std::vector<double>, empty, empty>
    type_vv_real_real_real_real_610;
typedef std::tuple<var, std::vector<double>, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_611;
typedef std::tuple<var, std::vector<double>, std::vector<double>, var, empty,
                   empty>
    type_vv_real_real_real_real_612;
typedef std::tuple<var, std::vector<double>, std::vector<double>,
                   std::vector<var>, empty, empty>
    type_vv_real_real_real_real_613;
typedef std::tuple<
    var, std::vector<double>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_614;
typedef std::tuple<var, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty,
                   empty>
    type_vv_real_real_real_real_615;
typedef std::tuple<var, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty, empty>
    type_vv_real_real_real_real_616;
typedef std::tuple<var, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_617;
typedef std::tuple<var, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty, empty>
    type_vv_real_real_real_real_618;
typedef std::tuple<var, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty, empty>
    type_vv_real_real_real_real_619;
typedef std::tuple<
    var, std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_620;
typedef std::tuple<var, std::vector<double>, var, double, empty, empty>
    type_vv_real_real_real_real_621;
typedef std::tuple<var, std::vector<double>, var, std::vector<double>, empty,
                   empty>
    type_vv_real_real_real_real_622;
typedef std::tuple<var, std::vector<double>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_623;
typedef std::tuple<var, std::vector<double>, var, var, empty, empty>
    type_vv_real_real_real_real_624;
typedef std::tuple<var, std::vector<double>, var, std::vector<var>, empty,
                   empty>
    type_vv_real_real_real_real_625;
typedef std::tuple<
    var, std::vector<double>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_626;
typedef std::tuple<var, std::vector<double>, std::vector<var>, double, empty,
                   empty>
    type_vv_real_real_real_real_627;
typedef std::tuple<var, std::vector<double>, std::vector<var>,
                   std::vector<double>, empty, empty>
    type_vv_real_real_real_real_628;
typedef std::tuple<var, std::vector<double>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_629;
typedef std::tuple<var, std::vector<double>, std::vector<var>, var, empty,
                   empty>
    type_vv_real_real_real_real_630;
typedef std::tuple<var, std::vector<double>, std::vector<var>, std::vector<var>,
                   empty, empty>
    type_vv_real_real_real_real_631;
typedef std::tuple<
    var, std::vector<double>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_632;
typedef std::tuple<
    var, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty, empty>
    type_vv_real_real_real_real_633;
typedef std::tuple<
    var, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty, empty>
    type_vv_real_real_real_real_634;
typedef std::tuple<
    var, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_635;
typedef std::tuple<
    var, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty,
    empty>
    type_vv_real_real_real_real_636;
typedef std::tuple<
    var, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty, empty>
    type_vv_real_real_real_real_637;
typedef std::tuple<
    var, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_638;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
                   double, empty, empty>
    type_vv_real_real_real_real_639;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
                   std::vector<double>, empty, empty>
    type_vv_real_real_real_real_640;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_641;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>, double, var,
                   empty, empty>
    type_vv_real_real_real_real_642;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
                   std::vector<var>, empty, empty>
    type_vv_real_real_real_real_643;
typedef std::tuple<
    var, Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_644;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, double, empty, empty>
    type_vv_real_real_real_real_645;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_646;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_647;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, var, empty, empty>
    type_vv_real_real_real_real_648;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_649;
typedef std::tuple<
    var, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_650;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty,
                   empty>
    type_vv_real_real_real_real_651;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty, empty>
    type_vv_real_real_real_real_652;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_653;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty, empty>
    type_vv_real_real_real_real_654;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty, empty>
    type_vv_real_real_real_real_655;
typedef std::tuple<
    var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_656;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>, var, double,
                   empty, empty>
    type_vv_real_real_real_real_657;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   std::vector<double>, empty, empty>
    type_vv_real_real_real_real_658;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_659;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>, var, var,
                   empty, empty>
    type_vv_real_real_real_real_660;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   std::vector<var>, empty, empty>
    type_vv_real_real_real_real_661;
typedef std::tuple<
    var, Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_662;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<var>, double, empty, empty>
    type_vv_real_real_real_real_663;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<var>, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_664;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   empty, empty>
    type_vv_real_real_real_real_665;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<var>, var, empty, empty>
    type_vv_real_real_real_real_666;
typedef std::tuple<var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<var>, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_667;
typedef std::tuple<
    var, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_668;
typedef std::tuple<
    var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty, empty>
    type_vv_real_real_real_real_669;
typedef std::tuple<
    var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty, empty>
    type_vv_real_real_real_real_670;
typedef std::tuple<
    var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_671;
typedef std::tuple<
    var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty,
    empty>
    type_vv_real_real_real_real_672;
typedef std::tuple<
    var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty, empty>
    type_vv_real_real_real_real_673;
typedef std::tuple<
    var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_674;
typedef std::tuple<var, var, double, double, empty, empty>
    type_vv_real_real_real_real_675;
typedef std::tuple<var, var, double, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_676;
typedef std::tuple<var, var, double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   empty, empty>
    type_vv_real_real_real_real_677;
typedef std::tuple<var, var, double, var, empty, empty>
    type_vv_real_real_real_real_678;
typedef std::tuple<var, var, double, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_679;
typedef std::tuple<
    var, var, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_680;
typedef std::tuple<var, var, std::vector<double>, double, empty, empty>
    type_vv_real_real_real_real_681;
typedef std::tuple<var, var, std::vector<double>, std::vector<double>, empty,
                   empty>
    type_vv_real_real_real_real_682;
typedef std::tuple<var, var, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_683;
typedef std::tuple<var, var, std::vector<double>, var, empty, empty>
    type_vv_real_real_real_real_684;
typedef std::tuple<var, var, std::vector<double>, std::vector<var>, empty,
                   empty>
    type_vv_real_real_real_real_685;
typedef std::tuple<
    var, var, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_686;
typedef std::tuple<var, var, Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
                   empty, empty>
    type_vv_real_real_real_real_687;
typedef std::tuple<var, var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty, empty>
    type_vv_real_real_real_real_688;
typedef std::tuple<var, var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_689;
typedef std::tuple<var, var, Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
                   empty, empty>
    type_vv_real_real_real_real_690;
typedef std::tuple<var, var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<var>, empty, empty>
    type_vv_real_real_real_real_691;
typedef std::tuple<
    var, var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_692;
typedef std::tuple<var, var, var, double, empty, empty>
    type_vv_real_real_real_real_693;
typedef std::tuple<var, var, var, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_694;
typedef std::tuple<var, var, var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   empty, empty>
    type_vv_real_real_real_real_695;
typedef std::tuple<var, var, var, var, empty, empty>
    type_vv_real_real_real_real_696;
typedef std::tuple<var, var, var, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_697;
typedef std::tuple<
    var, var, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_698;
typedef std::tuple<var, var, std::vector<var>, double, empty, empty>
    type_vv_real_real_real_real_699;
typedef std::tuple<var, var, std::vector<var>, std::vector<double>, empty,
                   empty>
    type_vv_real_real_real_real_700;
typedef std::tuple<var, var, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_701;
typedef std::tuple<var, var, std::vector<var>, var, empty, empty>
    type_vv_real_real_real_real_702;
typedef std::tuple<var, var, std::vector<var>, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_703;
typedef std::tuple<
    var, var, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_704;
typedef std::tuple<
    var, var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    double, empty, empty>
    type_vv_real_real_real_real_705;
typedef std::tuple<
    var, var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty, empty>
    type_vv_real_real_real_real_706;
typedef std::tuple<
    var, var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_707;
typedef std::tuple<
    var, var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    var, empty, empty>
    type_vv_real_real_real_real_708;
typedef std::tuple<
    var, var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty, empty>
    type_vv_real_real_real_real_709;
typedef std::tuple<
    var, var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_710;
typedef std::tuple<var, std::vector<var>, double, double, empty, empty>
    type_vv_real_real_real_real_711;
typedef std::tuple<var, std::vector<var>, double, std::vector<double>, empty,
                   empty>
    type_vv_real_real_real_real_712;
typedef std::tuple<var, std::vector<var>, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_713;
typedef std::tuple<var, std::vector<var>, double, var, empty, empty>
    type_vv_real_real_real_real_714;
typedef std::tuple<var, std::vector<var>, double, std::vector<var>, empty,
                   empty>
    type_vv_real_real_real_real_715;
typedef std::tuple<
    var, std::vector<var>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_716;
typedef std::tuple<var, std::vector<var>, std::vector<double>, double, empty,
                   empty>
    type_vv_real_real_real_real_717;
typedef std::tuple<var, std::vector<var>, std::vector<double>,
                   std::vector<double>, empty, empty>
    type_vv_real_real_real_real_718;
typedef std::tuple<var, std::vector<var>, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_719;
typedef std::tuple<var, std::vector<var>, std::vector<double>, var, empty,
                   empty>
    type_vv_real_real_real_real_720;
typedef std::tuple<var, std::vector<var>, std::vector<double>, std::vector<var>,
                   empty, empty>
    type_vv_real_real_real_real_721;
typedef std::tuple<
    var, std::vector<var>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_722;
typedef std::tuple<var, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty,
                   empty>
    type_vv_real_real_real_real_723;
typedef std::tuple<var, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty, empty>
    type_vv_real_real_real_real_724;
typedef std::tuple<var, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_725;
typedef std::tuple<var, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty, empty>
    type_vv_real_real_real_real_726;
typedef std::tuple<var, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty, empty>
    type_vv_real_real_real_real_727;
typedef std::tuple<
    var, std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_728;
typedef std::tuple<var, std::vector<var>, var, double, empty, empty>
    type_vv_real_real_real_real_729;
typedef std::tuple<var, std::vector<var>, var, std::vector<double>, empty,
                   empty>
    type_vv_real_real_real_real_730;
typedef std::tuple<var, std::vector<var>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_731;
typedef std::tuple<var, std::vector<var>, var, var, empty, empty>
    type_vv_real_real_real_real_732;
typedef std::tuple<var, std::vector<var>, var, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_733;
typedef std::tuple<
    var, std::vector<var>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_734;
typedef std::tuple<var, std::vector<var>, std::vector<var>, double, empty,
                   empty>
    type_vv_real_real_real_real_735;
typedef std::tuple<var, std::vector<var>, std::vector<var>, std::vector<double>,
                   empty, empty>
    type_vv_real_real_real_real_736;
typedef std::tuple<var, std::vector<var>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_737;
typedef std::tuple<var, std::vector<var>, std::vector<var>, var, empty, empty>
    type_vv_real_real_real_real_738;
typedef std::tuple<var, std::vector<var>, std::vector<var>, std::vector<var>,
                   empty, empty>
    type_vv_real_real_real_real_739;
typedef std::tuple<
    var, std::vector<var>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_740;
typedef std::tuple<
    var, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty, empty>
    type_vv_real_real_real_real_741;
typedef std::tuple<
    var, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty, empty>
    type_vv_real_real_real_real_742;
typedef std::tuple<
    var, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_743;
typedef std::tuple<
    var, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty,
    empty>
    type_vv_real_real_real_real_744;
typedef std::tuple<
    var, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty, empty>
    type_vv_real_real_real_real_745;
typedef std::tuple<
    var, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_746;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    double, double, empty, empty>
    type_vv_real_real_real_real_747;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    double, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_748;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    double, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_749;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    double, var, empty, empty>
    type_vv_real_real_real_real_750;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    double, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_751;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    empty, empty>
    type_vv_real_real_real_real_752;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, double, empty, empty>
    type_vv_real_real_real_real_753;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_754;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_755;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, var, empty, empty>
    type_vv_real_real_real_real_756;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_757;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_758;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty, empty>
    type_vv_real_real_real_real_759;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_760;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_761;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty, empty>
    type_vv_real_real_real_real_762;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_763;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_764;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    double, empty, empty>
    type_vv_real_real_real_real_765;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<double>, empty, empty>
    type_vv_real_real_real_real_766;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_767;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    var, empty, empty>
    type_vv_real_real_real_real_768;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<var>, empty, empty>
    type_vv_real_real_real_real_769;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_770;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, double, empty, empty>
    type_vv_real_real_real_real_771;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_772;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_773;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, var, empty, empty>
    type_vv_real_real_real_real_774;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_775;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_776;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty, empty>
    type_vv_real_real_real_real_777;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty, empty>
    type_vv_real_real_real_real_778;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_779;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty,
    empty>
    type_vv_real_real_real_real_780;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty, empty>
    type_vv_real_real_real_real_781;
typedef std::tuple<
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_782;
typedef std::tuple<std::vector<var>, double, double, double, empty, empty>
    type_vv_real_real_real_real_783;
typedef std::tuple<std::vector<var>, double, double, std::vector<double>, empty,
                   empty>
    type_vv_real_real_real_real_784;
typedef std::tuple<std::vector<var>, double, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_785;
typedef std::tuple<std::vector<var>, double, double, var, empty, empty>
    type_vv_real_real_real_real_786;
typedef std::tuple<std::vector<var>, double, double, std::vector<var>, empty,
                   empty>
    type_vv_real_real_real_real_787;
typedef std::tuple<
    std::vector<var>, double, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_788;
typedef std::tuple<std::vector<var>, double, std::vector<double>, double, empty,
                   empty>
    type_vv_real_real_real_real_789;
typedef std::tuple<std::vector<var>, double, std::vector<double>,
                   std::vector<double>, empty, empty>
    type_vv_real_real_real_real_790;
typedef std::tuple<std::vector<var>, double, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_791;
typedef std::tuple<std::vector<var>, double, std::vector<double>, var, empty,
                   empty>
    type_vv_real_real_real_real_792;
typedef std::tuple<std::vector<var>, double, std::vector<double>,
                   std::vector<var>, empty, empty>
    type_vv_real_real_real_real_793;
typedef std::tuple<
    std::vector<var>, double, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_794;
typedef std::tuple<std::vector<var>, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty,
                   empty>
    type_vv_real_real_real_real_795;
typedef std::tuple<std::vector<var>, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty, empty>
    type_vv_real_real_real_real_796;
typedef std::tuple<std::vector<var>, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_797;
typedef std::tuple<std::vector<var>, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty, empty>
    type_vv_real_real_real_real_798;
typedef std::tuple<std::vector<var>, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty, empty>
    type_vv_real_real_real_real_799;
typedef std::tuple<
    std::vector<var>, double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_800;
typedef std::tuple<std::vector<var>, double, var, double, empty, empty>
    type_vv_real_real_real_real_801;
typedef std::tuple<std::vector<var>, double, var, std::vector<double>, empty,
                   empty>
    type_vv_real_real_real_real_802;
typedef std::tuple<std::vector<var>, double, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_803;
typedef std::tuple<std::vector<var>, double, var, var, empty, empty>
    type_vv_real_real_real_real_804;
typedef std::tuple<std::vector<var>, double, var, std::vector<var>, empty,
                   empty>
    type_vv_real_real_real_real_805;
typedef std::tuple<
    std::vector<var>, double, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_806;
typedef std::tuple<std::vector<var>, double, std::vector<var>, double, empty,
                   empty>
    type_vv_real_real_real_real_807;
typedef std::tuple<std::vector<var>, double, std::vector<var>,
                   std::vector<double>, empty, empty>
    type_vv_real_real_real_real_808;
typedef std::tuple<std::vector<var>, double, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_809;
typedef std::tuple<std::vector<var>, double, std::vector<var>, var, empty,
                   empty>
    type_vv_real_real_real_real_810;
typedef std::tuple<std::vector<var>, double, std::vector<var>, std::vector<var>,
                   empty, empty>
    type_vv_real_real_real_real_811;
typedef std::tuple<
    std::vector<var>, double, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_812;
typedef std::tuple<
    std::vector<var>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty, empty>
    type_vv_real_real_real_real_813;
typedef std::tuple<
    std::vector<var>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty, empty>
    type_vv_real_real_real_real_814;
typedef std::tuple<
    std::vector<var>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_815;
typedef std::tuple<
    std::vector<var>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty,
    empty>
    type_vv_real_real_real_real_816;
typedef std::tuple<
    std::vector<var>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty, empty>
    type_vv_real_real_real_real_817;
typedef std::tuple<
    std::vector<var>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_818;
typedef std::tuple<std::vector<var>, std::vector<double>, double, double, empty,
                   empty>
    type_vv_real_real_real_real_819;
typedef std::tuple<std::vector<var>, std::vector<double>, double,
                   std::vector<double>, empty, empty>
    type_vv_real_real_real_real_820;
typedef std::tuple<std::vector<var>, std::vector<double>, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_821;
typedef std::tuple<std::vector<var>, std::vector<double>, double, var, empty,
                   empty>
    type_vv_real_real_real_real_822;
typedef std::tuple<std::vector<var>, std::vector<double>, double,
                   std::vector<var>, empty, empty>
    type_vv_real_real_real_real_823;
typedef std::tuple<
    std::vector<var>, std::vector<double>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_824;
typedef std::tuple<std::vector<var>, std::vector<double>, std::vector<double>,
                   double, empty, empty>
    type_vv_real_real_real_real_825;
typedef std::tuple<std::vector<var>, std::vector<double>, std::vector<double>,
                   std::vector<double>, empty, empty>
    type_vv_real_real_real_real_826;
typedef std::tuple<std::vector<var>, std::vector<double>, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_827;
typedef std::tuple<std::vector<var>, std::vector<double>, std::vector<double>,
                   var, empty, empty>
    type_vv_real_real_real_real_828;
typedef std::tuple<std::vector<var>, std::vector<double>, std::vector<double>,
                   std::vector<var>, empty, empty>
    type_vv_real_real_real_real_829;
typedef std::tuple<
    std::vector<var>, std::vector<double>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_830;
typedef std::tuple<std::vector<var>, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty,
                   empty>
    type_vv_real_real_real_real_831;
typedef std::tuple<std::vector<var>, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty, empty>
    type_vv_real_real_real_real_832;
typedef std::tuple<std::vector<var>, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_833;
typedef std::tuple<std::vector<var>, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty, empty>
    type_vv_real_real_real_real_834;
typedef std::tuple<std::vector<var>, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty, empty>
    type_vv_real_real_real_real_835;
typedef std::tuple<
    std::vector<var>, std::vector<double>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_836;
typedef std::tuple<std::vector<var>, std::vector<double>, var, double, empty,
                   empty>
    type_vv_real_real_real_real_837;
typedef std::tuple<std::vector<var>, std::vector<double>, var,
                   std::vector<double>, empty, empty>
    type_vv_real_real_real_real_838;
typedef std::tuple<std::vector<var>, std::vector<double>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_839;
typedef std::tuple<std::vector<var>, std::vector<double>, var, var, empty,
                   empty>
    type_vv_real_real_real_real_840;
typedef std::tuple<std::vector<var>, std::vector<double>, var, std::vector<var>,
                   empty, empty>
    type_vv_real_real_real_real_841;
typedef std::tuple<
    std::vector<var>, std::vector<double>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_842;
typedef std::tuple<std::vector<var>, std::vector<double>, std::vector<var>,
                   double, empty, empty>
    type_vv_real_real_real_real_843;
typedef std::tuple<std::vector<var>, std::vector<double>, std::vector<var>,
                   std::vector<double>, empty, empty>
    type_vv_real_real_real_real_844;
typedef std::tuple<std::vector<var>, std::vector<double>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_845;
typedef std::tuple<std::vector<var>, std::vector<double>, std::vector<var>, var,
                   empty, empty>
    type_vv_real_real_real_real_846;
typedef std::tuple<std::vector<var>, std::vector<double>, std::vector<var>,
                   std::vector<var>, empty, empty>
    type_vv_real_real_real_real_847;
typedef std::tuple<
    std::vector<var>, std::vector<double>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_848;
typedef std::tuple<
    std::vector<var>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty, empty>
    type_vv_real_real_real_real_849;
typedef std::tuple<
    std::vector<var>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty, empty>
    type_vv_real_real_real_real_850;
typedef std::tuple<
    std::vector<var>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_851;
typedef std::tuple<
    std::vector<var>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty,
    empty>
    type_vv_real_real_real_real_852;
typedef std::tuple<
    std::vector<var>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty, empty>
    type_vv_real_real_real_real_853;
typedef std::tuple<
    std::vector<var>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_854;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   double, double, empty, empty>
    type_vv_real_real_real_real_855;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   double, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_856;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   double, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty,
                   empty>
    type_vv_real_real_real_real_857;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   double, var, empty, empty>
    type_vv_real_real_real_real_858;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   double, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_859;
typedef std::tuple<
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_860;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, double, empty, empty>
    type_vv_real_real_real_real_861;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_862;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_863;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, var, empty, empty>
    type_vv_real_real_real_real_864;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_865;
typedef std::tuple<
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_866;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty,
                   empty>
    type_vv_real_real_real_real_867;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty, empty>
    type_vv_real_real_real_real_868;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_869;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty, empty>
    type_vv_real_real_real_real_870;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty, empty>
    type_vv_real_real_real_real_871;
typedef std::tuple<
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_872;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   var, double, empty, empty>
    type_vv_real_real_real_real_873;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   var, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_874;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   var, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_875;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   var, var, empty, empty>
    type_vv_real_real_real_real_876;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   var, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_877;
typedef std::tuple<
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_878;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<var>, double, empty, empty>
    type_vv_real_real_real_real_879;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<var>, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_880;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   empty, empty>
    type_vv_real_real_real_real_881;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<var>, var, empty, empty>
    type_vv_real_real_real_real_882;
typedef std::tuple<std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<var>, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_883;
typedef std::tuple<
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_884;
typedef std::tuple<
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty, empty>
    type_vv_real_real_real_real_885;
typedef std::tuple<
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty, empty>
    type_vv_real_real_real_real_886;
typedef std::tuple<
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_887;
typedef std::tuple<
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty,
    empty>
    type_vv_real_real_real_real_888;
typedef std::tuple<
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty, empty>
    type_vv_real_real_real_real_889;
typedef std::tuple<
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_890;
typedef std::tuple<std::vector<var>, var, double, double, empty, empty>
    type_vv_real_real_real_real_891;
typedef std::tuple<std::vector<var>, var, double, std::vector<double>, empty,
                   empty>
    type_vv_real_real_real_real_892;
typedef std::tuple<std::vector<var>, var, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_893;
typedef std::tuple<std::vector<var>, var, double, var, empty, empty>
    type_vv_real_real_real_real_894;
typedef std::tuple<std::vector<var>, var, double, std::vector<var>, empty,
                   empty>
    type_vv_real_real_real_real_895;
typedef std::tuple<
    std::vector<var>, var, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_896;
typedef std::tuple<std::vector<var>, var, std::vector<double>, double, empty,
                   empty>
    type_vv_real_real_real_real_897;
typedef std::tuple<std::vector<var>, var, std::vector<double>,
                   std::vector<double>, empty, empty>
    type_vv_real_real_real_real_898;
typedef std::tuple<std::vector<var>, var, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_899;
typedef std::tuple<std::vector<var>, var, std::vector<double>, var, empty,
                   empty>
    type_vv_real_real_real_real_900;
typedef std::tuple<std::vector<var>, var, std::vector<double>, std::vector<var>,
                   empty, empty>
    type_vv_real_real_real_real_901;
typedef std::tuple<
    std::vector<var>, var, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_902;
typedef std::tuple<std::vector<var>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty,
                   empty>
    type_vv_real_real_real_real_903;
typedef std::tuple<std::vector<var>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty, empty>
    type_vv_real_real_real_real_904;
typedef std::tuple<std::vector<var>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_905;
typedef std::tuple<std::vector<var>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty, empty>
    type_vv_real_real_real_real_906;
typedef std::tuple<std::vector<var>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty, empty>
    type_vv_real_real_real_real_907;
typedef std::tuple<
    std::vector<var>, var, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_908;
typedef std::tuple<std::vector<var>, var, var, double, empty, empty>
    type_vv_real_real_real_real_909;
typedef std::tuple<std::vector<var>, var, var, std::vector<double>, empty,
                   empty>
    type_vv_real_real_real_real_910;
typedef std::tuple<std::vector<var>, var, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_911;
typedef std::tuple<std::vector<var>, var, var, var, empty, empty>
    type_vv_real_real_real_real_912;
typedef std::tuple<std::vector<var>, var, var, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_913;
typedef std::tuple<
    std::vector<var>, var, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_914;
typedef std::tuple<std::vector<var>, var, std::vector<var>, double, empty,
                   empty>
    type_vv_real_real_real_real_915;
typedef std::tuple<std::vector<var>, var, std::vector<var>, std::vector<double>,
                   empty, empty>
    type_vv_real_real_real_real_916;
typedef std::tuple<std::vector<var>, var, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_917;
typedef std::tuple<std::vector<var>, var, std::vector<var>, var, empty, empty>
    type_vv_real_real_real_real_918;
typedef std::tuple<std::vector<var>, var, std::vector<var>, std::vector<var>,
                   empty, empty>
    type_vv_real_real_real_real_919;
typedef std::tuple<
    std::vector<var>, var, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_920;
typedef std::tuple<
    std::vector<var>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty, empty>
    type_vv_real_real_real_real_921;
typedef std::tuple<
    std::vector<var>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty, empty>
    type_vv_real_real_real_real_922;
typedef std::tuple<
    std::vector<var>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_923;
typedef std::tuple<
    std::vector<var>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty,
    empty>
    type_vv_real_real_real_real_924;
typedef std::tuple<
    std::vector<var>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty, empty>
    type_vv_real_real_real_real_925;
typedef std::tuple<
    std::vector<var>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_926;
typedef std::tuple<std::vector<var>, std::vector<var>, double, double, empty,
                   empty>
    type_vv_real_real_real_real_927;
typedef std::tuple<std::vector<var>, std::vector<var>, double,
                   std::vector<double>, empty, empty>
    type_vv_real_real_real_real_928;
typedef std::tuple<std::vector<var>, std::vector<var>, double,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_929;
typedef std::tuple<std::vector<var>, std::vector<var>, double, var, empty,
                   empty>
    type_vv_real_real_real_real_930;
typedef std::tuple<std::vector<var>, std::vector<var>, double, std::vector<var>,
                   empty, empty>
    type_vv_real_real_real_real_931;
typedef std::tuple<
    std::vector<var>, std::vector<var>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_932;
typedef std::tuple<std::vector<var>, std::vector<var>, std::vector<double>,
                   double, empty, empty>
    type_vv_real_real_real_real_933;
typedef std::tuple<std::vector<var>, std::vector<var>, std::vector<double>,
                   std::vector<double>, empty, empty>
    type_vv_real_real_real_real_934;
typedef std::tuple<std::vector<var>, std::vector<var>, std::vector<double>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_935;
typedef std::tuple<std::vector<var>, std::vector<var>, std::vector<double>, var,
                   empty, empty>
    type_vv_real_real_real_real_936;
typedef std::tuple<std::vector<var>, std::vector<var>, std::vector<double>,
                   std::vector<var>, empty, empty>
    type_vv_real_real_real_real_937;
typedef std::tuple<
    std::vector<var>, std::vector<var>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_938;
typedef std::tuple<std::vector<var>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty,
                   empty>
    type_vv_real_real_real_real_939;
typedef std::tuple<std::vector<var>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   std::vector<double>, empty, empty>
    type_vv_real_real_real_real_940;
typedef std::tuple<std::vector<var>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_941;
typedef std::tuple<std::vector<var>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty, empty>
    type_vv_real_real_real_real_942;
typedef std::tuple<std::vector<var>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
                   empty, empty>
    type_vv_real_real_real_real_943;
typedef std::tuple<
    std::vector<var>, std::vector<var>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_944;
typedef std::tuple<std::vector<var>, std::vector<var>, var, double, empty,
                   empty>
    type_vv_real_real_real_real_945;
typedef std::tuple<std::vector<var>, std::vector<var>, var, std::vector<double>,
                   empty, empty>
    type_vv_real_real_real_real_946;
typedef std::tuple<std::vector<var>, std::vector<var>, var,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_947;
typedef std::tuple<std::vector<var>, std::vector<var>, var, var, empty, empty>
    type_vv_real_real_real_real_948;
typedef std::tuple<std::vector<var>, std::vector<var>, var, std::vector<var>,
                   empty, empty>
    type_vv_real_real_real_real_949;
typedef std::tuple<
    std::vector<var>, std::vector<var>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_950;
typedef std::tuple<std::vector<var>, std::vector<var>, std::vector<var>, double,
                   empty, empty>
    type_vv_real_real_real_real_951;
typedef std::tuple<std::vector<var>, std::vector<var>, std::vector<var>,
                   std::vector<double>, empty, empty>
    type_vv_real_real_real_real_952;
typedef std::tuple<std::vector<var>, std::vector<var>, std::vector<var>,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_953;
typedef std::tuple<std::vector<var>, std::vector<var>, std::vector<var>, var,
                   empty, empty>
    type_vv_real_real_real_real_954;
typedef std::tuple<std::vector<var>, std::vector<var>, std::vector<var>,
                   std::vector<var>, empty, empty>
    type_vv_real_real_real_real_955;
typedef std::tuple<
    std::vector<var>, std::vector<var>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_956;
typedef std::tuple<
    std::vector<var>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty, empty>
    type_vv_real_real_real_real_957;
typedef std::tuple<
    std::vector<var>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty, empty>
    type_vv_real_real_real_real_958;
typedef std::tuple<
    std::vector<var>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_959;
typedef std::tuple<
    std::vector<var>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty,
    empty>
    type_vv_real_real_real_real_960;
typedef std::tuple<
    std::vector<var>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty, empty>
    type_vv_real_real_real_real_961;
typedef std::tuple<
    std::vector<var>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_962;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    double, empty, empty>
    type_vv_real_real_real_real_963;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    std::vector<double>, empty, empty>
    type_vv_real_real_real_real_964;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_965;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    var, empty, empty>
    type_vv_real_real_real_real_966;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    std::vector<var>, empty, empty>
    type_vv_real_real_real_real_967;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_968;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, double, empty, empty>
    type_vv_real_real_real_real_969;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_970;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_971;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, var, empty, empty>
    type_vv_real_real_real_real_972;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_973;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_974;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty, empty>
    type_vv_real_real_real_real_975;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_976;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_977;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty, empty>
    type_vv_real_real_real_real_978;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_979;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_980;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    double, empty, empty>
    type_vv_real_real_real_real_981;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<double>, empty, empty>
    type_vv_real_real_real_real_982;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_983;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, var,
    empty, empty>
    type_vv_real_real_real_real_984;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<var>, empty, empty>
    type_vv_real_real_real_real_985;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_986;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, double, empty, empty>
    type_vv_real_real_real_real_987;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_988;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_989;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, var, empty, empty>
    type_vv_real_real_real_real_990;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_991;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_992;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty, empty>
    type_vv_real_real_real_real_993;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty, empty>
    type_vv_real_real_real_real_994;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_995;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty,
    empty>
    type_vv_real_real_real_real_996;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty, empty>
    type_vv_real_real_real_real_997;
typedef std::tuple<
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_998;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    double, double, empty, empty>
    type_vv_real_real_real_real_999;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    double, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_1000;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    double, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_1001;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    double, var, empty, empty>
    type_vv_real_real_real_real_1002;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    double, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_1003;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    empty, empty>
    type_vv_real_real_real_real_1004;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    std::vector<double>, double, empty, empty>
    type_vv_real_real_real_real_1005;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    std::vector<double>, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_1006;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_1007;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    std::vector<double>, var, empty, empty>
    type_vv_real_real_real_real_1008;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    std::vector<double>, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_1009;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_1010;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty, empty>
    type_vv_real_real_real_real_1011;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_1012;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_1013;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty, empty>
    type_vv_real_real_real_real_1014;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_1015;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_1016;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    var, double, empty, empty>
    type_vv_real_real_real_real_1017;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    var, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_1018;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    var, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_1019;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    var, var, empty, empty>
    type_vv_real_real_real_real_1020;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    var, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_1021;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    var, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_1022;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    std::vector<var>, double, empty, empty>
    type_vv_real_real_real_real_1023;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    std::vector<var>, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_1024;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_1025;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    std::vector<var>, var, empty, empty>
    type_vv_real_real_real_real_1026;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    std::vector<var>, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_1027;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_1028;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty, empty>
    type_vv_real_real_real_real_1029;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty, empty>
    type_vv_real_real_real_real_1030;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_1031;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty,
    empty>
    type_vv_real_real_real_real_1032;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty, empty>
    type_vv_real_real_real_real_1033;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_1034;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, double, double, empty, empty>
    type_vv_real_real_real_real_1035;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, double, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_1036;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, double, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    empty, empty>
    type_vv_real_real_real_real_1037;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, double, var, empty, empty>
    type_vv_real_real_real_real_1038;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, double, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_1039;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_1040;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<double>, double, empty, empty>
    type_vv_real_real_real_real_1041;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<double>, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_1042;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<double>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_1043;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<double>, var, empty, empty>
    type_vv_real_real_real_real_1044;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<double>, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_1045;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_1046;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
    empty, empty>
    type_vv_real_real_real_real_1047;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    std::vector<double>, empty, empty>
    type_vv_real_real_real_real_1048;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_1049;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty,
    empty>
    type_vv_real_real_real_real_1050;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    std::vector<var>, empty, empty>
    type_vv_real_real_real_real_1051;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_1052;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, var, double, empty, empty>
    type_vv_real_real_real_real_1053;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, var, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_1054;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, var, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty,
    empty>
    type_vv_real_real_real_real_1055;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, var, var, empty, empty>
    type_vv_real_real_real_real_1056;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, var, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_1057;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_1058;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<var>, double, empty, empty>
    type_vv_real_real_real_real_1059;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<var>, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_1060;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<var>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_1061;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<var>, var, empty, empty>
    type_vv_real_real_real_real_1062;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<var>, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_1063;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_1064;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty, empty>
    type_vv_real_real_real_real_1065;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty, empty>
    type_vv_real_real_real_real_1066;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_1067;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty,
    empty>
    type_vv_real_real_real_real_1068;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty, empty>
    type_vv_real_real_real_real_1069;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_1070;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, double, empty, empty>
    type_vv_real_real_real_real_1071;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, std::vector<double>,
    empty, empty>
    type_vv_real_real_real_real_1072;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_1073;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, var, empty, empty>
    type_vv_real_real_real_real_1074;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, std::vector<var>, empty,
    empty>
    type_vv_real_real_real_real_1075;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_1076;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, double,
    empty, empty>
    type_vv_real_real_real_real_1077;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>,
    std::vector<double>, empty, empty>
    type_vv_real_real_real_real_1078;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_1079;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, var, empty,
    empty>
    type_vv_real_real_real_real_1080;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>,
    std::vector<var>, empty, empty>
    type_vv_real_real_real_real_1081;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_1082;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty, empty>
    type_vv_real_real_real_real_1083;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_1084;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_1085;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty, empty>
    type_vv_real_real_real_real_1086;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_1087;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_1088;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, double, empty, empty>
    type_vv_real_real_real_real_1089;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, std::vector<double>, empty,
    empty>
    type_vv_real_real_real_real_1090;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_1091;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, var, empty, empty>
    type_vv_real_real_real_real_1092;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, std::vector<var>, empty,
    empty>
    type_vv_real_real_real_real_1093;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_1094;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, double, empty,
    empty>
    type_vv_real_real_real_real_1095;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
    std::vector<double>, empty, empty>
    type_vv_real_real_real_real_1096;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_1097;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, var, empty,
    empty>
    type_vv_real_real_real_real_1098;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
    std::vector<var>, empty, empty>
    type_vv_real_real_real_real_1099;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_1100;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty, empty>
    type_vv_real_real_real_real_1101;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty, empty>
    type_vv_real_real_real_real_1102;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_1103;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty,
    empty>
    type_vv_real_real_real_real_1104;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty, empty>
    type_vv_real_real_real_real_1105;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_1106;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    double, double, empty, empty>
    type_vv_real_real_real_real_1107;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    double, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_1108;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    double, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_1109;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    double, var, empty, empty>
    type_vv_real_real_real_real_1110;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    double, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_1111;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    double, stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    empty, empty>
    type_vv_real_real_real_real_1112;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<double>, double, empty, empty>
    type_vv_real_real_real_real_1113;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<double>, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_1114;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_1115;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<double>, var, empty, empty>
    type_vv_real_real_real_real_1116;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<double>, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_1117;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_1118;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty, empty>
    type_vv_real_real_real_real_1119;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_1120;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_1121;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty, empty>
    type_vv_real_real_real_real_1122;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_1123;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_1124;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, var,
    double, empty, empty>
    type_vv_real_real_real_real_1125;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, var,
    std::vector<double>, empty, empty>
    type_vv_real_real_real_real_1126;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, var,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_1127;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, var,
    var, empty, empty>
    type_vv_real_real_real_real_1128;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, var,
    std::vector<var>, empty, empty>
    type_vv_real_real_real_real_1129;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_1130;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<var>, double, empty, empty>
    type_vv_real_real_real_real_1131;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<var>, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_1132;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_1133;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<var>, var, empty, empty>
    type_vv_real_real_real_real_1134;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<var>, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_1135;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_1136;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty, empty>
    type_vv_real_real_real_real_1137;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty, empty>
    type_vv_real_real_real_real_1138;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_1139;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty,
    empty>
    type_vv_real_real_real_real_1140;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty, empty>
    type_vv_real_real_real_real_1141;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_1142;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, double, double, empty, empty>
    type_vv_real_real_real_real_1143;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, double, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_1144;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, double, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty,
    empty>
    type_vv_real_real_real_real_1145;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, double, var, empty, empty>
    type_vv_real_real_real_real_1146;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, double, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_1147;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_1148;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<double>, double, empty, empty>
    type_vv_real_real_real_real_1149;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<double>, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_1150;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<double>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_1151;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<double>, var, empty, empty>
    type_vv_real_real_real_real_1152;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<double>, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_1153;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_1154;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty,
    empty>
    type_vv_real_real_real_real_1155;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    std::vector<double>, empty, empty>
    type_vv_real_real_real_real_1156;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_1157;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty,
    empty>
    type_vv_real_real_real_real_1158;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    std::vector<var>, empty, empty>
    type_vv_real_real_real_real_1159;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_1160;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, var, double, empty, empty>
    type_vv_real_real_real_real_1161;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, var, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_1162;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, var, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty,
    empty>
    type_vv_real_real_real_real_1163;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, var, var, empty, empty>
    type_vv_real_real_real_real_1164;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, var, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_1165;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_1166;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<var>, double, empty, empty>
    type_vv_real_real_real_real_1167;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<var>, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_1168;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<var>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_1169;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<var>, var, empty, empty>
    type_vv_real_real_real_real_1170;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<var>, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_1171;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_1172;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty, empty>
    type_vv_real_real_real_real_1173;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty, empty>
    type_vv_real_real_real_real_1174;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_1175;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty,
    empty>
    type_vv_real_real_real_real_1176;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty, empty>
    type_vv_real_real_real_real_1177;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_1178;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    double, empty, empty>
    type_vv_real_real_real_real_1179;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    std::vector<double>, empty, empty>
    type_vv_real_real_real_real_1180;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_1181;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    var, empty, empty>
    type_vv_real_real_real_real_1182;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    std::vector<var>, empty, empty>
    type_vv_real_real_real_real_1183;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_1184;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, double, empty, empty>
    type_vv_real_real_real_real_1185;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_1186;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_1187;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, var, empty, empty>
    type_vv_real_real_real_real_1188;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_1189;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_1190;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty, empty>
    type_vv_real_real_real_real_1191;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_1192;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_1193;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, var, empty, empty>
    type_vv_real_real_real_real_1194;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_1195;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_1196;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    double, empty, empty>
    type_vv_real_real_real_real_1197;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<double>, empty, empty>
    type_vv_real_real_real_real_1198;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_1199;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, var,
    empty, empty>
    type_vv_real_real_real_real_1200;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    std::vector<var>, empty, empty>
    type_vv_real_real_real_real_1201;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_1202;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, double, empty, empty>
    type_vv_real_real_real_real_1203;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<double>, empty, empty>
    type_vv_real_real_real_real_1204;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_1205;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, var, empty, empty>
    type_vv_real_real_real_real_1206;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, std::vector<var>, empty, empty>
    type_vv_real_real_real_real_1207;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_1208;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, double,
    empty, empty>
    type_vv_real_real_real_real_1209;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<double>, empty, empty>
    type_vv_real_real_real_real_1210;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty>
    type_vv_real_real_real_real_1211;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, var, empty,
    empty>
    type_vv_real_real_real_real_1212;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    std::vector<var>, empty, empty>
    type_vv_real_real_real_real_1213;
typedef std::tuple<
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
    stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty,
    empty>
    type_vv_real_real_real_real_1214;
