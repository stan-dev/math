#include <gtest/gtest.h>
#include <tuple>
#include <test/prob/test_fixture_distr.hpp>
#include <test/prob/test_fixture_cdf.hpp>
#include <test/prob/test_fixture_cdf_log.hpp>
#include <test/prob/test_fixture_ccdf_log.hpp>

typedef std::tuple<double, double, fvar<fvar<var> >, empty, empty, empty> type_ffv_real_real_real_0;
typedef std::tuple<double, double, std::vector<fvar<fvar<var> >>, empty, empty, empty> type_ffv_real_real_real_1;
typedef std::tuple<double, double, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_2;
typedef std::tuple<double, std::vector<double>, fvar<fvar<var> >, empty, empty, empty> type_ffv_real_real_real_3;
typedef std::tuple<double, std::vector<double>, std::vector<fvar<fvar<var> >>, empty, empty, empty> type_ffv_real_real_real_4;
typedef std::tuple<double, std::vector<double>, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_5;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>, fvar<fvar<var> >, empty, empty, empty> type_ffv_real_real_real_6;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<fvar<fvar<var> >>, empty, empty, empty> type_ffv_real_real_real_7;
typedef std::tuple<double, Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_8;
typedef std::tuple<double, fvar<fvar<var> >, double, empty, empty, empty> type_ffv_real_real_real_9;
typedef std::tuple<double, fvar<fvar<var> >, std::vector<double>, empty, empty, empty> type_ffv_real_real_real_10;
typedef std::tuple<double, fvar<fvar<var> >, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_11;
typedef std::tuple<double, fvar<fvar<var> >, fvar<fvar<var> >, empty, empty, empty> type_ffv_real_real_real_12;
typedef std::tuple<double, fvar<fvar<var> >, std::vector<fvar<fvar<var> >>, empty, empty, empty> type_ffv_real_real_real_13;
typedef std::tuple<double, fvar<fvar<var> >, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_14;
typedef std::tuple<double, std::vector<fvar<fvar<var> >>, double, empty, empty, empty> type_ffv_real_real_real_15;
typedef std::tuple<double, std::vector<fvar<fvar<var> >>, std::vector<double>, empty, empty, empty> type_ffv_real_real_real_16;
typedef std::tuple<double, std::vector<fvar<fvar<var> >>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_17;
typedef std::tuple<double, std::vector<fvar<fvar<var> >>, fvar<fvar<var> >, empty, empty, empty> type_ffv_real_real_real_18;
typedef std::tuple<double, std::vector<fvar<fvar<var> >>, std::vector<fvar<fvar<var> >>, empty, empty, empty> type_ffv_real_real_real_19;
typedef std::tuple<double, std::vector<fvar<fvar<var> >>, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_20;
typedef std::tuple<double, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, double, empty, empty, empty> type_ffv_real_real_real_21;
typedef std::tuple<double, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, std::vector<double>, empty, empty, empty> type_ffv_real_real_real_22;
typedef std::tuple<double, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_23;
typedef std::tuple<double, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, fvar<fvar<var> >, empty, empty, empty> type_ffv_real_real_real_24;
typedef std::tuple<double, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, std::vector<fvar<fvar<var> >>, empty, empty, empty> type_ffv_real_real_real_25;
typedef std::tuple<double, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_26;
typedef std::tuple<std::vector<double>, double, fvar<fvar<var> >, empty, empty, empty> type_ffv_real_real_real_27;
typedef std::tuple<std::vector<double>, double, std::vector<fvar<fvar<var> >>, empty, empty, empty> type_ffv_real_real_real_28;
typedef std::tuple<std::vector<double>, double, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_29;
typedef std::tuple<std::vector<double>, std::vector<double>, fvar<fvar<var> >, empty, empty, empty> type_ffv_real_real_real_30;
typedef std::tuple<std::vector<double>, std::vector<double>, std::vector<fvar<fvar<var> >>, empty, empty, empty> type_ffv_real_real_real_31;
typedef std::tuple<std::vector<double>, std::vector<double>, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_32;
typedef std::tuple<std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, fvar<fvar<var> >, empty, empty, empty> type_ffv_real_real_real_33;
typedef std::tuple<std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<fvar<fvar<var> >>, empty, empty, empty> type_ffv_real_real_real_34;
typedef std::tuple<std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_35;
typedef std::tuple<std::vector<double>, fvar<fvar<var> >, double, empty, empty, empty> type_ffv_real_real_real_36;
typedef std::tuple<std::vector<double>, fvar<fvar<var> >, std::vector<double>, empty, empty, empty> type_ffv_real_real_real_37;
typedef std::tuple<std::vector<double>, fvar<fvar<var> >, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_38;
typedef std::tuple<std::vector<double>, fvar<fvar<var> >, fvar<fvar<var> >, empty, empty, empty> type_ffv_real_real_real_39;
typedef std::tuple<std::vector<double>, fvar<fvar<var> >, std::vector<fvar<fvar<var> >>, empty, empty, empty> type_ffv_real_real_real_40;
typedef std::tuple<std::vector<double>, fvar<fvar<var> >, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_41;
typedef std::tuple<std::vector<double>, std::vector<fvar<fvar<var> >>, double, empty, empty, empty> type_ffv_real_real_real_42;
typedef std::tuple<std::vector<double>, std::vector<fvar<fvar<var> >>, std::vector<double>, empty, empty, empty> type_ffv_real_real_real_43;
typedef std::tuple<std::vector<double>, std::vector<fvar<fvar<var> >>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_44;
typedef std::tuple<std::vector<double>, std::vector<fvar<fvar<var> >>, fvar<fvar<var> >, empty, empty, empty> type_ffv_real_real_real_45;
typedef std::tuple<std::vector<double>, std::vector<fvar<fvar<var> >>, std::vector<fvar<fvar<var> >>, empty, empty, empty> type_ffv_real_real_real_46;
typedef std::tuple<std::vector<double>, std::vector<fvar<fvar<var> >>, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_47;
typedef std::tuple<std::vector<double>, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, double, empty, empty, empty> type_ffv_real_real_real_48;
typedef std::tuple<std::vector<double>, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, std::vector<double>, empty, empty, empty> type_ffv_real_real_real_49;
typedef std::tuple<std::vector<double>, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_50;
typedef std::tuple<std::vector<double>, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, fvar<fvar<var> >, empty, empty, empty> type_ffv_real_real_real_51;
typedef std::tuple<std::vector<double>, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, std::vector<fvar<fvar<var> >>, empty, empty, empty> type_ffv_real_real_real_52;
typedef std::tuple<std::vector<double>, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_53;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double, fvar<fvar<var> >, empty, empty, empty> type_ffv_real_real_real_54;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double, std::vector<fvar<fvar<var> >>, empty, empty, empty> type_ffv_real_real_real_55;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, double, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_56;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, fvar<fvar<var> >, empty, empty, empty> type_ffv_real_real_real_57;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, std::vector<fvar<fvar<var> >>, empty, empty, empty> type_ffv_real_real_real_58;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_59;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Matrix<double, Eigen::Dynamic, 1>, fvar<fvar<var> >, empty, empty, empty> type_ffv_real_real_real_60;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<fvar<fvar<var> >>, empty, empty, empty> type_ffv_real_real_real_61;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_62;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, fvar<fvar<var> >, double, empty, empty, empty> type_ffv_real_real_real_63;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, fvar<fvar<var> >, std::vector<double>, empty, empty, empty> type_ffv_real_real_real_64;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, fvar<fvar<var> >, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_65;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, fvar<fvar<var> >, fvar<fvar<var> >, empty, empty, empty> type_ffv_real_real_real_66;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, fvar<fvar<var> >, std::vector<fvar<fvar<var> >>, empty, empty, empty> type_ffv_real_real_real_67;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, fvar<fvar<var> >, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_68;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<fvar<fvar<var> >>, double, empty, empty, empty> type_ffv_real_real_real_69;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<fvar<fvar<var> >>, std::vector<double>, empty, empty, empty> type_ffv_real_real_real_70;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<fvar<fvar<var> >>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_71;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<fvar<fvar<var> >>, fvar<fvar<var> >, empty, empty, empty> type_ffv_real_real_real_72;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<fvar<fvar<var> >>, std::vector<fvar<fvar<var> >>, empty, empty, empty> type_ffv_real_real_real_73;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<fvar<fvar<var> >>, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_74;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, double, empty, empty, empty> type_ffv_real_real_real_75;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, std::vector<double>, empty, empty, empty> type_ffv_real_real_real_76;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_77;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, fvar<fvar<var> >, empty, empty, empty> type_ffv_real_real_real_78;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, std::vector<fvar<fvar<var> >>, empty, empty, empty> type_ffv_real_real_real_79;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_80;
typedef std::tuple<fvar<fvar<var> >, double, double, empty, empty, empty> type_ffv_real_real_real_81;
typedef std::tuple<fvar<fvar<var> >, double, std::vector<double>, empty, empty, empty> type_ffv_real_real_real_82;
typedef std::tuple<fvar<fvar<var> >, double, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_83;
typedef std::tuple<fvar<fvar<var> >, double, fvar<fvar<var> >, empty, empty, empty> type_ffv_real_real_real_84;
typedef std::tuple<fvar<fvar<var> >, double, std::vector<fvar<fvar<var> >>, empty, empty, empty> type_ffv_real_real_real_85;
typedef std::tuple<fvar<fvar<var> >, double, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_86;
typedef std::tuple<fvar<fvar<var> >, std::vector<double>, double, empty, empty, empty> type_ffv_real_real_real_87;
typedef std::tuple<fvar<fvar<var> >, std::vector<double>, std::vector<double>, empty, empty, empty> type_ffv_real_real_real_88;
typedef std::tuple<fvar<fvar<var> >, std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_89;
typedef std::tuple<fvar<fvar<var> >, std::vector<double>, fvar<fvar<var> >, empty, empty, empty> type_ffv_real_real_real_90;
typedef std::tuple<fvar<fvar<var> >, std::vector<double>, std::vector<fvar<fvar<var> >>, empty, empty, empty> type_ffv_real_real_real_91;
typedef std::tuple<fvar<fvar<var> >, std::vector<double>, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_92;
typedef std::tuple<fvar<fvar<var> >, Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty, empty, empty> type_ffv_real_real_real_93;
typedef std::tuple<fvar<fvar<var> >, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, empty, empty, empty> type_ffv_real_real_real_94;
typedef std::tuple<fvar<fvar<var> >, Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_95;
typedef std::tuple<fvar<fvar<var> >, Eigen::Matrix<double, Eigen::Dynamic, 1>, fvar<fvar<var> >, empty, empty, empty> type_ffv_real_real_real_96;
typedef std::tuple<fvar<fvar<var> >, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<fvar<fvar<var> >>, empty, empty, empty> type_ffv_real_real_real_97;
typedef std::tuple<fvar<fvar<var> >, Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_98;
typedef std::tuple<fvar<fvar<var> >, fvar<fvar<var> >, double, empty, empty, empty> type_ffv_real_real_real_99;
typedef std::tuple<fvar<fvar<var> >, fvar<fvar<var> >, std::vector<double>, empty, empty, empty> type_ffv_real_real_real_100;
typedef std::tuple<fvar<fvar<var> >, fvar<fvar<var> >, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_101;
typedef std::tuple<fvar<fvar<var> >, fvar<fvar<var> >, fvar<fvar<var> >, empty, empty, empty> type_ffv_real_real_real_102;
typedef std::tuple<fvar<fvar<var> >, fvar<fvar<var> >, std::vector<fvar<fvar<var> >>, empty, empty, empty> type_ffv_real_real_real_103;
typedef std::tuple<fvar<fvar<var> >, fvar<fvar<var> >, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_104;
typedef std::tuple<fvar<fvar<var> >, std::vector<fvar<fvar<var> >>, double, empty, empty, empty> type_ffv_real_real_real_105;
typedef std::tuple<fvar<fvar<var> >, std::vector<fvar<fvar<var> >>, std::vector<double>, empty, empty, empty> type_ffv_real_real_real_106;
typedef std::tuple<fvar<fvar<var> >, std::vector<fvar<fvar<var> >>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_107;
typedef std::tuple<fvar<fvar<var> >, std::vector<fvar<fvar<var> >>, fvar<fvar<var> >, empty, empty, empty> type_ffv_real_real_real_108;
typedef std::tuple<fvar<fvar<var> >, std::vector<fvar<fvar<var> >>, std::vector<fvar<fvar<var> >>, empty, empty, empty> type_ffv_real_real_real_109;
typedef std::tuple<fvar<fvar<var> >, std::vector<fvar<fvar<var> >>, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_110;
typedef std::tuple<fvar<fvar<var> >, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, double, empty, empty, empty> type_ffv_real_real_real_111;
typedef std::tuple<fvar<fvar<var> >, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, std::vector<double>, empty, empty, empty> type_ffv_real_real_real_112;
typedef std::tuple<fvar<fvar<var> >, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_113;
typedef std::tuple<fvar<fvar<var> >, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, fvar<fvar<var> >, empty, empty, empty> type_ffv_real_real_real_114;
typedef std::tuple<fvar<fvar<var> >, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, std::vector<fvar<fvar<var> >>, empty, empty, empty> type_ffv_real_real_real_115;
typedef std::tuple<fvar<fvar<var> >, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_116;
typedef std::tuple<std::vector<fvar<fvar<var> >>, double, double, empty, empty, empty> type_ffv_real_real_real_117;
typedef std::tuple<std::vector<fvar<fvar<var> >>, double, std::vector<double>, empty, empty, empty> type_ffv_real_real_real_118;
typedef std::tuple<std::vector<fvar<fvar<var> >>, double, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_119;
typedef std::tuple<std::vector<fvar<fvar<var> >>, double, fvar<fvar<var> >, empty, empty, empty> type_ffv_real_real_real_120;
typedef std::tuple<std::vector<fvar<fvar<var> >>, double, std::vector<fvar<fvar<var> >>, empty, empty, empty> type_ffv_real_real_real_121;
typedef std::tuple<std::vector<fvar<fvar<var> >>, double, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_122;
typedef std::tuple<std::vector<fvar<fvar<var> >>, std::vector<double>, double, empty, empty, empty> type_ffv_real_real_real_123;
typedef std::tuple<std::vector<fvar<fvar<var> >>, std::vector<double>, std::vector<double>, empty, empty, empty> type_ffv_real_real_real_124;
typedef std::tuple<std::vector<fvar<fvar<var> >>, std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_125;
typedef std::tuple<std::vector<fvar<fvar<var> >>, std::vector<double>, fvar<fvar<var> >, empty, empty, empty> type_ffv_real_real_real_126;
typedef std::tuple<std::vector<fvar<fvar<var> >>, std::vector<double>, std::vector<fvar<fvar<var> >>, empty, empty, empty> type_ffv_real_real_real_127;
typedef std::tuple<std::vector<fvar<fvar<var> >>, std::vector<double>, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_128;
typedef std::tuple<std::vector<fvar<fvar<var> >>, Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty, empty, empty> type_ffv_real_real_real_129;
typedef std::tuple<std::vector<fvar<fvar<var> >>, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, empty, empty, empty> type_ffv_real_real_real_130;
typedef std::tuple<std::vector<fvar<fvar<var> >>, Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_131;
typedef std::tuple<std::vector<fvar<fvar<var> >>, Eigen::Matrix<double, Eigen::Dynamic, 1>, fvar<fvar<var> >, empty, empty, empty> type_ffv_real_real_real_132;
typedef std::tuple<std::vector<fvar<fvar<var> >>, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<fvar<fvar<var> >>, empty, empty, empty> type_ffv_real_real_real_133;
typedef std::tuple<std::vector<fvar<fvar<var> >>, Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_134;
typedef std::tuple<std::vector<fvar<fvar<var> >>, fvar<fvar<var> >, double, empty, empty, empty> type_ffv_real_real_real_135;
typedef std::tuple<std::vector<fvar<fvar<var> >>, fvar<fvar<var> >, std::vector<double>, empty, empty, empty> type_ffv_real_real_real_136;
typedef std::tuple<std::vector<fvar<fvar<var> >>, fvar<fvar<var> >, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_137;
typedef std::tuple<std::vector<fvar<fvar<var> >>, fvar<fvar<var> >, fvar<fvar<var> >, empty, empty, empty> type_ffv_real_real_real_138;
typedef std::tuple<std::vector<fvar<fvar<var> >>, fvar<fvar<var> >, std::vector<fvar<fvar<var> >>, empty, empty, empty> type_ffv_real_real_real_139;
typedef std::tuple<std::vector<fvar<fvar<var> >>, fvar<fvar<var> >, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_140;
typedef std::tuple<std::vector<fvar<fvar<var> >>, std::vector<fvar<fvar<var> >>, double, empty, empty, empty> type_ffv_real_real_real_141;
typedef std::tuple<std::vector<fvar<fvar<var> >>, std::vector<fvar<fvar<var> >>, std::vector<double>, empty, empty, empty> type_ffv_real_real_real_142;
typedef std::tuple<std::vector<fvar<fvar<var> >>, std::vector<fvar<fvar<var> >>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_143;
typedef std::tuple<std::vector<fvar<fvar<var> >>, std::vector<fvar<fvar<var> >>, fvar<fvar<var> >, empty, empty, empty> type_ffv_real_real_real_144;
typedef std::tuple<std::vector<fvar<fvar<var> >>, std::vector<fvar<fvar<var> >>, std::vector<fvar<fvar<var> >>, empty, empty, empty> type_ffv_real_real_real_145;
typedef std::tuple<std::vector<fvar<fvar<var> >>, std::vector<fvar<fvar<var> >>, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_146;
typedef std::tuple<std::vector<fvar<fvar<var> >>, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, double, empty, empty, empty> type_ffv_real_real_real_147;
typedef std::tuple<std::vector<fvar<fvar<var> >>, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, std::vector<double>, empty, empty, empty> type_ffv_real_real_real_148;
typedef std::tuple<std::vector<fvar<fvar<var> >>, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_149;
typedef std::tuple<std::vector<fvar<fvar<var> >>, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, fvar<fvar<var> >, empty, empty, empty> type_ffv_real_real_real_150;
typedef std::tuple<std::vector<fvar<fvar<var> >>, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, std::vector<fvar<fvar<var> >>, empty, empty, empty> type_ffv_real_real_real_151;
typedef std::tuple<std::vector<fvar<fvar<var> >>, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_152;
typedef std::tuple<Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, double, double, empty, empty, empty> type_ffv_real_real_real_153;
typedef std::tuple<Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, double, std::vector<double>, empty, empty, empty> type_ffv_real_real_real_154;
typedef std::tuple<Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, double, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_155;
typedef std::tuple<Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, double, fvar<fvar<var> >, empty, empty, empty> type_ffv_real_real_real_156;
typedef std::tuple<Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, double, std::vector<fvar<fvar<var> >>, empty, empty, empty> type_ffv_real_real_real_157;
typedef std::tuple<Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, double, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_158;
typedef std::tuple<Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, std::vector<double>, double, empty, empty, empty> type_ffv_real_real_real_159;
typedef std::tuple<Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, std::vector<double>, std::vector<double>, empty, empty, empty> type_ffv_real_real_real_160;
typedef std::tuple<Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, std::vector<double>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_161;
typedef std::tuple<Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, std::vector<double>, fvar<fvar<var> >, empty, empty, empty> type_ffv_real_real_real_162;
typedef std::tuple<Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, std::vector<double>, std::vector<fvar<fvar<var> >>, empty, empty, empty> type_ffv_real_real_real_163;
typedef std::tuple<Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, std::vector<double>, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_164;
typedef std::tuple<Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, Eigen::Matrix<double, Eigen::Dynamic, 1>, double, empty, empty, empty> type_ffv_real_real_real_165;
typedef std::tuple<Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<double>, empty, empty, empty> type_ffv_real_real_real_166;
typedef std::tuple<Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_167;
typedef std::tuple<Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, Eigen::Matrix<double, Eigen::Dynamic, 1>, fvar<fvar<var> >, empty, empty, empty> type_ffv_real_real_real_168;
typedef std::tuple<Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<fvar<fvar<var> >>, empty, empty, empty> type_ffv_real_real_real_169;
typedef std::tuple<Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_170;
typedef std::tuple<Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, fvar<fvar<var> >, double, empty, empty, empty> type_ffv_real_real_real_171;
typedef std::tuple<Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, fvar<fvar<var> >, std::vector<double>, empty, empty, empty> type_ffv_real_real_real_172;
typedef std::tuple<Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, fvar<fvar<var> >, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_173;
typedef std::tuple<Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, fvar<fvar<var> >, fvar<fvar<var> >, empty, empty, empty> type_ffv_real_real_real_174;
typedef std::tuple<Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, fvar<fvar<var> >, std::vector<fvar<fvar<var> >>, empty, empty, empty> type_ffv_real_real_real_175;
typedef std::tuple<Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, fvar<fvar<var> >, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_176;
typedef std::tuple<Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, std::vector<fvar<fvar<var> >>, double, empty, empty, empty> type_ffv_real_real_real_177;
typedef std::tuple<Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, std::vector<fvar<fvar<var> >>, std::vector<double>, empty, empty, empty> type_ffv_real_real_real_178;
typedef std::tuple<Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, std::vector<fvar<fvar<var> >>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_179;
typedef std::tuple<Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, std::vector<fvar<fvar<var> >>, fvar<fvar<var> >, empty, empty, empty> type_ffv_real_real_real_180;
typedef std::tuple<Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, std::vector<fvar<fvar<var> >>, std::vector<fvar<fvar<var> >>, empty, empty, empty> type_ffv_real_real_real_181;
typedef std::tuple<Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, std::vector<fvar<fvar<var> >>, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_182;
typedef std::tuple<Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, double, empty, empty, empty> type_ffv_real_real_real_183;
typedef std::tuple<Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, std::vector<double>, empty, empty, empty> type_ffv_real_real_real_184;
typedef std::tuple<Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_185;
typedef std::tuple<Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, fvar<fvar<var> >, empty, empty, empty> type_ffv_real_real_real_186;
typedef std::tuple<Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, std::vector<fvar<fvar<var> >>, empty, empty, empty> type_ffv_real_real_real_187;
typedef std::tuple<Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, empty, empty, empty> type_ffv_real_real_real_188;

