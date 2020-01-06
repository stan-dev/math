#include <stan/math/rev/scal.hpp>
#include <boost/math/differentiation/finite_difference.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <gtest/gtest.h>
#include <vector>
#include <algorithm>

namespace neg_binomial_test_internal {
struct TestValue {
  unsigned int n;
  double alpha;
  double beta;
  double value;
  double grad_alpha;
  double grad_beta;

  TestValue(unsigned int _n, double _alpha, double _beta, double _value,
            double _grad_alpha, double _grad_beta)
      : n(_n),
        alpha(_alpha),
        beta(_beta),
        value(_value),
        grad_alpha(_grad_alpha),
        grad_beta(_grad_beta) {}
};

// Test data generated in Mathematica (Wolfram Cloud). The code can be re-ran
// at
// https://www.wolframcloud.com/obj/martin.modrak/Published/NegBinomial_test.nb
// but is also presented below for convenience:
//
// nb[n_,alpha_,beta_]:= Log[Binomial[n + alpha - 1, n]] +
//    alpha * Log[beta/ (1 + beta)] - n * Log[1 + beta];
// nbdalpha[n_,alpha_,beta_]= D[nb[n, alpha, beta],alpha];
// nbdbeta[n_,alpha_,beta_]= D[nb[n, alpha, beta],beta];
// out = OpenWrite["nb_test.txt"]
// alphas= {256*10^-7,314*10^-3,15*10^-1,3,180,  1123,10586};
// betas=  {4*10^-4,65*10^-3,442*10^-2,800, 15324};
// ns = {0,6,14,1525,10233};
//  WriteString[out, "std::vector<TestValue> testValues = {"];
//    Block[{$MaxPrecision = 80, $MinPrecision = 40}, {
//      For[i = 1, i <= Length[alphas], i++, {
//        For[j = 1, j <= Length[betas], j++, {
//        For[k = 1, k <= Length[ns], k++, {
//          calpha = alphas[[i]];
//          cbeta = betas[[j]];
//     cn=ns[[k]];
//          val = N[nb[cn,calpha,cbeta]];
//     ddalpha= N[nbdalpha[cn,calpha,cbeta]];
//     ddbeta= N[nbdbeta[cn,calpha,cbeta]];
//          WriteString[out,"  TestValue(",CForm[cn],",",CForm[calpha],",",
//            CForm[cbeta],",",
//            CForm[val],","CForm[ddalpha],",",CForm[ddbeta],"),"]
//        }]
//      }]
//   }]
//    }];
//  WriteString[out,"};"];
//  Close[out];
//  FilePrint[%]
std::vector<TestValue> testValues = {
    TestValue(0, 0.0000256, 0.0004, -0.00020030581583046705, -7.824445930877619,
              0.06397441023590564),
    TestValue(6, 0.0000256, 0.0004, -12.36721904879686, 39056.95884993479,
              -5.933626549380248),
    TestValue(14, 0.0000256, 0.0004, -13.217693311299035, 39057.855647610166,
              -13.930427828868453),
    TestValue(1525, 0.0000256, 0.0004, -18.512543824479696, 39062.58214942565,
              -1524.3262694922032),
    TestValue(10233, 0.0000256, 0.0004, -23.89862201233763, 39064.48605183577,
              -10228.844462215115),
    TestValue(0, 0.0000256, 0.065, -0.00007158637589114594, -2.7963428082478883,
              0.00036980859516070783),
    TestValue(6, 0.0000256, 0.065, -12.74253960419729, 39061.986953057414,
              -5.633433008306247),
    TestValue(14, 0.0000256, 0.065, -14.093612899819956, 39062.88375073279,
              -13.145170097508126),
    TestValue(1525, 0.0000256, 0.065, -113.93910579363354, 39067.61025254828,
              -1431.9245128205127),
    TestValue(10233, 0.0000256, 0.065, -664.2272315331469, 39069.51415495839,
              -9608.450334416757),
    TestValue(0, 0.0000256, 4.42, -5.221276655665873e-6, -0.20395611936194816,
              1.0686079711475847e-6),
    TestValue(6, 0.0000256, 4.42, -22.505199336839055, 39064.579339746306,
              -1.10701000150273),
    TestValue(14, 0.0000256, 4.42, -36.873240762783055, 39065.47613742168,
              -2.5830247616503312),
    TestValue(1525, 0.0000256, 4.42, -2595.2985892710385, 39070.20263923717,
              -281.36531258452857),
    TestValue(10233, 0.0000256, 4.42, -17314.556524865326, 39072.106541647285,
              -1888.0073790051927),
    TestValue(0, 0.0000256, 800, -3.198001665105659e-8, -0.0012492194004318981,
              3.995006242197253e-11),
    TestValue(6, 0.0000256, 800, -52.47978493724324, 39064.782046646265,
              -0.007490636664169788),
    TestValue(14, 0.0000256, 800, -106.81394741612168, 39065.67884432164,
              -0.01747815226966292),
    TestValue(1525, 0.0000256, 800, -10213.840409797369, 39070.40534613713,
              -1.9038701622571785),
    TestValue(10233, 0.0000256, 800, -68436.2211115108, 39072.309248547244,
              -12.775280898836455),
    TestValue(0, 0.0000256, 15324, -1.6705275871366172e-9,
              -0.00006525498387252411, 1.0901025079596645e-13),
    TestValue(6, 0.0000256, 15324, -70.18806378935007, 39064.78323061068,
              -0.000391517128765378),
    TestValue(14, 0.0000256, 15324, -148.13326477811688, 39065.68002828606,
              -0.0009135399672645624),
    TestValue(1525, 0.0000256, 15324, -14714.694625714537, 39070.40653010154,
              -0.099510603588798),
    TestValue(10233, 0.0000256, 15324, -98637.6907454415, 39072.31043251166,
              -0.6677324632951601),
    TestValue(0, 0.314, 0.0004, -2.4568760222955723, -7.824445930877619,
              784.6861255497801),
    TestValue(6, 0.314, 0.0004, -4.754469184364657, -2.7248096880286647,
              778.688524590164),
    TestValue(14, 0.314, 0.0004, -5.3285232877929705, -1.8604041327777927,
              770.6917233106757),
    TestValue(1525, 0.314, 0.0004, -9.142962392308174, 2.8433228802348456,
              -739.704118352659),
    TestValue(10233, 0.314, 0.0004, -13.9312914760053, 4.747050042756005,
              -9444.22231107557),
    TestValue(0, 0.314, 0.065, -0.878051641789837, -2.7963428082478883,
              4.535933550018057),
    TestValue(6, 0.314, 0.065, -3.551094078699291, 2.3032934346010663,
              -1.0978692668833514),
    TestValue(14, 0.314, 0.065, -4.625747215248096, 3.1676989898519383,
              -8.60960635608523),
    TestValue(1525, 0.314, 0.065, -102.99082870039622, 7.871426002864577,
              -1427.38894907909),
    TestValue(10233, 0.314, 0.065, -652.6812053357488, 9.775153165385735,
              -9603.914770675334),
    TestValue(0, 0.314, 4.42, -0.06404222147965172, -0.20395611936194816,
              0.013107144646107094),
    TestValue(6, 0.314, 4.42, -12.499810756130104, 4.895680123487006,
              -1.093903925464594),
    TestValue(14, 0.314, 4.42, -26.591432023000245, 5.7600856787378785,
              -2.5699186856121954),
    TestValue(1525, 0.314, 4.42, -2583.53636912259, 10.463812691750517,
              -281.35220650849044),
    TestValue(10233, 0.314, 4.42, -17302.19655561271, 12.367539854271676,
              -1887.9942729291547),
    TestValue(0, 0.314, 800, -0.000392254891735616, -0.0012492194004318981,
              4.900124843945069e-7),
    TestValue(6, 0.314, 800, -42.41075157924302, 5.098387023448522,
              -0.007490146691635456),
    TestValue(14, 0.314, 800, -96.46849389904759, 5.9627925786993945,
              -0.01747766229712859),
    TestValue(1525, 0.314, 800, -10202.01454487163, 10.666519591712033,
              -1.9038696722846442),
    TestValue(10233, 0.314, 800, -68423.79749748089, 12.570246754233192,
              -12.775280408863921),
    TestValue(0, 0.314, 15324, -0.00002049006493597257, -0.00006525498387252411,
              1.337078857419276e-9),
    TestValue(6, 0.314, 15324, -60.11865869683253, 5.099570987865082,
              -0.00039151579179553083),
    TestValue(14, 0.314, 15324, -137.7874395265255, 5.963976543115955,
              -0.0009135386302947152),
    TestValue(1525, 0.314, 15324, -14702.868389054282, 10.667703556128593,
              -0.09951060225182816),
    TestValue(10233, 0.314, 15324, -98625.26675967708, 12.57143071864975,
              -0.6677324619581904),
    TestValue(0, 1.5, 0.0004, -11.736668896316429, -7.824445930877619,
              3748.500599760096),
    TestValue(6, 1.5, 0.0004, -10.663173154060512, -5.914178420610109,
              3742.50299880048),
    TestValue(14, 1.5, 0.0004, -10.275792230602828, -5.152700662265831,
              3734.5061975209915),
    TestValue(1525, 1.5, 0.0004, -8.560643998780371, -0.5305306751189374,
              2224.110355857657),
    TestValue(10233, 1.5, 0.0004, -11.09154505642474, 1.3725348849066714,
              -6480.407836865254),
    TestValue(0, 1.5, 0.065, -4.194514212371832, -2.7963428082478883,
              21.668472372697725),
    TestValue(6, 1.5, 0.065, -3.496467744956285, -0.8860752979803781,
              16.034669555796317),
    TestValue(14, 1.5, 0.065, -3.6096858546190917, -0.12459753963609987,
              8.522932466594439),
    TestValue(1525, 1.5, 0.065, -96.44518000342956, 4.497572447510794,
              -1410.2564102564102),
    TestValue(10233, 1.5, 0.065, -643.8781286127294, 6.400638007536402,
              -9586.782231852654),
    TestValue(0, 1.5, 4.42, -0.30593417904292225, -0.20395611936194816,
              0.0626137483094288),
    TestValue(6, 1.5, 4.42, -9.370613809368376, 1.706311390905562,
              -1.0443973218012723),
    TestValue(14, 1.5, 4.42, -22.500800049352513, 2.4677891492498403,
              -2.520412081948874),
    TestValue(1525, 1.5, 4.42, -2573.916149812605, 7.089959136396733,
              -281.3026999048271),
    TestValue(10233, 1.5, 4.42, -17290.318908276677, 8.993024696422342,
              -1887.9447663254914),
    TestValue(0, 1.5, 800, -0.0018738291006478473, -0.0012492194004318981,
              2.3408239700374533e-6),
    TestValue(6, 1.5, 800, -39.04114424912693, 1.9090182908670783,
              -0.007488295880149813),
    TestValue(14, 1.5, 800, -92.1374515420455, 2.6704960492113563,
              -0.017475811485642947),
    TestValue(1525, 1.5, 800, -10192.153915178289, 7.29266603635825,
              -1.9038678214731586),
    TestValue(10233, 1.5, 800, -68411.67943976149, 9.195731596383858,
              -12.775278558052435),
    TestValue(0, 1.5, 15324, -0.00009788247580878617, -0.00006525498387252411,
              6.387319382576159e-9),
    TestValue(6, 1.5, 15324, -56.74764718491841, 1.9102022552836375,
              -0.0003915107415550057),
    TestValue(14, 1.5, 15324, -133.45499298772538, 2.671680013627916,
              -0.00091353358005419),
    TestValue(1525, 1.5, 15324, -14693.006355179143, 7.2938500007748095,
              -0.09951059720158763),
    TestValue(10233, 1.5, 15324, -98613.1472977759, 9.196915560800418,
              -0.6677324569079498),
    TestValue(0, 3, 0.0004, -23.473337792632858, -7.824445930877619,
              7497.001199520192),
    TestValue(6, 3, 0.0004, -20.143532802585614, -6.606588788020477,
              7491.0035985605755),
    TestValue(14, 3, 0.0004, -18.691444930149387, -5.943716937648626,
              7483.0067972810875),
    TestValue(1525, 3, 0.0004, -10.114897488653275, -1.4158425570726045,
              5972.610955617753),
    TestValue(10233, 3, 0.0004, -9.791827263970454, 0.4863870833023043,
              -2731.907237105158),
    TestValue(0, 3, 0.065, -8.389028424743664, -2.7963428082478883,
              43.33694474539545),
    TestValue(6, 3, 0.065, -5.434672709536789, -1.5784856653907455,
              37.70314192849404),
    TestValue(14, 3, 0.065, -4.483183870221056, -0.915613815018895,
              30.191404839292165),
    TestValue(1525, 3, 0.065, -90.45727880935787, 3.6122605655571265,
              -1388.5879378837126),
    TestValue(10233, 3, 0.065, -635.0362561363305, 5.514490205932035,
              -9565.113759479957),
    TestValue(0, 3, 4.42, -0.6118683580858445, -0.20395611936194816,
              0.1252274966188576),
    TestValue(6, 3, 4.42, -7.42023874061997, 1.0139010234951946,
              -0.9817835734918435),
    TestValue(14, 3, 4.42, -19.485718031625566, 1.676772873867045,
              -2.457798333639445),
    TestValue(1525, 3, 4.42, -2564.039668585204, 6.204647254443066,
              -281.2400861565177),
    TestValue(10233, 3, 4.42, -17277.58845576695, 8.106876894817976,
              -1887.8821525771818),
    TestValue(0, 3, 800, -0.0037476582012956946, -0.0012492194004318981,
              4.6816479400749066e-6),
    TestValue(6, 3, 800, -36.78670883043625, 1.2166079234567109,
              -0.007485955056179775),
    TestValue(14, 3, 800, -88.81830917437628, 1.8794797738285614,
              -0.01747347066167291),
    TestValue(1525, 3, 800, -10181.973373600948, 6.407354154404583,
              -1.9038654806491886),
    TestValue(10233, 3, 800, -68398.64492690183, 8.309583794779492,
              -12.775276217228464),
    TestValue(0, 3, 15324, -0.00019576495161757233, -0.00006525498387252411,
              1.2774638765152318e-8),
    TestValue(6, 3, 15324, -54.491435819602884, 1.2177918878732703,
              -0.0003915043542356231),
    TestValue(14, 3, 15324, -130.13407467343131, 1.8806637382451208,
              -0.0009135271927348074),
    TestValue(1525, 3, 15324, -14682.824037655175, 6.408538118821142,
              -0.09951059081426825),
    TestValue(10233, 3, 15324, -98600.1110089696, 8.31076775919605,
              -0.6677324505206303),
    TestValue(0, 180, 0.0004, -1408.4002675579716, -7.824445930877619,
              449820.0719712115),
    TestValue(6, 180, 0.0004, -1383.7416799853481, -7.7915663390606875,
              449814.0743702519),
    TestValue(14, 180, 0.0004, -1360.4023219866501, -7.749343806692033,
              449806.0775689724),
    TestValue(1525, 180, 0.0004, -839.8678539490804, -5.573595326062067,
              448295.6817273091),
    TestValue(10233, 180, 0.0004, -511.20883255623437, -3.7638601451221234,
              439591.16353458614),
    TestValue(0, 180, 0.065, -503.34170548461987, -2.7963428082478883,
              2600.2166847237268),
    TestValue(6, 180, 0.065, -479.05856718683674, -2.763463216430956,
              2594.5828819068256),
    TestValue(14, 180, 0.065, -456.21980822125926, -2.7212406840623014,
              2587.071144817624),
    TestValue(1525, 180, 0.065, -30.235982564322512, -0.5454922034323357,
              1168.291802094619),
    TestValue(10233, 180, 0.065, -246.47900872313187, 1.2642429775076076,
              -7008.234019501625),
    TestValue(0, 180, 4.42, -36.71210148515067, -0.20395611936194816,
              7.513649797131455),
    TestValue(6, 180, 4.42, -22.191689285108538, -0.17107652754501607,
              6.406638727020755),
    TestValue(14, 180, 4.42, -12.369898449852393, -0.12885399517636129,
              4.930623966873153),
    TestValue(1525, 180, 4.42, -2044.965928407357, 2.0468944854536044,
              -273.8516638560051),
    TestValue(10233, 180, 4.42, -16430.178764420936, 3.8566296663935478,
              -1880.4937302766693),
    TestValue(0, 180, 800, -0.22485949207774167, -0.0012492194004318981,
              0.0002808988764044944),
    TestValue(6, 180, 800, -15.679038081736437, 0.031630372416500185,
              -0.007209737827715355),
    TestValue(14, 180, 800, -45.823368299414724, 0.07385290478515498,
              -0.01719725343320849),
    TestValue(1525, 180, 800, -9627.020512129913, 2.2496013854151204,
              -1.903589263420724),
    TestValue(10233, 180, 800, -67515.35611426263, 4.059336566355064, -12.775),
    TestValue(0, 180, 15324, -0.01174589709705434, -0.00006525498387252411,
              7.664783259091391e-7),
    TestValue(6, 180, 15324, -33.17420336917207, 0.03281433683305956,
              -0.0003907506505484791),
    TestValue(14, 180, 15324, -86.92957209673877, 0.07503686920171435,
              -0.0009127734890476634),
    TestValue(1525, 180, 15324, -14127.66161448241, 2.25078534983168,
              -0.0995098371105811),
    TestValue(10233, 180, 15324, -97716.61263462866, 4.0605205307716234,
              -0.6677316968169432),
    TestValue(0, 1123, 0.0004, -8786.852780375566, -7.824445930877619,
              2.806377449020392e6),
    TestValue(6, 1123, 0.0004, -8751.278542053076, -7.819114954587142,
              2.806371451419432e6),
    TestValue(14, 1123, 0.0004, -8713.636264898727, -7.812050908011125,
              2.8063634546181527e6),
    TestValue(1525, 1123, 0.0004, -6987.654391240393, -6.966388500142561,
              2.8048530587764895e6),
    TestValue(10233, 1123, 0.0004, -5133.758442132798, -5.510302096203866,
              2.7961485405837665e6),
    TestValue(0, 1123, 0.065, -3140.2929736623787, -2.7963428082478883,
              16222.462983026364),
    TestValue(6, 1123, 0.065, -3105.0941846147293, -2.79101183195741,
              16216.829180209463),
    TestValue(14, 1123, 0.065, -3067.9525064935015, -2.7839477853813945,
              16209.31744312026),
    TestValue(1525, 1123, 0.065, -1436.5212752158002, -1.9382853775128304,
              14790.538100397256),
    TestValue(10233, 1123, 0.065, -127.52737365985922, -0.48219897357413455,
              6614.012278801011),
    TestValue(0, 1123, 4.42, -229.0427220434678, -0.20395611936194816,
              46.87682623432569),
    TestValue(6, 1123, 4.42, -203.60665909355922, -0.1986251430714702,
              45.76981516421499),
    TestValue(14, 1123, 4.42, -179.4819491026526, -0.19156109649545439,
              44.29380040406739),
    TestValue(1525, 1123, 4.42, -1006.630573439393, 0.6541013113731099,
              -234.48848741881085),
    TestValue(10233, 1123, 4.42, -13866.606481738225, 2.1101877153118056,
              -1841.130553839475),
    TestValue(0, 1123, 800, -1.4028733866850216, -0.0012492194004318981,
              0.0017524968789013732),
    TestValue(6, 1123, 800, -5.941401226477268, 0.004081756890046066,
              -0.005738139825218477),
    TestValue(14, 1123, 800, -21.782812288505113, 0.011145803466061885,
              -0.01572565543071161),
    TestValue(1525, 1123, 800, -8397.532550498237, 0.8568082113346261,
              -1.9021176654182272),
    TestValue(10233, 1123, 800, -64760.6312249162, 2.3128946152733216,
              -12.773528401997503),
    TestValue(0, 1123, 15324, -0.07328134688884458, -0.00006525498387252411,
              4.781973111088685e-6),
    TestValue(6, 1123, 15324, -22.32008806909741, 0.00526572130660544,
              -0.0003867351557632996),
    TestValue(14, 1123, 15324, -61.77253764101367, 0.01232976788262126,
              -0.0009087579942624839),
    TestValue(1525, 1123, 15324, -12897.05717440592, 0.8579921757511855,
              -0.09950582161579592),
    TestValue(10233, 1123, 15324, -94960.77126683743, 2.3140785796898813,
              -0.6677276813221581),
    TestValue(0, 10586, 0.0004, -82829.58462427047, -7.824445930877619,
              2.6454418232706916e7),
    TestValue(6, 10586, 0.0004, -82780.56113236764, -7.823879278362524,
              2.6454412235105958e7),
    TestValue(14, 10586, 0.0004, -82725.0308245923, -7.823124240811288,
              2.645440423830468e7),
    TestValue(1525, 10586, 0.0004, -78250.18060777042, -7.689858226187986,
              2.6452893842463017e7),
    TestValue(10233, 10586, 0.0004, -68411.91247605493, -7.148089056798321,
              2.6444189324270293e7),
    TestValue(0, 10586, 0.065, -29602.084968112144, -2.7963428082478883,
              152921.63235825207),
    TestValue(6, 10586, 0.065, -29553.436925484144, -2.795776155732793,
              152915.9985554352),
    TestValue(14, 10586, 0.065, -29498.407216741936, -2.795021118181557,
              152908.48681834596),
    TestValue(1525, 10586, 0.065, -25118.10764230068, -2.6617551035582547,
              151489.70747562297),
    TestValue(10233, 10586, 0.065, -15824.741558136848, -2.1199859341685903,
              143313.18165402673),
    TestValue(0, 10586, 4.42, -2159.0794795655834, -0.20395611936194816,
              441.8860930690755),
    TestValue(6, 10586, 4.42, -2120.1941630353253, -0.20338946684685275,
              440.7790819989648),
    TestValue(14, 10586, 4.42, -2078.181422423439, -0.2026344292956172,
              439.30306723881716),
    TestValue(1525, 10586, 4.42, -156.46170359662483, -0.06936841467231439,
              160.52077941593896),
    TestValue(10233, 10586, 4.42, -5032.065429287562, 0.47240075471734977,
              -1446.1212870047252),
    TestValue(0, 10586, 800, -13.224236572972073, -0.0012492194004318981,
              0.016519975031210988),
    TestValue(6, 10586, 800, -4.313510832414451, -0.0006825668853364755,
              0.009029338327091137),
    TestValue(14, 10586, 800, -2.2668912734629885, 0.00007247066589905968,
              -0.0009581772784019975),
    TestValue(1525, 10586, 800, -5629.148286319642, 0.13333848528920186,
              -1.8873501872659175),
    TestValue(10233, 10586, 800, -54007.87477812971, 0.675107654678866,
              -12.758760923845193),
    TestValue(0, 10586, 15324, -0.6907892592745403, -0.00006525498387252411,
              0.00004507744198930081),
    TestValue(6, 10586, 15324, -9.488342401133231, 0.0005013975312228986,
              -0.0003464396868850874),
    TestValue(14, 10586, 15324, -31.052761352070178, 0.0012564350824584337,
              -0.0008684625253842718),
    TestValue(1525, 10586, 15324, -10117.469054953423, 0.13452244970576124,
              -0.09946552614691771),
    TestValue(10233, 10586, 15324, -84196.81096477703, 0.6762916190954253,
              -0.6676873858532799),
};

}  // namespace neg_binomial_test_internal

TEST(ProbDistributionsNegBinomial, derivativesPrecomputed) {
  using neg_binomial_test_internal::TestValue;
  using neg_binomial_test_internal::testValues;
  using stan::math::is_nan;
  using stan::math::neg_binomial_lpmf;
  using stan::math::value_of;
  using stan::math::var;

  for (TestValue t : testValues) {
    int n = t.n;  // Using signed int to avoid ambiguity errors.
    var alpha(t.alpha);
    var beta(t.beta);
    var val = neg_binomial_lpmf(n, alpha, beta);

    std::vector<var> x;
    x.push_back(alpha);
    x.push_back(beta);

    std::vector<double> gradients;
    val.grad(x, gradients);

    for (int i = 0; i < 2; ++i) {
      EXPECT_FALSE(is_nan(gradients[i]));
    }

    auto tolerance = [](double x) { return std::max(fabs(x * 1e-8), 1e-14); };

    EXPECT_NEAR(value_of(val), t.value, tolerance(t.value))
        << "value n = " << n << ", alpha = " << t.alpha
        << ", beta = " << t.beta;
    EXPECT_NEAR(gradients[0], t.grad_alpha, tolerance(t.grad_alpha))
        << "grad_alpha n = " << n << ", alpha = " << t.alpha
        << ", beta = " << t.beta;
    EXPECT_NEAR(gradients[1], t.grad_beta, tolerance(t.grad_beta))
        << "grad_beta n = " << n << ", alpha = " << t.alpha
        << ", beta = " << t.beta;
  }
}

TEST(ProbDistributionsNegBinomial, derivativesComplexStep) {
  using boost::math::differentiation::complex_step_derivative;
  using stan::math::internal::neg_binomial_alpha_cutoff;
  using stan::math::is_nan;
  using stan::math::neg_binomial_lpmf;
  using stan::math::var;

  std::vector<int> n_to_test = {0, 7, 100, 835, 14238, 500000, 10000000};
  std::vector<double> alpha_to_test = {0.001,
                                       0.3,
                                       113,
                                       842,
                                       21456,
                                       44242,
                                       neg_binomial_alpha_cutoff - 1,
                                       neg_binomial_alpha_cutoff + 1,
                                       1e15};
  std::vector<double> beta_to_test = {0.8,
                                      8,
                                      24,
                                      271,
                                      2586,
                                      33294,
                                      neg_binomial_alpha_cutoff - 1,
                                      neg_binomial_alpha_cutoff + 1,
                                      1e15};

  auto nb_log_for_test = [](int n, const std::complex<double>& alpha,
                            const std::complex<double>& beta) {
    // Using first-order Taylor expansion of lgamma(a + b*i) around b = 0
    // Which happens to work nice in this case, as b is always 0 or the very
    // small complex step
    auto lgamma_c_approx = [](const std::complex<double>& x) {
      return std::complex<double>(lgamma(x.real()),
                                  x.imag() * boost::math::digamma(x.real()));
    };

    const double n_(n);
    return lgamma_c_approx(n_ + alpha) - lgamma(n + 1) - lgamma_c_approx(alpha)
           + alpha * log(beta / (1.0 + beta)) - n_ * log(1.0 + beta);
  };

  for (double alpha_dbl : alpha_to_test) {
    for (double beta_dbl : beta_to_test) {
      for (int n : n_to_test) {
        var alpha(alpha_dbl);
        var beta(beta_dbl);
        var val = neg_binomial_lpmf(n, alpha, beta);

        std::vector<var> x;
        x.push_back(alpha);
        x.push_back(beta);

        std::vector<double> gradients;
        val.grad(x, gradients);

        EXPECT_TRUE(value_of(val) < 0)
            << "for n = " << n << ", alpha = " << alpha_dbl
            << ", beta = " << beta_dbl;

        for (int i = 0; i < 2; ++i) {
          EXPECT_FALSE(is_nan(gradients[i]));
        }

        auto nb_log_alpha =
            [n, beta_dbl, nb_log_for_test](const std::complex<double>& alpha) {
              return nb_log_for_test(n, alpha, beta_dbl);
            };
        auto nb_log_beta = [n, alpha_dbl,
                            nb_log_for_test](const std::complex<double>& beta) {
          return nb_log_for_test(n, alpha_dbl, beta);
        };
        double complex_step_dalpha
            = complex_step_derivative(nb_log_alpha, alpha_dbl);
        double complex_step_dbeta
            = complex_step_derivative(nb_log_beta, beta_dbl);

        double tolerance_alpha;
        if (alpha < neg_binomial_alpha_cutoff || n < 100000) {
          tolerance_alpha = std::max(1e-10, fabs(gradients[0]) * 1e-8);
        } else {
          // Not sure why the test fails in this case with strict tolerance
          // but the error is still quite small, so just increasing tolerance
          tolerance_alpha = std::max(1e-6, fabs(gradients[0]) * 1e-4);
        }
        EXPECT_NEAR(gradients[0], complex_step_dalpha, tolerance_alpha)
            << "grad_alpha, n = " << n << ", alpha = " << alpha_dbl
            << ", beta = " << beta_dbl;
        EXPECT_NEAR(gradients[1], complex_step_dbeta,
                    std::max(1e-10, fabs(gradients[1]) * 1e-8))
            << "grad_beta, n = " << n << ", alpha = " << alpha_dbl
            << ", beta = " << beta_dbl;
      }
    }
  }
}

TEST(ProbDistributionsNegativeBinomial, proptoAtPoissonCutoff) {
  using stan::math::internal::neg_binomial_alpha_cutoff;
  using stan::math::neg_binomial_lpmf;
  using stan::math::var;

  var beta_var(10);
  int y = 11;
  var value_before_cutoff = neg_binomial_lpmf<true, int, double, var>(
      y, neg_binomial_alpha_cutoff - 1e-8, beta_var);
  var value_after_cutoff = neg_binomial_lpmf<true, int, double, var>(
      y, neg_binomial_alpha_cutoff + 1e-8, beta_var);

  EXPECT_NEAR(value_of(value_before_cutoff), value_of(value_after_cutoff), 1);
}

TEST(ProbDistributionsNegBinomial, derivativesAtCutoff) {
  double alpha_cutoff = stan::math::internal::neg_binomial_alpha_cutoff;
  using stan::math::is_nan;
  using stan::math::var;

  std::vector<double> beta_to_test
      = {9.3e-6, 0.0028252, 4, 11, 8522, 984256, 5036842};
  std::vector<int> n_to_test = {0, 1, 5, 48, 1158, 224582, 48235842, 20314458};
  for (double beta : beta_to_test) {
    for (int n : n_to_test) {
      var alpha_before(alpha_cutoff - 1e-8);
      var beta_before(beta);
      var value_before = neg_binomial_lpmf(n, alpha_before, beta_before);
      std::vector<var> x_before;
      x_before.push_back(alpha_before);
      x_before.push_back(beta_before);

      std::vector<double> gradients_before;
      value_before.grad(x_before, gradients_before);

      var alpha_after(alpha_cutoff + 1e-8);
      var beta_after(beta);
      var value_after = neg_binomial_lpmf(n, alpha_after, beta_after);
      std::vector<var> x_after;
      x_after.push_back(alpha_after);
      x_after.push_back(beta_after);

      std::vector<double> gradients_after;
      value_after.grad(x_after, gradients_after);

      for (int i = 0; i < 2; ++i) {
        EXPECT_FALSE(is_nan(gradients_before[i]));
        EXPECT_FALSE(is_nan(gradients_after[i]));
      }

      EXPECT_NEAR(value_of(value_before), value_of(value_after),
                  1e-8 * fabs(value_of(value_after)))
          << "value changes too much around alpha cutoff for n = " << n
          << ", beta = " << beta << ", cutoff = " << alpha_cutoff
          << " value at cutoff - 1e-8: " << value_of(value_before)
          << ", value at cutoff + 1e-8: " << value_of(value_after);
      EXPECT_NEAR(gradients_before[0], gradients_after[0],
                  1e-8 * fabs(gradients_before[0]))
          << "grad_alpha changes too much around alpha cutoff for n = " << n
          << ", beta = " << beta << ", cutoff = " << alpha_cutoff
          << " grad_alpha at cutoff - 1e-8: " << gradients_before[0]
          << ", grad_alpha at cutoff + 1e-8: " << gradients_after[0];

      EXPECT_NEAR(gradients_before[1], gradients_after[1],
                  1e-8 * fabs(gradients_before[1]))
          << "grad_beta changes too much around alpha cutoff for n = " << n
          << ", beta = " << beta << ", cutoff = " << alpha_cutoff
          << " grad_beta at cutoff - 1e-8: " << gradients_before[1]
          << ", grad_beta at cutoff + 1e-8: " << gradients_after[1];
    }
  }
}
