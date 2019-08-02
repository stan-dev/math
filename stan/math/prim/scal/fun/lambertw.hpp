#ifndef STAN_MATH_PRIM_SCAL_FUN_LAMBERTW_HPP
#define STAN_MATH_PRIM_SCAL_FUN_LAMBERTW_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <cmath>
#include <limits>

namespace stan {
namespace math {

template <int branch>
class lambert {
private:
  const std::vector<double> lambert_small_approx_consts {
      -1,
      1,
      -0.333333333333333333,
      +0.152777777777777778,
      -0.0796296296296296296,
      +0.0445023148148148148,
      -0.0259847148736037625,
      +0.0156356325323339212,
      -0.00961689202429943171,
      +0.00601454325295611786,
      -0.00381129803489199923,
      +0.00244087799114398267,
      -0.00157693034468678425,
      +0.00102626332050760715,
      -0.000672061631156136204,
      +0.000442473061814620910,
      -0.000292677224729627445,
      +0.000194387276054539318,
      -0.000129574266852748819,
      +0.0000866503580520812717,
      -0.0000581136075044138168};

  const std::vector<double> lambert_large_approx_consts {
  1,
  1,
  1.5,
  2.6666666666666666667,
  5.2083333333333333333,
  10.8,
  23.343055555555555556,
  52.012698412698412698,
  118.62522321428571429,
  275.57319223985890653,
  649.78717234347442681,
  1551.1605194805194805,
  3741.4497029592385495,
  9104.5002411580189358,
  22324.308512706601434,
  55103.621972903835338,
  136808.86090394293563};

  inline double lambertw_series(const double val) {
    const double abs_val = std::abs(val);
    double series_sum = 0.0;
    auto approx_iters = 20;
    if (abs_val < 0.01159) {
      approx_iters = 7;
    } else if (abs_val < 0.0766) {
      approx_iters = 11;
    }
    for (auto i = 0; i < approx_iters; i++) {
      series_sum += lambert_small_approx_consts[i] * pow(val, i);
    }
    return series_sum;
  }

  inline double lambert_w0_zero_series(const double val) {
    double series_sum = 0.0;
    for (auto i = 0; i < lambert_large_approx_consts.size(); i += 2) {
        series_sum += lambert_large_approx_consts[i] * pow(val, (i+1));
    }
    for (auto i = 1; i < lambert_large_approx_consts.size(); i += 2) {
        series_sum -= lambert_large_approx_consts[i] * pow(val, (i + 1));
    }
    return series_sum;
  }

  inline double final_result(const double w, const double y) {
    auto top = -102 * pow(w, 2) * y - 3 * pow(w, 2) + 21 * w * pow(y, 3);
    top -= 6 * w * pow(y, 2) - 52 * w * y - 3 * w + pow(y, 4) - pow(y, 3) - pow(y, 2) - 6 * y;
    auto bottom = 34 * w * y + w - 7 * pow(y, 3) + 2 * pow(y, 2) + 17 * y + 1;
    return top/bottom;
  }

  inline double lambertw0_final_approx(const std::vector<double> a_base,
    const std::vector<double> b_base, const std::vector<double> e_base,
    const double val, const double n) {
    int jmax = 8;
    if (val <= -0.36) {
      jmax = 12;
    } else if (val <= -0.3) {
      jmax = 11;
    } else if (n < 0) {
      jmax = 10;
    } else if (n < 1) {
      jmax = 9;
    }
    double y = val * e_base[n];
    double w = n - 1;
    for (int j = 0; j < jmax; ++j) {
      const double wj = w + b_base[j];
      const double yj = y * a_base[j];
      if (wj < yj) {
        w = wj;
        y = yj;
      }
    }
    return final_result(w, y);
  }

  inline double lambertw1_final_approx(const std::vector<double> a_base,
    const std::vector<double> b_base, const std::vector<double> e_base,
    const double val, const double n) {
      int jmax = 11;
      if (n > 8) {
        jmax = 8;
      } else if (n > 3) {
        jmax = 9;
      } else if (n > 2) {
        jmax = 10;
      }
      double w = -(n - 1);
      double y = val * e_base[n - 2];
      for (int j = 0; j < jmax; ++j) {
        const double wj = w - b_base[j];
        const double yj = y * a_base[j];
        if (wj < yj) {
          w = wj;
          y = yj;
        }
    }
    return final_result(w, y);
  }
  inline double lambert_w0(const double val) {
    std::vector<double> e_base(66);
    std::vector<double> g_base(65);
    std::vector<double> a_base(12);
    std::vector<double> b_base(12);
    double e_iter = 1;
    const double e_inv = 1 / E;
    e_base[0] = E;
    e_base[1] = 1;
    g_base[0] = 0;
    a_base[0] = sqrt(e_inv);
    b_base[0] = 0.5;
    for (auto i = 2; i < e_base.size(); i++) {
      e_iter *= E;
      e_base[i] = e_base[i - 1] * e_inv;
      g_base[i - 1] = (i - 1) * e_iter;
    }
    for (auto i = 1; i < a_base.size(); i++) {
      a_base[i] = sqrt(a_base[i - 1]);
      b_base[i] = b_base[i - 1] * 0.5;
    }
    if (abs(val) < 0.05) {
      return lambert_w0_zero_series(val);
    } else if (val < -0.35) {
      const double quick_approx = 2 * (E * val + 1);
      if (quick_approx > 0) {
        return lambertw_series(sqrt(quick_approx));
      } else if (quick_approx == 0) {
        return -1.0;
      }
    }
    for (auto i = 0; i <= 2; i++) {
      if (g_base[i] > val) {
        return lambertw0_final_approx(a_base, b_base, e_base, val, i);
      }
    }
    auto n = 2;
    for (int j = 1; j <= 5; ++j) {
      n *= 2;
      if (g_base[n] > val) {
        int nh = n / 2;
        for (int j = 1; j <= 5; ++j) {
          nh /= 2;
          if (nh <= 0)
            break;
          if (g_base[n-nh] > val)
            n -= nh;
        }
        return lambertw0_final_approx(a_base, b_base, e_base, val, n);
      }
    }
    // This only happens if val is too large
                  printf("%s", "Gets Herew00");
    return std::numeric_limits<double>::quiet_NaN();
  }

  double lambert_w1(const double val) {
    if (val >= 0) {
              printf("%s", "Gets Here0");
      return std::numeric_limits<double>::quiet_NaN();
    }

    std::vector<double> e_base(64);
    std::vector<double> g_base(64);
    std::vector<double> a_base(12);
    std::vector<double> b_base(12);


    const double e_first = 1 / E;
    double e_iter = e_first;
    e_base[0] = E;
    g_base[0] = -e_first;
    a_base[0] = sqrt(E);
    b_base[0] = 0.5;
    for (int i = 1; i < e_base.size(); i++) {
      e_iter *= e_first;
      e_base[i] = e_base[i - 1] * E;
      g_base[i] = -(i + 1) * e_iter;
    }
    for (int i = 1; i < a_base.size(); i++) {
      a_base[i] = sqrt(a_base[i - 1]);
      b_base[i] = b_base[i - 1] * 0.5;
    }
    if (val < -0.35) {
      const double quick_iter = 2 * (E * val + 1);
      if (quick_iter > 0) {
        return lambertw_series(-sqrt(quick_iter));
      } else if (quick_iter == 0) {
        return -1;
      } 
    }
    int n = 2;
    if (g_base[1] > val)
      return lambertw1_final_approx(a_base, b_base, e_base, val, n);
    for (int j = 1; j <= 5; j++) {
      n *= 2;
      if (g_base[n - 1] > val) {
        int nh = n / 2;
        for (int k = 1; k <= 5; k++) {
          nh /= 2;
          if (nh <= 0)
            break;
          if (g_base[n-nh - 1] > val)
            n -= nh;
        }
        return lambertw1_final_approx(a_base, b_base, e_base, val, n);
      }
    }
    return std::numeric_limits<double>::quiet_NaN();

  }
public:
  double operator()(const double val) {
    if (branch == 1) {
      return lambert_w1(val);
    } else if (branch == 0) {
      return lambert_w0(val);
    }
  }
};

auto lambert_w1 = lambert<1>();
auto lambert_w0 = lambert<0>();

}
}
#endif
