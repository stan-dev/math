#ifndef STAN_MATH_PRIM_SCAL_FUN_LAMBERTW_HPP
#define STAN_MATH_PRIM_SCAL_FUN_LAMBERTW_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <cmath>
#include <limits>

namespace stan {
namespace math {

class lambert {
 private:
  const std::vector<double> lagrange_inv_coeffs{-1,
                                                +1,
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

  inline double quick_series_approx(const double p) {
    const double ap = std::abs(p);
    double approx_val = lagrange_inv_coeffs[0] + lagrange_inv_coeffs[1] * p;
    // We could move this loop to set the bound and have one, but bounded loops
    // Are better targets for unrolling
    if (ap < 0.01159) {
      for (int i = 2; i < 7; i++) {
        approx_val += lagrange_inv_coeffs[i] * std::pow(p, i);
      }
    } else if (ap < 0.0766) {
      for (int i = 2; i < 11; i++) {
        approx_val += lagrange_inv_coeffs[i] * std::pow(p, i);
      }
    } else {
      for (int i = 2; i < 21; i++) {
        approx_val += lagrange_inv_coeffs[i] * std::pow(p, i);
      }
    }
    return approx_val;
  }
  // Used when input value is near zero on branch 0
  const std::vector<double> lagrange_inv_small_val_coeffs{
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

  inline double small_value_branch0_approx(const double val) {
    double approx_w = 0;
    for (int i = 1; i < 18; i++, i++) {
        approx_w += lagrange_inv_small_val_coeffs[i - 1] * std::pow(val, i);
    }
    for (int i = 2; i < 18; i++, i++) {
        approx_w -= lagrange_inv_small_val_coeffs[i - 1] * std::pow(val, i);
    }
    return approx_w;
  }

  inline double schroder_method(const double w, const double y) {
    const double f0 = w - y;
    const double f1 = 1 + y;
    const double f00 = f0 * f0;
    const double f11 = f1 * f1;
    const double f0y = f0 * y;
    return w
           - 4 * f0 * (6 * f1 * (f11 + f0y) + f00 * y)
                 / (f11 * (24 * f11 + 36 * f0y)
                    + f00 * (6 * y * y + 8 * f1 * y + f0y));
  }

  inline double lambert_w0_final(const double val, const int nn, const std::vector<double>& e, const std::vector<double>& a, const std::vector<double>& b) {
    int n = nn - 1;
    int jmax = 8;
    if (val <= -0.36)
      jmax = 12;
    else if (val <= -0.3)
      jmax = 11;
    else if (n <= 0)
      jmax = 10;
    else if (n <= 1)
      jmax = 9;
    double y = val * e[n + 1];
    double w = n;
    for (int j = 0; j < jmax; ++j) {
      const double wj = w + b[j];
      const double yj = y * a[j];
      if (wj < yj) {
        w = wj;
        y = yj;
      }
    }
    return schroder_method(w, y);
  }

  inline double lambert_w1_final(const double val, const int nn, const std::vector<double>& e, const std::vector<double>& a, const std::vector<double>& b) {
  int n = nn - 1;
  int jmax = 11;
  if (n >= 8)
    jmax = 8;
  else if (n >= 3)
    jmax = 9;
  else if (n >= 2)
    jmax = 10;
  double w = -n;
  double y = val * e[n - 1];
  for (int j = 0; j < jmax; ++j) {
    const double wj = w - b[j];
    const double yj = y * a[j];
    if (wj < yj) {
      w = wj;
      y = yj;
    }
  }
  return schroder_method(w, y);
}



  double lambert_w0(const double val) {
    if (abs(val) < 0.05) {
      return small_value_branch0_approx(val);
    }
    if (val < -0.35) {
      const double p2 = 2 * (E * val + 1);
      if (p2 > 0)
        return quick_series_approx(sqrt(p2));
      if (p2 == 0)
        return -1;
      return std::numeric_limits<double>::quiet_NaN();
    }

    std::vector<double> e(66);
    std::vector<double> g(65);
    std::vector<double> a(12);
    std::vector<double> b(12);

    e[0] = E;
    e[1] = 1;
    g[0] = 0;
    a[0] = sqrt(INV_E);
    b[0] = 0.5;
    double ej = 1;
    for (int jj = 2; jj < 66; ++jj) {
      const int j = jj - 1;
      ej *= E;
      e[jj] = e[j] * INV_E;
      g[j] = j * ej;
    }
    for (int j = 0, jj = 1; jj < 12; ++jj) {
      a[jj] = sqrt(a[j]);
      b[jj] = b[j] * 0.5;
      j = jj;
    }
    int n;
    for (n = 0; n <= 2; ++n) {
      if (g[n] > val) {
        printf("Jump to line 1 happened\n");
        return lambert_w0_final(val, n, e, a, b);
      }
    }
    n = 2;
    for (int j = 1; j <= 5; ++j) {
      n *= 2;
      if (g[n] > val) {
        int nh = n / 2;
        for (int j = 1; j <= 5; ++j) {
          nh /= 2;
          if (nh <= 0)
            break;
          if (g[n - nh] > val)
            n -= nh;
        }
        return lambert_w0_final(val, n, e, a, b);
      }
    }
    return std::numeric_limits<double>::quiet_NaN();
  }


  double lambert_wm1(const double val) {
    if (val >= 0) {
      return std::numeric_limits<double>::quiet_NaN();
    }
    if (val < -0.35) {
      const double p2 = 2 * (E * val + 1);
      if (p2 > 0)
        return quick_series_approx(-sqrt(p2));
      if (p2 == 0)
        return -1;
      return std::numeric_limits<double>::quiet_NaN();
    }
    std::vector<double> e(64);
    std::vector<double> g(64);
    std::vector<double> a(12);
    std::vector<double> b(12);

    e[0] = E;
    g[0] = -INV_E;
    a[0] = sqrt(E);
    b[0] = 0.5;
    double ej = INV_E;
    for (int jj = 1; jj < 64; ++jj) {
      const int j = jj - 1;
      ej *= INV_E;
      e[jj] = e[j] * E;
      g[jj] = -(jj + 1) * ej;
    }
    for (int j = 0, jj = 1; jj < 12; ++jj) {
      a[jj] = sqrt(a[j]);
      b[jj] = b[j] * 0.5;
      j = jj;
    }

    int n = 2;
    if (g[n - 1] > val) {
      return lambert_w1_final(val, n, e, a, b);
    }
    for (int j = 1; j <= 5; ++j) {
      n *= 2;
      if (g[n - 1] > val) {
        int nh = n / 2;
        for (int j = 1; j <= 5; ++j) {
          nh /= 2;
          if (nh <= 0)
            break;
          if (g[n - nh - 1] > val)
            n -= nh;
        }
      return lambert_w1_final(val, n, e, a, b);
      }
    }
    return std::numeric_limits<double>::quiet_NaN();
  }

 public:
  double operator()(const double val, const int branch = 0) {
    if (branch == 0) {
      return lambert_w0(val);
    } else if (branch == -1) {
      return lambert_wm1(val);
    } else {
      return std::numeric_limits<double>::quiet_NaN();
    }
  }
};

lambert lambert_w;
}  // namespace math
}  // namespace stan
#endif
