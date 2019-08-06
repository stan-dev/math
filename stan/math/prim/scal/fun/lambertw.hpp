#ifndef STAN_MATH_PRIM_SCAL_FUN_LAMBERTW_HPP
#define STAN_MATH_PRIM_SCAL_FUN_LAMBERTW_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <cmath>
#include <limits>

namespace stan {
namespace math {

// Contains methods for doing a quick approximation of the Lambert W function.
template <int branch>
class lambert {
 private:
  // Pre-computed series of the quick series approximation
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
  /**
   * Quick approximation of the lambert W
   * @param p Value transfored by \code{2*(e * input_value + 1)}
   * @return The solution of the lambert W function
   */
  inline double quick_series_approx(const double p) {
    const double ap = std::abs(p);
    double approx_val = lagrange_inv_coeffs[0] + lagrange_inv_coeffs[1] * p;
    // We could move this loop to set the bound dynamically, but bounded loops
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

  /** Approximation of Lambert W on the zero branch when value is small
   *
   * @param val The value to be approximated
   * @return The result of the lambert W transform.
   */
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
  /**
   * Approximate Solution of Lambert W function through Schroder's fifth-order
   * update formula.
   * @param w a first estimate for the value of the Lambert W
   * @param y Approximation of \code{val*e^-W}
   * @return Estimate of the lambert W function's output
   */
  inline double schroder_method(const double w, const double y) {
    const double f0 = w - y;  // base function
    const double f1 = 1 + y;  // first derivative wrt w
    const double f00 = f0 * f0;
    const double f11 = f1 * f1;
    const double f0y = f0 * y;
    return w
           - 4 * f0 * (6 * f1 * (f11 + f0y) + f00 * y)
                 / (f11 * (24 * f11 + 36 * f0y)
                    + f00 * (6 * y * y + 8 * f1 * y + f0y));
  }

  /**
   * Get the heuristic number of bisection iterations
   * @param val The value going through the lambert W approx.
   * @param n
   * @param branch
   * @return The number of times the bisection method should go on for.
   */
  inline int number_of_bisections(const double val, const int n) {
    int jmax;
    if (branch == 0) {
      jmax = 8;
      if (val <= -0.36) {
        jmax = 12;
      } else if (val <= -0.3) {
        jmax = 11;
      } else if (n <= 0) {
        jmax = 10;
      } else if (n <= 1) {
        jmax = 9;
      }
    } else if (branch == -1) {
      jmax = 11;
      if (n >= 8) {
        jmax = 8;
      } else if (n >= 3) {
        jmax = 9;
      } else if (n >= 2) {
        jmax = 10;
      }
    }
    return jmax;
  };

  /**
   * Approximation method with mix of series approx and finish with schroder for
   * branch -1
   */
  inline double lambert_final(const double val, const int n,
                                 const std::vector<double>& e) {
    std::vector<double> a(12);
    std::vector<double> b(12);
    double y;
    double w;
    if (branch == 0) {
      a[0] = sqrt(INV_E);
      w = n - 1;  // here
    } else if (branch == -1) {
      a[0] = sqrt(E);
      w = -(n - 1);  // here
    }
    y = val * e[n + (2 * branch)]; // b = 0 is no shift and b = -1 is n - 2
    b[0] = 0.5;

    const int jmax = number_of_bisections(val, n - 1);
    for (int j = 0; j < jmax; ++j) {
      double wj;
      if (branch == 0) {
        wj = w + b[j];
      } else if (branch == -1) {
        wj = w - b[j];  // here
      }
      const double yj = y * a[j];
      if (wj < yj) {
        w = wj;
        y = yj;
      }
      const int jj = j + 1;
      a[jj] = sqrt(a[j]);
      b[jj] = b[j] * 0.5;
    }
    return schroder_method(w, y);
  }

  inline int find_midpoint(const double val, const std::vector<double>& g) {
    int n = 2;
    for (int j = 1; j <= 5; ++j) {
      n *= 2;
      if (g[n + branch] > val) {
        int m = n / 2;
        for (int j = 1; j <= 5; ++j) {
          m /= 2;
          if (m <= 0) {
            break;
          }
          if (g[n - m + branch] > val) {
            n -= m;
          }
        }
        break;
      }
    }
    return n;
  }

  inline double lambert_w0(const double val) {
    // For certain values we can do  a quick series approx and be fine
    if (abs(val) < 0.05) {
      return small_value_branch0_approx(val);
    } else if (val < -0.35) {  // very near limit of -INV_E
      const double p2 = 2 * (E * val + 1);
      if (p2 > 0) {
        return quick_series_approx(sqrt(p2));
      } else if (p2 == 0) {
        return -1;
      } else if (val < -INV_E) {
        return std::numeric_limits<double>::quiet_NaN();
      }
    }
    std::vector<double> e(66);
    std::vector<double> g(66);

    e[0] = E;
    e[1] = 1;
    g[0] = 0;
    double ej = 1;
    for (int jj = 2; jj < 66; ++jj) {
      const int j = jj - 1;
      ej *= E;
      e[jj] = e[j] * INV_E;
      g[j] = j * ej;  // In paper this is F_k
    }
    int n;
    for (n = 0; n <= 2; ++n) {
      if (g[n] > val) {
        printf("Jump to line 1 happened\n");
        return lambert_final(val, n, e);
      }
    }
    n = find_midpoint(val, g);
    return lambert_final(val, n, e);
  }

  inline double lambert_wm1(const double val) {
    if (val >= 0) {
      return std::numeric_limits<double>::quiet_NaN();
    } else if (val <= -INV_E) {
      return std::numeric_limits<double>::quiet_NaN();
    }
    if (val < -0.35) {
      const double p2 = 2 * (E * val + 1);
      if (p2 > 0)
        return quick_series_approx(-sqrt(p2));
      if (p2 == 0)
        return -1;
    }
    std::vector<double> e(64);
    std::vector<double> g(64);

    e[0] = E;
    g[0] = -INV_E;
    double ej = INV_E;
    for (int jj = 1; jj < 64; ++jj) {
      const int j = jj - 1;
      ej *= INV_E;
      e[jj] = e[j] * E;
      g[jj] = -(jj + 1) * ej;
    }

    if (g[1] > val) {
      return lambert_final(val, 2, e);
    }
    const int n = find_midpoint(val, g);
    return lambert_final(val, n, e);
  }

 public:
  double operator()(const double val) {
    if (branch == 0) {
      return lambert_w0(val);
    } else if (branch == -1) {
      return lambert_wm1(val);
    } else {
      return std::numeric_limits<double>::quiet_NaN();
    }
  }
};

lambert<-1> lambert_w1;
lambert<0> lambert_w0;
}  // namespace math
}  // namespace stan
#endif
