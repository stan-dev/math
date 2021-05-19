#ifndef DOUBLE_D_HPP
#define DOUBLE_D_HPP

#include <string>
#include <cmath>

#define STRINGIFY2(...) \
#__VA_ARGS__;         \
  __VA_ARGS__
namespace stan {
namespace math {
namespace internal {
// \cond
const char* double_d_src = STRINGIFY2(
    // \endcond
    typedef struct double_d {
      double high, low;

#ifdef __cplusplus
      double_d() {}
      double_d(double a) : high(a), low(0) {}
      double_d(double h, double l) : high(h), low(l) {}
#endif
    } double_d;

    inline double_d add_quick_d_d(double a, double b) {
      double_d res;
      res.high = a + b;
      res.low = b - (res.high - a);
      return res;
    }

    inline double_d add_d_d(double a, double b) {
#ifdef __cplusplus
      using std::isfinite;
#endif
      double_d res;
      res.high = a + b;
      if (isfinite(res.high)) {
        double v = res.high - a;
        res.low = (a - (res.high - v)) + (b - v);
      } else {
        res.low = 0;
      }
      return res;
    }

    inline double_d mul_d_d(double a, double b) {
      double_d res;
      res.high = a * b;
      res.low = fma(a, b, -res.high);
      return res;
    }

    inline double_d neg(double_d a) {
      double_d res;
      res.high = -a.high;
      res.low = -a.low;
      return res;
    }

    inline double_d add_dd_dd(double_d a, double_d b) {
#ifdef __cplusplus
      using std::isfinite;
#endif
      //      double_d high = add_d_d(a.high, b.high);
      //      if (isfinite(high.high)) {
      //        double_d low = add_d_d(a.low, b.low);
      //        double_d mid = add_d_d(high.low, low.high);
      //        mid.low += low.high;
      //        double_d high2 = add_d_d(high.high, mid.high);
      //        double_d mid2 = add_d_d(mid.low, high2.high);
      //        double_d low2 = add_d_d(high2.low, mid2.low);
      //        return { mid2.high, low2.high }

      //        //        double_d v = add_quick_d_d(high.high, high.low +
      //        low.high);
      //        //        return add_quick_d_d(high.high, low.low + v.low);
      //      } else {
      //        return high;
      //      }
      double_d res;
      double r = a.high + b.high;
      if (isfinite(r)) {
        double s;
        if (fabs(a.high) > fabs(b.high)) {
          s = a.high - r + b.high;
        } else {
          s = b.high - r + a.high;
        }
        s += a.low + b.low;
        res.high = r + s;
        res.low = r - res.high + s;
      } else {
        res.high = r;
        res.low = 0;
      }
      return res;
    }

    inline double_d add_dd_d(double_d a, double b) {
#ifdef __cplusplus
      using std::isfinite;
#endif
      //      double_d high = add_d_d(a.high, b);
      //      if (isfinite(high.high)) {
      //        double_d low = add_d_d(a.low, b);
      //        double_d high2 = add_d_d(high.high, low.high);
      //        double_d mid2 = add_d_d(low.low, high2.high);
      //        double_d low2 = add_d_d(high2.low, mid2.low);
      //        return {mid2.high, low2.high};
      //      } else {
      //        return high;
      //      }

      //      //    double_d v = add_quick_d_d(high.high, high.low + a.low);
      //      //    return add_quick_d_d(high.high, v.low);
      double_d res;
      double r = a.high + b;
      if (isfinite(r)) {
        double s;
        if (fabs(a.high) > fabs(b)) {
          s = a.high - r + b;
        } else {
          s = b - r + a.high;
        }
        s += a.low;
        res.high = r + s;
        res.low = r - res.high + s;
      } else {
        res.high = r;
        res.low = 0;
      }
      return res;
    }

    inline double_d sub_dd_dd(double_d a, double_d b) {
      return add_dd_dd(a, neg(b));
    }

    inline double_d sub_dd_d(double_d a, double b) { return add_dd_d(a, -b); }

    inline double_d sub_d_dd(double a, double_d b) {
      return add_dd_d(neg(b), a);
    }

    inline double_d mul_dd_dd(double_d a, double_d b) {
#ifdef __cplusplus
      using std::isfinite;
#endif
      double_d high = mul_d_d(a.high, b.high);
      double_d res;
      if (isfinite(high.high)) {
        high.low += a.high * b.low + a.low * b.high;
        res.high = high.high + high.low;
        res.low = high.high - res.high + high.low;
      } else {
        high.low = 0;
        return high;
      }
      return res;
    }

    inline double_d mul_dd_d(double_d a, double b) {
#ifdef __cplusplus
      using std::isfinite;
#endif
      double_d high = mul_d_d(a.high, b);
      double_d res;
      if (isfinite(high.high)) {
        high.low += a.low * b;
        res.high = high.high + high.low;
        res.low = high.high - res.high + high.low;
        return res;
      } else {
        high.low = 0;
        return high;
      }
    }

    inline double_d div_dd_dd(double_d a, double_d b) {
#ifdef __cplusplus
      using std::isfinite;
#endif
      double high = a.high / b.high;
      double_d res;
      if (isfinite(high)) {
        double_d u = mul_d_d(high, b.high);
        double low = (a.high - u.high - u.low + a.low - high * b.low) / b.high;
        res.high = high + low;
        res.low = high - res.high + low;
      } else {
        res.high = high;
        res.low = 0;
      }
      return res;
    }

    inline double_d div_dd_d(double_d a, double b) {
#ifdef __cplusplus
      using std::isfinite;
#endif
      double high = a.high / b;
      double_d res;
      if (isfinite(high)) {
        double_d u = mul_d_d(high, b);
        double low = (a.high - u.high - u.low + a.low) / b;
        res.high = high + low;
        res.low = high - res.high + low;
      } else {
        res.high = high;
        res.low = 0;
      }
      return res;
    }

    inline double_d div_d_dd(double a, double_d b) {
#ifdef __cplusplus
      using std::isfinite;
#endif
      double high = a / b.high;
      double_d res;
      if (isfinite(high)) {
        double_d u = mul_d_d(high, b.high);
        double low = (a - u.high - u.low - high * b.low) / b.high;
        res.high = high + low;
        res.low = high - res.high + low;
      } else {
        res.high = high;
        res.low = 0;
      }
      return res;
    }

    inline double copysign(double a, double_d b) {
#ifdef __cplusplus
      using std::copysign;
#endif
      return copysign(a, b.high);
    }

    inline double_d fabs(double_d a) {
#ifdef __cplusplus
      using std::fabs;
#endif
      return a.high > 0 ? a : neg(a);
    }

    inline bool isnan(double_d a) {
#ifdef __cplusplus
      using std::isnan;
#endif
      return isnan(a.high);
    }

    inline bool isinf(double_d a) {
#ifdef __cplusplus
      using std::isinf;
#endif
      return isinf(a.high);
    }

    inline bool lt_dd_dd(double_d a, double_d b) {
      return a.high < b.high || (a.high == b.high && a.low < b.low);
    }

    inline bool lt_d_dd(double a, double_d b) {
      return lt_dd_dd({a, 0}, b);
    }

    inline bool lt_dd_d(double_d a, double b) {
      return lt_dd_dd(a, {b, 0});
    }

    inline bool gt_dd_dd(double_d a, double_d b) {
      return a.high > b.high || (a.high == b.high && a.low > b.low);
    }

    inline bool gt_d_dd(double a, double_d b) {
      return gt_dd_dd({a, 0}, b);
    }

    inline bool gt_dd_d(double_d a, double b) {
      return gt_dd_dd(a, {b, 0});
    }

    inline bool le_dd_dd(double_d a, double_d b) {
      return !gt_dd_dd(a, b);
    } inline bool le_d_dd(double a, double_d b) {
      return !gt_d_dd(a, b);
    } inline bool le_dd_d(double_d a, double b) { return !gt_dd_d(a, b); }

    inline bool ge_dd_dd(double_d a, double_d b) {
      return !lt_dd_dd(a, b);
    } inline bool ge_d_dd(double a, double_d b) {
      return !lt_d_dd(a, b);
    } inline bool ge_dd_d(double_d a, double b) { return !lt_dd_d(a, b); }

    // \cond
);
// \endcond

inline double_d operator+(double_d a, double b) { return add_dd_d(a, b); }
inline double_d operator+(double a, double_d b) { return add_dd_d(b, a); }
inline double_d operator+(double_d a, double_d b) { return add_dd_dd(a, b); }

inline double_d operator-(double a, double_d b) { return sub_d_dd(a, b); }
inline double_d operator-(double_d a, double b) { return sub_dd_d(a, b); }
inline double_d operator-(double_d a, double_d b) { return sub_dd_dd(a, b); }

inline double_d operator*(double_d a, double b) { return mul_dd_d(a, b); }
inline double_d operator*(double a, double_d b) { return mul_dd_d(b, a); }
inline double_d operator*(double_d a, double_d b) { return mul_dd_dd(a, b); }

inline double_d operator/(double_d a, double b) { return div_dd_d(a, b); }
inline double_d operator/(double a, double_d b) { return div_dd_d(b, a); }
inline double_d operator/(double_d a, double_d b) { return div_dd_dd(a, b); }

inline double_d operator-(double_d a) { return neg(a); }

inline bool operator<(double_d a, double_d b) { return lt_dd_dd(a, b); }
inline bool operator<(double a, double_d b) { return lt_d_dd(a, b); }
inline bool operator<(double_d a, double b) { return lt_dd_d(a, b); }
inline bool operator>(double_d a, double_d b) { return gt_dd_dd(a, b); }
inline bool operator>(double a, double_d b) { return gt_d_dd(a, b); }
inline bool operator>(double_d a, double b) { return gt_dd_d(a, b); }
inline bool operator<=(double_d a, double_d b) { return le_dd_dd(a, b); }
inline bool operator<=(double a, double_d b) { return le_d_dd(a, b); }
inline bool operator<=(double_d a, double b) { return le_dd_d(a, b); }
inline bool operator>=(double_d a, double_d b) { return ge_dd_dd(a, b); }
inline bool operator>=(double a, double_d b) { return ge_d_dd(a, b); }
inline bool operator>=(double_d a, double b) { return ge_dd_d(a, b); }

}  // namespace internal
}  // namespace math
}  // namespace stan

#endif  // DOUBLE_D_HPP
