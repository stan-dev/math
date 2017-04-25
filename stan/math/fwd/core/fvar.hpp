#ifndef STAN_MATH_FWD_CORE_FVAR_HPP
#define STAN_MATH_FWD_CORE_FVAR_HPP

#include <stan/math/prim/scal/meta/likely.hpp>
#include <stan/math/prim/scal/fun/is_nan.hpp>
#include <stan/math/fwd/scal/meta/ad_promotable.hpp>
#include <boost/utility/enable_if.hpp>
#include <ostream>

namespace stan {
  namespace math {

    template <typename T>
    struct fvar {
      T val_;  // value
      T d_;    // tangent (aka derivative)

      T val() const { return val_; }
      T tangent() const { return d_; }

      typedef fvar value_type;

      /**
       * Construct a forward variable with zero value and zero
       * tangent.
       */
      fvar() : val_(0), d_(0) { }

      /**
       * Copy construct a forward variable by copying its value and
       * tangent.
       *
       * @param x forward variable to copy
       */
      fvar(const fvar<T>& x) : val_(x.val_), d_(x.d_) {  }

      /**
       * Construct a forward variable with specified value and zero
       * tangent.  If the value is not-a-number, tangent will be set to
       * not-a-number.
       *
       * @param v value of constructed forward variable
       */
      fvar(double v) : val_(v), d_(0) {  // NOLINT
        if (unlikely(is_nan(v)))
          d_ = v;
      }

      /**
       * Construct a forward variable with the specified value and
       * tangent value.  If the value is not-a-number, tangent will be
       * set to not-a-number.
       *
       * @tparam V type of value, must be assignable to T
       * @param v value
       * @param d tangent, defaulting to 0
       */
      template <typename V>
      fvar(const V& v, double d = 0) : val_(v), d_(d) {  // NOLINT
        if (unlikely(is_nan(v)))
          d_ = v;
      }

      /**
       * Construct a forward variable with the specified value and
       * tangent value.  If the value is not-a-number, tangent will be
       * set to not-a-number.
       *
       * @tparam TV type of value, must be assignable to T
       * @tparam TD type of tangent, must be assignable to T, defaults
       * to 0
       * @param v value
       * @param d tangent
       */
      // TV and TD must be assignable to T
      template <typename TV, typename TD>
      fvar(const TV& v, const TD& d = 0.0) : val_(v), d_(d) {  // NOLINT
          d_ = val;
      }

      inline fvar<T>& operator+=(const fvar<T>& x2) {
        val_ += x2.val_;
        d_ += x2.d_;
        return *this;
      }

      inline fvar<T>& operator+=(double x2) {
        val_ += x2;
        return *this;
      }

      inline fvar<T>& operator-=(const fvar<T>& x2) {
        val_ -= x2.val_;
        d_ -= x2.d_;
        return *this;
      }

      inline fvar<T>& operator-=(double x2) {
        val_ -= x2;
        return *this;
      }

      inline fvar<T>& operator*=(const fvar<T>& x2) {
        d_ = d_ * x2.val_ + val_ * x2.d_;
        val_ *= x2.val_;
        return *this;
      }

      inline fvar<T>& operator*=(double x2) {
        val_ *= x2;
        d_ *= x2;
        return *this;
      }

      inline fvar<T>& operator/=(const fvar<T>& x2) {
        d_ = (d_ * x2.val_ - val_ * x2.d_) / (x2.val_ * x2.val_);
        val_ /= x2.val_;
        return *this;
      }

      inline fvar<T>& operator/=(double x2) {
        val_ /= x2;
        d_ /= x2;
        return *this;
      }

      inline fvar<T>& operator++() {
        ++val_;
        return *this;
      }

      inline fvar<T> operator++(int /*dummy*/) {
        fvar<T> result(val_, d_);
        ++val_;
        return result;
      }

      inline fvar<T>& operator--() {
        --val_;
        return *this;
      }

      inline fvar<T> operator--(int /*dummy*/) {
        fvar<T> result(val_, d_);
        --val_;
        return result;
      }

      friend std::ostream& operator<<(std::ostream& os, const fvar<T>& v) {
        return os << v.val_;
      }
    };
  }
}
#endif
