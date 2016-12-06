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

      fvar() : val_(0.0), d_(0.0) { }

      fvar(const fvar<T>& x) : val_(x.val_), d_(x.d_) {  }

      fvar(const T& v) : val_(v), d_(0.0) {  // NOLINT(runtime/explicit)
        if (is_nan(v))
          d_ = v;
      }

      template <typename V>
      fvar(const V& v,
           typename boost::enable_if_c<ad_promotable<V, T>::value>::type*
           dummy = 0)
        : val_(v), d_(0.0) {
        if (is_nan(v))
          d_ = v;
      }

      // TV and TD must be assignable to T
      template <typename TV, typename TD>
      fvar(const TV& val, const TD& deriv) : val_(val), d_(deriv) {
        if (unlikely(is_nan(val)))
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
