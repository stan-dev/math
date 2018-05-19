#ifndef STAN_MATH_CPLX_COMPLEX_HPP
#define STAN_MATH_CPLX_COMPLEX_HPP

#include <boost/math/tools/promotion.hpp>

#include <complex>
#include <type_traits>

// this file applies to both forward and reverse mode, but isn't under mix
// because it should be included from var.hpp and fvar.hpp

namespace stan {
namespace math {
namespace internal {

template <class T>
class complex;  // forward declaration of stan's complex

/// trait to see if the template parameter is std's or stan's complex
template <class>
struct is_cplx_helper : std::false_type {};
template <class T>
struct is_cplx_helper<std::complex<T>> : std::true_type {};
template <class T>
struct is_cplx_helper<complex<T>> : std::true_type {};
template <class T>
struct is_cplx : is_cplx_helper<std::decay_t<T>> {};
template <class T>
inline constexpr bool is_cplx_v = is_cplx<T>::value;

///trait check for forward or reverse variable (helper must be specialized)
template <class>
struct is_fr_var_helper: std::false_type {};
template <class T>
struct is_fr_var: is_fr_var_helper<std::decay_t<T>> {};
template <class T>
inline constexpr bool is_fr_var_v = is_fr_var<T>::value;

///trait check if any parameter is_fr_var
template <class,class,class = void>
struct any_fr_var: std::false_type {};
template <class T,class U>
struct any_fr_var<T,U,std::enable_if_t<is_fr_var_v<T>||
 is_fr_var_v<U>>>: std::true_type {};
template <class T,class U>
inline constexpr bool any_fr_var_v = any_fr_var<T,U>::value;

/// trait check for arithmetic types
template <class T>
struct is_arith : std::integral_constant<bool, is_fr_var_v<T> ||
 std::is_arithmetic<std::decay_t<T>>::value> {};
template <class T>
inline constexpr bool is_arith_v = is_arith<T>::value;

/// trait to see if the template parameter is complex or arithmetic
template <class T>
struct is_cplx_or_arith :
 std::integral_constant<bool, is_cplx_v<T> || is_arith_v<T>> {};
template <class T>
inline constexpr bool is_cplx_or_arith_v = is_cplx_or_arith<T>::value;

/// trait to remove the complex wrapper around a type
template <class T>
struct rm_cplx_helper {
 typedef T type;
};
template <class T>
struct rm_cplx_helper<std::complex<T>> {
 typedef T type;
};
template <class T>
struct rm_cplx_helper<complex<T>> {
 typedef T type;
};
template <class T>
struct rm_cplx: rm_cplx_helper<std::decay_t<T>> {};
template <class T>
using rm_cplx_t = typename rm_cplx<T>::type;

/// trait to enforce var and related return types, not zvars
template <class T>
struct to_arith_helper {//helper must be specialized by zvars
 typedef T type;
};
template <class T,
 std::enable_if_t<is_arith_v<std::decay_t<T>>>* =nullptr>
struct to_arith {
 typedef typename to_arith_helper<std::decay_t<T>>::type type;
};
template <class T>
using to_arith_t = typename to_arith<T>::type;

/// trait to enforce std::complex<var> and related return types
template <class T>
struct to_cplx {
 typedef std::enable_if_t<is_cplx_v<T>,
  std::complex<to_arith_t<rm_cplx_t<T>>>> type;
};
template <class T>
using to_cplx_t = typename to_cplx<T>::type;

/// This class exists purely to forward the interface of std::complex
template <class T>
struct complex : std::complex<T> {
 using std::complex<T>::complex;  ///< inherit all std::complex ctors
 complex(const std::complex<T>& t)
     : std::complex<T>(t){};  ///< allow promotion
 /**without this converting ctor, copy initialization of std::complex<z_(f)var>
  * from an (f)var fails, because the default implementation in the
  * std::complex template requires the argument type to match the template
  * type. This is the same reason why std::complex<double> can't be multiplied
  * by an int.
  * tparam R type of real argument. Won't fire on non-arithmetic types.
  * tparam I type of imaginary argument.
  * param[in] real, the pure real component of the complex number.
  * param[in] imag, the pure imaginary component of the complex number*/
 template <class R = T, class I = T,
           std::enable_if_t<is_arith_v<R>>* = nullptr>
 complex(const R real = 0, const I imag = 0) : std::complex<T>(real, imag) {}

};

///complex promotion when std::complex<double> is combined with a var
//always set to fire because it also does decay conversions
//could be specialized later to short circuit 
template<class T, class U, class AT=to_arith_t<T>, class AU=to_arith_t<U>,
 class AP=typename boost::math::tools::promote_args<AT,AU>::type>
auto cplx_promote(std::complex<T>const&t){
 return std::complex<AP>(t.real(),t.imag());
}

}  // namespace internal
}  // namespace math
}  // namespace stan


namespace std {

// The free functions that follow extend std::complex. This is necessary because
// the std::complex template from which stan inherits does not allow simple
// operations such as std::complex<double>()*2, because 2 is an int, and the
// functions that handle normally std::complex only work if the  template type T
// (double in this case) appears in both function arguments. The compiler doesn't
// promote the int to double because the template type isn't found in time for
// that to happen.

// Since stan has vars and similar classes, it would also suffer the same problem
// when the other argument is anything except a var without extending the binary
// functions as we have below.

// Additionally, most functions receive and deliver std::complex
// rather than stan's complex because templates (which match names exactly) and
// some namespace ADL features in other packages like Eigen depend on it.

template <class T, class U, std::enable_if_t<
 stan::math::internal::any_fr_var_v<T,U> &&
 stan::math::internal::is_arith_v<U>>* =nullptr>
inline auto operator+(std::complex<T>const&t, U const&u) {
 auto r(stan::math::internal::cplx_promote<T,U>(t));
 r+=u;
 return r;
}

template <class T, class U,std::enable_if_t<
 stan::math::internal::any_fr_var_v<T,U> &&
 stan::math::internal::is_arith_v<U>>* =nullptr>
inline auto operator+(U const&u, std::complex<T>const&t) {
 auto r(stan::math::internal::cplx_promote<T,U>(t));
 r+=u;
 return r;
}

template <class T, class U, std::enable_if_t<
 stan::math::internal::any_fr_var_v<T,U> &&
 stan::math::internal::is_arith_v<U>>* =nullptr>
inline auto operator*(std::complex<T>const&t, U const&u) {
 auto r(stan::math::internal::cplx_promote<T,U>(t));
 r*=u;
 return r;
}

template <class T, class U, std::enable_if_t<
 stan::math::internal::any_fr_var_v<T,U> && 
 stan::math::internal::is_arith_v<U>>* =nullptr>
inline auto operator*(U const&u, std::complex<T>const&t) {
 auto r(stan::math::internal::cplx_promote<T,U>(t));
 r*=u;
 return r;
}

template <class T, class U,std::enable_if_t<
 stan::math::internal::any_fr_var_v<T,U> &&
 stan::math::internal::is_arith_v<U>>* =nullptr>
inline auto operator-(std::complex<T>const&t, U const&u) {
 auto r(stan::math::internal::cplx_promote<T,U>(t));
 r-=u;
 return r;
}

template <class T, class U,std::enable_if_t<
 stan::math::internal::any_fr_var_v<T,U> &&
 stan::math::internal::is_arith_v<U>>* =nullptr>
inline auto operator-(U const&u, std::complex<T>const&t) {
 auto r(stan::math::internal::cplx_promote<T,U>(t));
 r-=u;
 return -r;
}

template <class T, class U, std::enable_if_t<
 stan::math::internal::any_fr_var_v<T,U> &&
 stan::math::internal::is_arith_v<U>>* =nullptr>
inline auto operator/(std::complex<T>const&t, U const&u) {
 auto r(stan::math::internal::cplx_promote<T,U>(t));
 r/=u;
 return r;
}

template <class T, class U, std::enable_if_t<
 stan::math::internal::any_fr_var_v<T,U> &&
 stan::math::internal::is_arith_v<U>>* =nullptr>
inline auto operator/(U const&u, std::complex<T>const&t) {
 auto r(stan::math::internal::cplx_promote<T,U>(t));
 r/=u;
 return decltype(r)(1.0)/r;
}

/*template <class T, class U>
inline bool
operator==(std::complex<T>const&t, std::complex<U>const&u) {
  return t.real() == u.real() && t.imag() == u.imag();
}

template <class T, class U>
inline bool
operator==(std::complex<T>const&t, std::complex<U>const&u) {
  return complex<T>(u) == t;
}

template <class T, class U>
inline bool
operator==(std::complex<U>const&u, std::complex<T>const&t) {
  return complex<T>(u) == t;
}

template <class T, class U,
 std::enable_if_t<is_cplx_or_arith_v<U>>* =nullptr>
inline bool
operator==(std::complex<T>const&t, U u) {
  return t.real() == u && t.imag() == T();
}

template <class T, class U,
 std::enable_if_t<is_arith_v<U>>* =nullptr>
inline bool
operator==(U u, std::complex<T>const&t) {
  return t.real() == u && t.imag() == T();
}

template <class T, class U,
 std::enable_if_t<is_cplx_or_arith_v<U>>* =nullptr>
inline bool
operator!=(std::complex<T>const&t, U u) {
  return !(t == u);
}

template <class T, class U,
 std::enable_if_t<is_arith_v<U>>* =nullptr>
inline bool
operator!=(U u, std::complex<T>const&t) {
  return !(t == u);
}

template <class T, class U>
inline to_cplx_t<complex<T>>
pow(std::complex<T> t, std::complex<U> u) {
  using std::exp;
  using std::log;
  return t == T() ? T() : exp(u * log(t));
}

template <class T>
inline to_arith_t<T> abs(complex<T>const&t){
 using std::abs;
 return abs(std::complex<T>(t));
}

template <class T>
inline std::complex<T> sqrt(complex<T>const&t){
 using std::sqrt;
 return sqrt(std::complex<T>(t));
}*/

} // namespace std

namespace Eigen {

/// Eigen scalar op traits specialization for complex variables.
template<class T1,class T2,template<class,class>class OP>
struct ScalarBinaryOpTraits<T1,std::enable_if_t<
 stan::math::internal::is_cplx_or_arith_v<T1> && //!is_eigen VectorBlock
 stan::math::internal::is_cplx_or_arith_v<T2> && //!is_eigen VectorBlock
 !std::is_same_v<T1,T2> && //avoid Eigen's template
 ((stan::math::internal::is_cplx_v<T1>&& //next boolean avoids Eigen
    !std::is_same_v<stan::math::internal::rm_cplx_t<T1>,T2>) ||
  (stan::math::internal::is_cplx_v<T2>&& //next boolean avoids Eigen
    !std::is_same_v<T1,stan::math::internal::rm_cplx_t<T2>>)),T2>,
 OP<T1,T2>>{
  typedef std::complex<stan::math::internal::to_arith_t<
   typename boost::math::tools::promote_args<
    stan::math::internal::rm_cplx_t<std::decay_t<T1>>,
    stan::math::internal::rm_cplx_t<std::decay_t<T2>>>::type>>
   ReturnType;
};

}  // namespace Eigen

#endif
