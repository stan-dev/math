#ifndef STAN_MATH_PRIM_MAT_PROB_ORDERED_LOGISTIC_LOG_HPP
#define STAN_MATH_PRIM_MAT_PROB_ORDERED_LOGISTIC_LOG_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/prob/ordered_logistic_lpmf.hpp>
#include <boost/math/tools/promotion.hpp>
#include <vector>

namespace stan {
  namespace math {

    /**
     * Returns the (natural) log probability of the specified integer
     * outcome given the continuous location and specified cutpoints
     * in an ordered logistic model.
     *
     * <p>Typically the continous location
     * will be the dot product of a vector of regression coefficients
     * and a vector of predictors for the outcome.
     *
     * @deprecated use <code>ordered_logistic_lpmf</code>
     *
     * @tparam propto True if calculating up to a proportion.
     * @tparam T_y Type of y variable (should be int).
     * @tparam T_loc Location type.
     * @tparam T_cut Cut-point type.
     * @param y Outcome.
     * @param lambda Location.
     * @param c Positive increasing vector of cutpoints.
     * @return Log probability of outcome given location and
     * cutpoints.
     * @throw std::domain_error If the outcome is not between 1 and
     * the number of cutpoints plus 2; if the cutpoint vector is
     * empty; if the cutpoint vector contains a non-positive,
     * non-finite value; or if the cutpoint vector is not sorted in
     * ascending order.
     */
    template <bool propto, typename T_y, typename T_loc, typename T_cut>
    typename boost::math::tools::promote_args<T_loc, T_cut>::type
    ordered_logistic_log(const T_y& y, const T_loc& lambda,
                         const Eigen::Matrix<T_cut, Eigen::Dynamic, 1>& c) {
      return ordered_logistic_lpmf<propto, T_loc, T_cut>(y, lambda, c);
    }

    /**
     * @deprecated use <code>ordered_logistic_lpmf</code>
     */
    template <typename T_y, typename T_loc, typename T_cut>
    typename boost::math::tools::promote_args<T_loc, T_cut>::type
    ordered_logistic_log(const T_y& y, const T_loc& lambda,
                         const Eigen::Matrix<T_cut, Eigen::Dynamic, 1>& c) {
      return ordered_logistic_lpmf<T_loc, T_cut>(y, lambda, c);
    }

    /**
     * Returns the (natural) log probability of the specified array 
     * of integers given the vector of continuous locations and 
     * specified cutpoints in an ordered logistic model.
     *
     * <p>Typically the continous location
     * will be the dot product of a vector of regression coefficients
     * and a vector of predictors for the outcome.
     *
     * @tparam propto True if calculating up to a proportion.
     * @tparam T_y Type of y variable (should be std::vector<int>).
     * @tparam T_loc Location type.
     * @tparam T_cut Cut-point type.
     * @param y Array of integers
     * @param lambda Vector of continuous location variables.
     * @param c Positive increasing vector of cutpoints.
     * @return Log probability of outcome given location and
     * cutpoints.
     * @throw std::domain_error If the outcome is not between 1 and
     * the number of cutpoints plus 2; if the cutpoint vector is
     * empty; if the cutpoint vector contains a non-positive,
     * non-finite value; or if the cutpoint vector is not sorted in
     * ascending order.
     * @throw std::invalid_argument If y and lambda are different
     * lengths.
     */
    template <bool propto, typename T_y, typename T_loc, typename T_cut>
    typename boost::math::tools::promote_args<T_loc, T_cut>::type
    ordered_logistic_log(const T_y& y,
                         const Eigen::Matrix<T_loc, Eigen::Dynamic, 1>& lambda,
                         const Eigen::Matrix<T_cut, Eigen::Dynamic, 1>& c) {
      return ordered_logistic_lpmf<propto, T_loc, T_cut>(y, lambda, c);
    }

    /**
     * @deprecated use <code>ordered_logistic_lpmf</code>
     */
    template <typename T_y, typename T_loc, typename T_cut>
    typename boost::math::tools::promote_args<T_loc, T_cut>::type
    ordered_logistic_log(const T_y& y,
                         const Eigen::Matrix<T_loc, Eigen::Dynamic, 1>& lambda,
                         const Eigen::Matrix<T_cut, Eigen::Dynamic, 1>& c) {
      return ordered_logistic_lpmf<T_loc, T_cut>(y, lambda, c);
    }

    /**
     * Returns the (natural) log probability of the specified array 
     * of integers given the vector of continuous locations and 
     * array of specified cutpoints in an ordered logistic model.
     *
     * <p>Typically the continous location
     * will be the dot product of a vector of regression coefficients
     * and a vector of predictors for the outcome.
     *
     * @tparam propto True if calculating up to a proportion.
     * @tparam T_y Type of y variable (should be std::vector<int>).
     * @tparam T_loc Location type.
     * @tparam T_cut Cut-point type.
     * @param y Array of integers
     * @param lambda Vector of continuous location variables.
     * @param c array of Positive increasing vectors of cutpoints.
     * @return Log probability of outcome given location and
     * cutpoints.
     * @throw std::domain_error If the outcome is not between 1 and
     * the number of cutpoints plus 2; if the cutpoint vector is
     * empty; if the cutpoint vector contains a non-positive,
     * non-finite value; or if the cutpoint vector is not sorted in
     * ascending order.
     * @throw std::invalid_argument If y and lambda are different
     * lengths.
     */
    template <bool propto, typename T_y, typename T_loc, typename T_cut>
    typename boost::math::tools::promote_args<T_loc, T_cut>::type
    ordered_logistic_log(const T_y& y,
                         const Eigen::Matrix<T_loc, Eigen::Dynamic, 1>& lambda,
                         const std::vector<Eigen::Matrix<
                                                T_cut, Eigen::Dynamic, 1>>& c) {
      return ordered_logistic_lpmf<propto, T_loc, T_cut>(y, lambda, c);
    }

    /**
     * @deprecated use <code>ordered_logistic_lpmf</code>
     */
    template <typename T_y, typename T_loc, typename T_cut>
    typename boost::math::tools::promote_args<T_loc, T_cut>::type
    ordered_logistic_log(const T_y& y,
                         const Eigen::Matrix<T_loc, Eigen::Dynamic, 1>& lambda,
                         const std::vector<Eigen::Matrix<
                                                T_cut, Eigen::Dynamic, 1>>& c) {
      return ordered_logistic_lpmf<T_loc, T_cut>(y, lambda, c);
    }
  }
}
#endif
