#ifndef STAN_MATH_PRIM_FUN_INTERP1_HPP
#define STAN_MATH_PRIM_FUN_INTERP1_HPP

#include <stan/math/prim/meta.hpp>
using namespace std;
namespace stan {
    namespace math {
    /**
     * Return the one-dimensional interpolation of the second argument of points corresponding to the first argument.
     * @tparam T_true type of the true argument
     * @tparam T_false type of the false argument
     * @param c Boolean condition value.
     * @param y_true Value to return if condition is true.
     * @param y_false Value to return if condition is false.
     */
    inline Eigen::VectorXd
    interp1(const Eigen::VectorXd& xData, const Eigen::VectorXd& yData, const Eigen::VectorXd& xVector)
    {

        Eigen::VectorXd yVals;
        int size_x = xData.size();
        int size_xData = xData.size();

        for(int j=0; j < size_x; j++) {
            double x = xVector[j];
            int i = 0;
            if (x >= xData[size_xData - 2])                                           // special case: beyond right end
            {
                i = size_xData - 2;
            } else {
                while (x > xData[i + 1])
                    i++;                                   // find left end of interval for interpolation
            }
            double xL = xData[i], yL = yData[i], xR = xData[i + 1], yR = yData[i + 1]; // points on either side (unless beyond ends)
            if (x < xL) yR = yL;                                                 // if beyond ends of array
            if (x > xR) yL = yR;
            double dydx = (yR - yL) / (xR - xL);                               // gradient
            yVals[j] = yL + dydx * (x - xL);
        }
        return yVals;                                         // linear interpolation
    }
    }  // namespace math
}  // namespace stan

#endif
