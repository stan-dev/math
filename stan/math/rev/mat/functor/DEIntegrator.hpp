/*
 * Copyright (c) 2015, John D. Cook
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 * Neither the name of Columbia University nor the names of its
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 * ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
*/
#ifndef DEINTEGRATOR_H
#define DEINTEGRATOR_H

#include "DEIntegrationConstants.h"
#include <cmath>
#include <cfloat>
namespace stan {

  namespace math {

    template<class TFunctionObject>
    class DEIntegrator {
    public:
        /*! Numerical integration in one dimension using the double
         * expontial method of M. Mori.
         *
         * @tparam TFunctionObject Type of f.
         * @param f a functor
         * @param a lower limit of integration, must be double type.
         * @param b upper limit of integration, must be double type.
         * @param targetRelativeError target relative error.
         * @param targetAbsoluteError target absolute error.
         * @param numFunctionEvaluations number of function evaluations used
         * @param errorEstimate estimated error in integration.
         * @return numeric integral of function f.
         */
        static double Integrate
        (const TFunctionObject& f,
         double a,
         double b,
         double targetRelativeError,
         double targetAbsoluteError,
         int& numFunctionEvaluations,
         double& errorEstimate) {
            // Apply the linear change of variables x = ct + d
            // $$\int_a^b f(x) dx = c \int_{-1}^1 f( ct + d ) dt$$
            // c = (b-a)/2, d = (a+b)/2

            double c = 0.5*(b - a);
            double d = 0.5*(a + b);

            return IntegrateCore
            (f,
             c,
             d,
             targetRelativeError,
             targetAbsoluteError,
             numFunctionEvaluations,
             errorEstimate,
             doubleExponentialAbcissas,
             doubleExponentialWeights);
        }

        /*! Numerical integration in one dimension using the double
         * expontial method of M. Mori.
         *
         * @tparam TFunctionObject Type of f.
         * @param f a functor
         * @param a lower limit of integration, must be double type.
         * @param b upper limit of integration, must be double type.
         * @param targetRelativeError target relative error.
         * @param targetAbsoluteError target absolute error.
         * @return numeric integral of function f.
         */
        static double Integrate
        (const TFunctionObject& f,
         double a,
         double b,
         double targetRelativeError,
         double targetAbsoluteError) {
            int numFunctionEvaluations;
            double errorEstimate;
            return Integrate
            (f,
             a,
             b,
             targetRelativeError,
             targetAbsoluteError,
             numFunctionEvaluations,
             errorEstimate);
        }



    private:
        // Integrate f(cx + d) with the given integration constants
        static double IntegrateCore
        (const TFunctionObject& f,
         double c,   // slope of change of variables
         double d,   // intercept of change of variables
         double targetRelativeError,
         double targetAbsoluteError,
         int& numFunctionEvaluations,
         double& errorEstimate,
         const double* abcissas,
         const double* weights) {
            using std::log;
            using std::fabs;

            targetAbsoluteError /= c;

            // Offsets to where each level's integration constants
            // start. The last element is not a beginning but an end.
            int offsets[] = {1, 4, 7, 13, 25, 49, 97, 193};
            int numLevels = sizeof(offsets) / sizeof(int) - 1;

            double newContribution = 0.0;
            double integral = 0.0;
            double previousIntegral = 0.0;
            errorEstimate = DBL_MAX;
            double h = 1.0;
            double previousDelta, currentDelta = DBL_MAX;

            integral = f(c*abcissas[0] + d) * weights[0];
            int i;
            for (i = offsets[0]; i != offsets[1]; ++i)
                integral += weights[i]*(f(c*abcissas[i] + d)
                                        + f(-c*abcissas[i] + d));

            for (int level = 1; level != numLevels; ++level) {
                h *= 0.5;
                newContribution = 0.0;
                for (i = offsets[level]; i != offsets[level+1]; ++i)
                    newContribution += weights[i]*(f(c*abcissas[i] + d)
                                       + f(-c*abcissas[i] + d));
                newContribution *= h;

                // difference in consecutive integral estimates
                previousDelta = currentDelta;
                currentDelta = fabs(0.5*integral - newContribution);
                previousIntegral = integral;
                integral = 0.5*integral + newContribution;

                // Once convergence kicks in, error is approximately
                // squared at each step.
                // Determine whether we're in the convergent region by
                // looking at the trend in the error.
                if (level == 1)
                    // previousDelta meaningless, so cannot check
                    // convergence.
                    continue;

                // Exact comparison with zero is harmless here.
                // Could possibly be replaced with
                // a small positive upper limit on the size of
                // currentDelta, but determining
                // that upper limit would be difficult. At worse, the
                // loop is executed more times than necessary.
                // But no infinite loop can result since there is
                // an upper bound on the loop variable.
                if (currentDelta == 0.0)
                    break;
                // previousDelta != 0 or would have been kicked out
                // previously
                double r = log(currentDelta) / log(previousDelta);

                if (r > 1.9 && r < 2.1) {
                    // If convergence theory applied perfectly,
                    // r would be 2 in the convergence region.
                    // r close to 2 is good enough.
                    // We expect the difference between this integral
                    // estimate
                    // and the next one to be roughly delta^2.
                    errorEstimate = currentDelta*currentDelta;
                } else {
                    // Not in the convergence region.
                    // Assume only that error is decreasing.
                    errorEstimate = currentDelta;
                }

                if (errorEstimate < 0.1*targetAbsoluteError ||
                    fabs(1 - integral / previousIntegral) <
                    targetRelativeError)
                    break;
            }

            numFunctionEvaluations = 2*i - 1;
            errorEstimate *= c;
            return c*integral;
      }
    };
  }
}
#endif
