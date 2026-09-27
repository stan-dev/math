import numpy as np
from scipy.interpolate import CubicHermiteSpline

############################################
# Ground truth for interp1CubicHermiteSetupCoefficients test
############################################
print("\n=== Ground truth for interp1CubicHermiteSetupCoefficients ===")
xk_setup = np.array([0.0, 1.0, 2.5, 4.0])
yk_setup = np.array([10.0, 12.0, 9.0, 11.0])
dydxk_setup = np.array([1.0, 0.5, -0.5, 1.5])

spline_setup = CubicHermiteSpline(xk_setup, yk_setup, dydxk_setup)

def print_setup_coefficients(i):
    coeffs = spline_setup.c[:, i]
    print(f"Interval [{xk_setup[i]}, {xk_setup[i+1]}] coefficients:")
    print(f"a0 = {coeffs[3]}")
    print(f"a1 = {coeffs[2]}")
    print(f"a2 = {coeffs[1]}")
    print(f"a3 = {coeffs[0]}")
    print()

print_setup_coefficients(0)
print_setup_coefficients(1)
print_setup_coefficients(2)

############################################
# Ground truth for interp1CubicHermiteEvaluation test 
############################################
print("\n=== Ground truth for interp1CubicHermiteEvaluation ===")
xk_eval = np.array([0.0, 1.0, 3.0])
yk_eval = np.array([2.0, 4.0, 8.0])
dydxk_eval = np.array([1.0, 2.0, 1.5])

spline_eval = CubicHermiteSpline(xk_eval, yk_eval, dydxk_eval)

print("Evaluation results:")
points = [0.5, 1.2, 2.0]
for point in points:
    print(f"x = {point}: {spline_eval(point):.8f}")

############################################
# Ground truth for interp1CubicHermiteEdgeCases test
############################################
print("\n=== Ground truth for interp1CubicHermiteEdgeCases ===")
# Tightly spaced knots
xk_tight = np.array([0.0, 1e-10])
yk_tight = np.array([1.0, 1.0])
dydxk_tight = np.array([0.0, 0.0])

spline_tight = CubicHermiteSpline(xk_tight, yk_tight, dydxk_tight)
print("Tightly spaced knots:")
print(f"x = 0.5e-10: {spline_tight(0.5e-10)}")
print(f"x = 0.9e-10: {spline_tight(0.9e-10)}")

# Large values
xk_large = np.array([1e10, 2e10, 3e10])
yk_large = np.array([1e20, 2e20, 3e20])
dydxk_large = np.array([1e15, 2e15, 3e15])

spline_large = CubicHermiteSpline(xk_large, yk_large, dydxk_large)
print("\nLarge values:")
print(f"x = 1.5e10: {spline_large(1.5e10)}")
print(f"x = 2.9e10: {spline_large(2.9e10)}")