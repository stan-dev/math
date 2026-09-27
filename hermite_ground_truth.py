import numpy as np
from scipy.interpolate import CubicHermiteSpline
# generate ground truth for interp1_cubic_hermite_test 
# Can be omitted once tests have been established.

# Shared test inputs
xk = np.array([0.0, 1.0, 2.5, 4.0])
yk = np.array([10.0, 12.0, 9.0, 11.0])
dydxk = np.array([1.0, 0.5, -0.5, 1.5])

# Create spline
spline = CubicHermiteSpline(xk, yk, dydxk)

def print_coefficients(i):
    # Get coefficients for the i-th interval
    # SciPy stores coefficients in reverse order (a3, a2, a1, a0)
    coeffs = spline.c[:, i]
    print(f"Interval [{xk[i]}, {xk[i+1]}] coefficients:")
    print(f"a0 = {coeffs[3]}")  # Constant term
    print(f"a1 = {coeffs[2]}")  # Linear term
    print(f"a2 = {coeffs[1]}")  # Quadratic term
    print(f"a3 = {coeffs[0]}")  # Cubic term
    print()

# Print coefficients for each interval
print_coefficients(0)  # [0.0, 1.0]
print_coefficients(1)  # [1.0, 2.5]
print_coefficients(2)  # [2.5, 4.0]

# Print evaluation points
print("Evaluation results:")
print(f"x = 1.5: {spline(1.5):.8f}")
print(f"x = 3.0: {spline(3.0):.8f}")
# print(f"x = -1.0: {spline(-1.0):.8f} (extrapolated)")
# print(f"x = 5.0: {spline(5.0):.8f} (extrapolated)")