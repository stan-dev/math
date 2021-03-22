#!/usr/bin/env python
# ------------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
# ------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2021, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ------------------------------------------------------------------------------
# matplotlib-based plotting script for cvPendulum_dns.c example
# ------------------------------------------------------------------------------

# imports
import argparse
import numpy as np
import matplotlib.pyplot as plt

# command line options
parser = argparse.ArgumentParser(description='Plots cvPendulum_dns output')
parser.add_argument('sfile', type=str,
                    help='solution output file to read')

# parse inputs
args = parser.parse_args()

# read solution output file
data = np.loadtxt(args.sfile, dtype=np.double)

# extract times, positions, and velocities
t  = data[:, 0]
x  = data[:, 1]
y  = data[:, 2]
vx = data[:, 3]
vy = data[:, 4]

# lower half of unit circle
tt = np.linspace(np.pi, 2*np.pi, 10000)
xt = np.cos(tt)
yt = np.sin(tt)

# plot solution in xy plane
fig, ax = plt.subplots()
ax.axhline(y=0, color='black', linestyle='--')
ax.axvline(x=0, color='black', linestyle='--')
plt.plot(xt, yt, color='black', linestyle='--')
plt.scatter(x, y, color='red')

plt.xlabel('x')
plt.ylabel('y')
plt.title('Pendulum')
ax.set_aspect('equal')

# plot position over time
fig, ax = plt.subplots()
ax.axhline(y=0, color='black', linestyle='--')
plt.plot(t, x, label='x')
plt.plot(t, y, label='y')

plt.xlabel('t')
plt.ylabel('position')
plt.title('Pendulum Position')
plt.legend()

# plot velocity over time
fig, ax = plt.subplots()
ax.axhline(y=0, color='black', linestyle='--')
plt.plot(t, vx, label='$v_x$')
plt.plot(t, vy, label='$v_y$')

plt.xlabel('t')
plt.ylabel('velocity')
plt.title('Pendulum Velocity')
plt.legend()

# display plots
plt.show()

##### end of script #####
