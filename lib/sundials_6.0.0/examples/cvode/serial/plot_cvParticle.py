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
# matplotlib-based plotting script for cvPraticle_dns.c example
# ------------------------------------------------------------------------------

# imports
import argparse
import numpy as np
import matplotlib.pyplot as plt

# command line options
parser = argparse.ArgumentParser(description='Plots cvPraticle_dns output')
parser.add_argument('--sfile', type=str,
                    default='cvParticle_solution.txt',
                    help='solution output file to read')
parser.add_argument('--efile', type=str,
                    default='cvParticle_error.txt',
                    help='error output file to read')
parser.add_argument('--alpha', type=float, nargs=1,
                    default=1.0,
                    help='set a non-default alpha value')
parser.add_argument('--slim', type=float, nargs=2,
                    help='x and y limits for solution plot')
parser.add_argument('--eylim', type=float, nargs=2,
                    help='y limits for error plot')

# parse inputs
args = parser.parse_args()

# read solution output file
data = np.loadtxt(args.sfile, dtype=np.double)

# extract times and positions
t = data[:, 0]
x = data[:, 1]
y = data[:, 2]

# unit circle
tt = np.linspace(0,np.pi*2,10000)
xt = np.cos(tt)
yt = np.sin(tt)

# plot solution
fig, ax = plt.subplots()
plt.plot(xt, yt, color='black', linestyle='--')
plt.scatter(x, y, color='red')

if (args.slim):
    plt.xlim((args.slim[0], args.slim[1]))
    plt.ylim((args.slim[0], args.slim[1]))

plt.xlabel('x')
plt.ylabel('y')
plt.title('Solution')
ax.set_aspect('equal')

# true solution
xt = np.cos(args.alpha * t)
yt = np.sin(args.alpha * t)

# plot solution
fig, ax = plt.subplots()
plt.plot(t, x, linestyle='-', label='x')
plt.plot(t, xt, linestyle='--', label='x true')
plt.plot(t, y, linestyle='-', label='y')
plt.plot(t, yt, linestyle='--', label='y true')

plt.xlabel('t')
plt.ylabel('position')
plt.title('Particle Position Over Time')
plt.legend(loc='lower right')

# read error output file
data = np.loadtxt(args.efile, dtype=np.double)

# extract times, position errors, and constraint error
t = data[:, 0]
xerr = np.absolute(data[:, 1])
yerr = np.absolute(data[:, 2])
cerr = np.absolute(data[:, 3])

# plot solution
fig, ax = plt.subplots()
plt.semilogy(t, xerr, label='x err')
plt.semilogy(t, yerr, label='y err')
plt.semilogy(t, cerr, label='c err')

if (args.eylim):
    plt.ylim((args.eylim[0], args.eylim[1]))

plt.xlabel('time')
plt.ylabel('error')
plt.legend(loc='lower right')
plt.title('Error in position and constraint')
plt.grid()

# display plots
plt.show()

##### end of script #####
