#!/usr/bin/env python
# ----------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ SMU
# ----------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2019, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ----------------------------------------------------------------
# matplotlib-based plotting script for heat1D.c example

# imports
import sys
import pylab as plt
import numpy as np

# load mesh data file as list of NumPy arrays
inp = open('heat_mesh.txt').readlines()
mesh = []
for line in inp:
    mesh.append(np.array(str.split(line), dtype=np.double))


# load solution data file as list of NumPy arrays
inp = open('heat1D.txt').readlines()
data = []
for line in inp:
    data.append(np.array(str.split(line), dtype=np.double))

# determine number of time steps
nt  = len(mesh)
nt2 = len(data)
if (nt != nt2):
    sys.exit('plot_heat1D_adapt.py error: data and mesh files have different numbers of time steps')

# determine minimum/maximum temperature
mintemp = 0.0
maxtemp = 0.0
for tstep in range(nt):
    mx = data[tstep].max()
    if (mx > maxtemp):
        maxtemp = mx
    mn = data[tstep].min()
    if (mn < mintemp):
        mintemp = mn
if (maxtemp > 0.0):
    maxtemp *= 1.1
else:
    maxtemp *= 0.9
if (mintemp > 0.0):
    mintemp *= 0.9
else:
    mintemp *= 1.1


# generate plots of results
for tstep in range(nt):

    # set string constants for output plots, current time, mesh size
    pname = 'heat1d.' + repr(tstep).zfill(3) + '.png'
    tstr  = repr(tstep)
    nxstr = repr(len(data[tstep]))

    # plot current solution and save to disk
    plt.figure(1)
    plt.plot(mesh[tstep],data[tstep],'-o')
    plt.xlabel('x')
    plt.ylabel('solution')
    plt.title('u(x) at output ' + tstr + ', mesh = ' + nxstr)
    plt.axis((0.0, 1.0, mintemp, maxtemp))
    plt.grid()
    plt.savefig(pname)
    plt.close()


##### end of script #####
