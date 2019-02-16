#!/usr/bin/env python
# ------------------------------------------------------------
# Programmer(s):  Daniel R. Reynolds @ SMU
# ------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2019, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ------------------------------------------------------------
# matplotlib-based plotting script for brusselator1D.c example

# imports
import sys
import pylab as plt
import numpy as np

# load mesh data file
mesh = np.loadtxt('bruss_FEM_mesh.txt', dtype=np.double)

# load solution data files
udata = np.loadtxt('bruss_FEM_u.txt', dtype=np.double)
vdata = np.loadtxt('bruss_FEM_v.txt', dtype=np.double)
wdata = np.loadtxt('bruss_FEM_w.txt', dtype=np.double)

# determine number of time steps, mesh size
nt,nx = np.shape(udata)

# determine min/max values
umin = 0.9*udata.min()
umax = 1.1*udata.max()
vmin = 0.9*vdata.min()
vmax = 1.1*vdata.max()
wmin = 0.9*wdata.min()
wmax = 1.1*wdata.max()
minval = np.array([umin, vmin, wmin]).min()
maxval = np.array([umax, vmax, wmax]).max()

# plot the mesh
plt.figure(1)
plt.plot(mesh,0.0*mesh,'o')
plt.xlabel('x')
plt.title('FEM mesh')
plt.savefig('brusselator1D_FEM_mesh.png')

# generate plots of results
for tstep in range(nt):

    # set string constants for output plots, current time, mesh size
    pname = 'brusselator1D_FEM.' + repr(tstep).zfill(3) + '.png'
    tstr  = repr(tstep)
    nxstr = repr(nx)

    # plot current solution and save to disk
    plt.figure(1)
    plt.plot(mesh,udata[tstep,:],label="u")
    plt.plot(mesh,vdata[tstep,:],label="v")
    plt.plot(mesh,wdata[tstep,:],label="w")
    plt.xlabel('x')
    plt.ylabel('solution')
    plt.title('Solutions at output ' + tstr + ', mesh = ' + nxstr)
    plt.axis((0.0, 1.0, minval, maxval))
    plt.grid()
    plt.legend(loc='upper right', shadow=True)
    plt.savefig(pname)
    plt.close()


##### end of script #####
