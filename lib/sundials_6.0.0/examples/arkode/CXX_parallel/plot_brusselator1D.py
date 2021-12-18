#!/usr/bin/env python
# ----------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ SMU
# ----------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2021, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ----------------------------------------------------------------
# matplotlib-based plotting script
# ----------------------------------------------------------------

# imports
import glob
import sys
import pylab as plt
import numpy as np

# load mesh data file
mesh = np.loadtxt('mesh.txt', dtype=np.double)

# load output time file
times = np.loadtxt('t.000000.txt', dtype=np.double)

# load solution data files
ufiles = glob.glob('u.' + ('[0-9]'*6) + '.txt'); ufiles.sort()
vfiles = glob.glob('v.' + ('[0-9]'*6) + '.txt'); vfiles.sort()
wfiles = glob.glob('w.' + ('[0-9]'*6) + '.txt'); wfiles.sort()
udata = np.loadtxt(ufiles[0], dtype=np.double)
vdata = np.loadtxt(vfiles[0], dtype=np.double)
wdata = np.loadtxt(wfiles[0], dtype=np.double)
for idx in range(1,len(ufiles)):
    udata = np.hstack((udata, np.loadtxt(ufiles[idx], dtype=np.double)))
    vdata = np.hstack((vdata, np.loadtxt(vfiles[idx], dtype=np.double)))
    wdata = np.hstack((wdata, np.loadtxt(wfiles[idx], dtype=np.double)))

# determine number of time steps, mesh size
nt,nx = np.shape(udata)

# determine min/max values
umin = 0.9*udata.min()
umax = 1.1*udata.max()
vmin = 0.9*vdata.min()
vmax = 1.1*vdata.max()
wmin = 0.9*wdata.min()
wmax = 1.1*wdata.max()
xmax = mesh.max()
minval = np.array([umin, vmin, wmin]).min()
maxval = np.array([umax, vmax, wmax]).max()

# generate plots of results
for tstep in range(nt):

    # set string constants for output plots, current time, mesh size
    pname = 'solution.' + repr(tstep).zfill(3) + '.png'
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
    plt.axis((0.0, xmax, minval, maxval))
    plt.grid()
    plt.legend(loc='upper right', shadow=True)
    plt.savefig(pname)
    plt.close()

# set string constants for output plots, current time, mesh size
pname = 'solution_at_x0.png'
xstr = repr(mesh[0])

# plot current solution and save to disk
plt.figure(1)
plt.plot(times,udata[:,0],label="u")
plt.plot(times,vdata[:,0],label="v")
plt.plot(times,wdata[:,0],label="w")
plt.xlabel('t')
plt.ylabel('solution')
plt.title('Solutions at output at x = '+xstr)
plt.axis((times[0], times[-1], minval, maxval))
plt.grid()
plt.legend(loc='upper right', shadow=True)
plt.savefig(pname)
plt.close()


##### end of script #####
