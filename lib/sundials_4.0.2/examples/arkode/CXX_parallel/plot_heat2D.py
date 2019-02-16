#!/usr/bin/env python
# ----------------------------------------------------------------
# Programmer(s):  Daniel R. Reynolds @ SMU
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
# matplotlib-based plotting script for heat2D.cpp example
# ----------------------------------------------------------------

# imports
import sys
import numpy as np
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt

# determine the number of MPI processes used
nprocs=1
for i in range(1000):
    sname = 'heat2d_subdomain.' + repr(i).zfill(3) + '.txt'
    try:
        f = open(sname,'r')
        f.close()
    except IOError:
        nprocs = i
        break

# load subdomain information, store in table
subdomains = np.zeros((nprocs,4), dtype=np.int)
for i in range(nprocs):
    sname = 'heat2d_subdomain.' + repr(i).zfill(3) + '.txt'
    subd = np.loadtxt(sname, dtype=np.int)
    if (i == 0):
        nx = subd[0]
        ny = subd[1]
    else:
        if ((subd[0] != nx) or (subd[1] != ny)):
            sys.exit("error: subdomain files incompatible (clean up and re-run test)")
    subdomains[i,:] = subd[2:6]
    
# load first processor's data, and determine total number of time steps
data = np.loadtxt('heat2d.000.txt', dtype=np.double)
nt = np.shape(data)[0]

# create empty array for all solution data
results = np.zeros((nt,ny,nx))

# insert first processor's data into results array
istart = subdomains[0,0]
iend = subdomains[0,1]
jstart = subdomains[0,2]
jend = subdomains[0,3]
nxl = iend-istart+1
nyl = jend-jstart+1
for i in range(nt):
    results[i,jstart:jend+1,istart:iend+1] = np.reshape(data[i,:], (nyl,nxl))
    
# iterate over remaining data files, inserting into output
if (nprocs > 1):
    for isub in range(1,nprocs):
        data = np.loadtxt('heat2d.' + repr(isub).zfill(3) + '.txt', dtype=np.double)
        # check that subdomain has correct number of time steps
        if (np.shape(data)[0] != nt):
            sys.exit('error: subdomain ' + isub + ' has an incorrect number of time steps')
        istart = subdomains[isub,0]
        iend = subdomains[isub,1]
        jstart = subdomains[isub,2]
        jend = subdomains[isub,3]
        nxl = iend-istart+1
        nyl = jend-jstart+1
        for i in range(nt):
            results[i,jstart:jend+1,istart:iend+1] = np.reshape(data[i,:], (nyl,nxl))

# determine extents of plots
maxtemp = 1.1*results.max()
mintemp = 0.9*results.min()

# generate plots of results
kx = 0.5
ky = 0.75
kprod = (kx+4.0*ky)*np.pi**2
dt = 0.015
for tstep in range(nt):

    # set string constants for output plots, current time, mesh size
    pname = 'heat2d_surf.' + repr(tstep).zfill(3) + '.png'
    cname = 'heat2d_err.' + repr(tstep).zfill(3) + '.png'
    tstr  = repr(tstep)
    nxstr = repr(nx)
    nystr = repr(ny)

    # set x and y meshgrid objects
    xspan = np.linspace(0.0, 1.0, nx)
    yspan = np.linspace(0.0, 1.0, ny)
    X,Y = np.meshgrid(xspan,yspan)

    # plot current solution as a surface, and save to disk
    fig = plt.figure(1)
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, results[tstep,:,:], rstride=1, cstride=1, 
                    cmap=cm.jet, linewidth=0, antialiased=True, shade=True)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlim((mintemp, maxtemp))
    ax.view_init(20,45)
    title('u(x,y) at output ' + tstr + ', mesh = ' + nxstr + 'x' + nystr)
    savefig(pname)
    plt.close()

    # plot error in current solution (as a contour, and save to disk)
    t = tstep*dt;
    at = (1.0 - np.exp(-t*kprod))/kprod
    utrue = at*np.sin(np.pi*X)*np.sin(2.0*np.pi*Y);
    uerr = np.abs(utrue - results[tstep,:,:])
    plt.contourf(xspan,yspan,uerr,15, cmap=plt.cm.jet)
    plt.colorbar()
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Error at output ' + tstr + ', mesh = ' + nxstr + 'x' + nystr)
    plt.savefig(cname)
    plt.close()




##### end of script #####
