#!/usr/bin/env python
# -----------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2021, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -----------------------------------------------------------------------------

# imports
import glob
import sys
import pylab as plt
import pandas as pd
import numpy as np

# load mesh data file
mesh = np.loadtxt('mesh.txt', dtype=np.double)
# X,Y,Z = np.meshgrid(mesh[0,:], mesh[1,:], mesh[2,:])

# calculate h
hx = mesh[0,1] - mesh[0,0]
hy = mesh[1,1] - mesh[1,0]
hz = mesh[2,1] - mesh[2,0]
nx = len(mesh[0,:])
ny = len(mesh[1,:])
nz = len(mesh[2,:])

print("nx, ny, nz = %d, %d, %d" % (nx, ny, nz))
print("hx, hy, hz = %g, %g, %g" % (hx, hy, hz))

# load output time file
times = np.loadtxt('t.000000.txt', dtype=np.double)

# load solution data files
ufiles = glob.glob('u.' + ('[0-9]'*6) + '.txt'); ufiles.sort()
vfiles = glob.glob('v.' + ('[0-9]'*6) + '.txt'); vfiles.sort()
wfiles = glob.glob('w.' + ('[0-9]'*6) + '.txt'); wfiles.sort()
udata = []
vdata = []
wdata = []

sys.stdout.write("reading 1/%d...\r" % len(ufiles))
sys.stdout.flush()
for idx in range(0,len(ufiles)):
    sys.stdout.write("reading %d/%d...\r" % (idx+1,len(ufiles)))
    sys.stdout.flush()
    udata.append(pd.read_csv(ufiles[idx], header=None, delimiter=' ', skipinitialspace=True, dtype=np.double))
    vdata.append(pd.read_csv(vfiles[idx], header=None, delimiter=' ', skipinitialspace=True, dtype=np.double))
    wdata.append(pd.read_csv(wfiles[idx], header=None, delimiter=' ', skipinitialspace=True, dtype=np.double))
sys.stdout.write("\n")
sys.stdout.flush()

print("stacking...")
udata = pd.concat(udata, axis=1).to_numpy()
vdata = pd.concat(vdata, axis=1).to_numpy()
wdata = pd.concat(wdata, axis=1).to_numpy()

# reshape data into time,x,y,z arrays
print("reshaping...")
nt = len(times)
udata = np.reshape(udata, (nt, nx, ny, nz))
vdata = np.reshape(vdata, (nt, nx, ny, nz))
wdata = np.reshape(wdata, (nt, nx, ny, nz))

# save data to pickle
print("saving...")
np.savez_compressed('output-with-h-%.2e.npz' % hx, t=times, u=udata, v=vdata, w=wdata, mesh=mesh)

