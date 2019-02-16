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
# matplotlib-based plotting script for ODE examples

# imports
import sys
import pylab as plt
import numpy as np

# load solution data file
data = np.loadtxt('solution.txt', dtype=np.double)

# determine number of time steps, number of fields
nt,nv = np.shape(data)

# extract time array
times = data[:,0]

# parse comment line to determine solution names
f = open('solution.txt', 'r')
commentline = f.readline()
commentsplit = commentline.split()
names = commentsplit[2:]

# create plot
plt.figure()

# add curves to figure
for i in range(nv-1):
    plt.plot(times,data[:,i+1],label=names[i])
plt.xlabel('t')
if (nv > 2):
    plt.ylabel('solutions')
else:
    plt.ylabel('solution')
plt.legend(loc='upper right', shadow=True)
plt.grid()
plt.savefig('solution.png')




##### end of script #####
