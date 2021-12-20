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
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


# load pickled data
def load_data(file):
    data = np.load(file)
    m = data['mesh']
    t = data['t']
    u = data['u']
    v = data['v']
    w = data['w']

    hx = m[0,1] - m[0,0]
    hy = m[1,1] - m[1,0]
    hz = m[2,1] - m[2,0]

    return { 'm': m, 'h': (hx,hy,hz), 't': t, 'u': u, 'v': v, 'w': w }


# grid function norm
def norm_3Dgrid(h, x, q=1):
    hx,hy,hz = h
    return (hx*hy*hz*np.sum(np.abs(x)**q, axis=(1,2,3)))**(1/q)


# computer order of accuracy p
def calc_order(h1, Eh1, h2, Eh2):
    return np.log( Eh1/Eh2 ) / np.log( np.prod(h1)/np.prod(h2) )


# load data files
h_over_8 = load_data('middle-h/output-with-h-1.04e-02.npz')
h_over_4 = load_data('large-h/output-with-h-2.08e-02.npz')
# h_over_2 = load_data('larger-h/output-with-h-4.16e-02.npz')
h_over_1 = load_data('largest-h/output-with-h-8.33e-02.npz')

for component in ['u', 'v', 'w']:
    # Restrict reference to the coarsest grid
    ref = h_over_8[component][:,::8,::8,::8]

    # Now compute E(h) = ||U(h) - \bar{U}(h)|| using the grid-function norm
    Eh_over_4 = norm_3Dgrid(h_over_4['h'], h_over_4[component][:,::4,::4,::4] - ref)
    Eh_over_1 = norm_3Dgrid(h_over_1['h'], h_over_1[component][:,:,:,:] - ref)

    # Compute order p as in O(h^p)
    p = calc_order(h_over_1['h'], Eh_over_1, h_over_4['h'], Eh_over_4)
    print('min p for %s component: %.4f' % (component, np.min(p)))

    # Plot error across time
    plt.figure()
    plt.plot(h_over_8['t'], Eh_over_4, 'r-')
    plt.plot(h_over_8['t'], Eh_over_1, 'b-')
    plt.ylabel('||E(hx,hy,hz)||')
    plt.xlabel('time')
    plt.savefig('error-in-time-plot-%s.png' % component)

    # Plot error norm with respect to h
    plt.figure()
    x = np.array([np.prod(h_over_4['h']), np.prod(h_over_1['h'])])
    plt.plot(x, x, 'k-')
    plt.plot(x, x**2, 'k-')
    plt.plot(x, [np.linalg.norm(Eh_over_4, np.Inf), np.linalg.norm(Eh_over_1, np.Inf)], 'r-')
    plt.legend(['1st order', '2nd order', 'actual'])
    plt.ylabel('|| ||E(hx,hy,hz)|| ||_inf')
    plt.xlabel('hx * hy * hz')
    plt.yscale('log')
    plt.xscale('log')
    plt.savefig('error-plot-%s.png' % component)
