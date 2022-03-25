#!/usr/bin/env python
# -----------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2022, Lawrence Livermore National Security
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
from mpl_toolkits.mplot3d import Axes3D
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
    s = np.shape(x)
    return (hx*hy*hz*np.sum(np.abs(x)**q, axis=(1,2,3)))**(1./q)


# load data files
np111 = load_data('np-111/output-with-h-8.33e-02.npz')
np211 = load_data('np-211/output-with-h-8.33e-02.npz')
np311 = load_data('np-311/output-with-h-8.33e-02.npz')
np131 = load_data('np-131/output-with-h-8.33e-02.npz')
np113 = load_data('np-113/output-with-h-8.33e-02.npz')
np911 = load_data('np-911/output-with-h-8.33e-02.npz')
# np133 = load_data('np-133/output-with-h-8.33e-02.npz')
np313 = load_data('np-313/output-with-h-8.33e-02.npz')
np331 = load_data('np-331/output-with-h-8.33e-02.npz')
np333 = load_data('np-333/output-with-h-8.33e-02.npz')
# np666 = load_data('np-666/output-with-h-8.33e-02.npz')

for component in ['u', 'v', 'w']:
    # Reference solution
    ref = np111[component]

    # Now compute E(h) = ||U(h) - \bar{U}(h)|| using the grid-function norm
    E_np211 = norm_3Dgrid(np211['h'], np211[component] - ref)
    E_np311 = norm_3Dgrid(np311['h'], np311[component] - ref)
    E_np131 = norm_3Dgrid(np131['h'], np131[component] - ref)
    E_np113 = norm_3Dgrid(np113['h'], np113[component] - ref)
    E_np911 = norm_3Dgrid(np911['h'], np911[component] - ref)
    # E_np133 = norm_3Dgrid(np133['h'], np133[component] - ref)
    E_np313 = norm_3Dgrid(np313['h'], np313[component] - ref)
    E_np331 = norm_3Dgrid(np331['h'], np331[component] - ref)
    E_np333 = norm_3Dgrid(np333['h'], np333[component] - ref)
    # E_np666 = norm_3Dgrid(np666['h'], np666[component] - ref)

    # Plot error across time
    X, Y = np.meshgrid(np111['m'][0,:], np111['t'])
    # fig = plt.figure()
    # ax = plt.subplot(311, projection='3d')
    # ax.plot_surface(X, Y, np.abs(np911[component][:,:,0,0] - ref[:,:,0,0]))
    # ax = plt.subplot(312, projection='3d')
    # ax.plot_surface(X, Y, np.abs(np911[component][:,0,:,0] - ref[:,0,:,0]))
    # ax = plt.subplot(313, projection='3d')
    # ax.plot_surface(X, Y, np.abs(np911[component][:,0,0,:] - ref[:,0,0,:]))
    plt.plot(np111['t'], E_np211)
    plt.plot(np111['t'], E_np131)
    plt.plot(np111['t'], E_np113)
    plt.plot(np111['t'], E_np911)
    # plt.plot(np111['t'], E_np133)
    plt.plot(np111['t'], E_np313)
    plt.plot(np111['t'], E_np331)
    plt.plot(np111['t'], E_np333)
    # plt.plot(np111['t'], E_np666)
    # plt.legend(['2 1 1', '3 1 1', '1 3 3', '3 1 3', '3 3 1', '3 3 3', '6 6 6'])
    # plt.legend(['3 1 1', '1 3 1', '1 1 3', '9 1 1', '1 3 3', '3 1 3', '3 3 1'])
    plt.ylabel('||E(hx,hy,hz)||')
    plt.xlabel('time')
    plt.savefig('compare-error-plot-%s.png' % component)
