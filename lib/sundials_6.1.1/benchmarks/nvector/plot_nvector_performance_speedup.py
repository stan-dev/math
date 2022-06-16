#!/usr/bin/env python
# -----------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
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
# This script plots the output from test_nector_performance_* and assumes:
#   1. vector lengths are powers of two starting from 0, and
#   2. output files are named: output_nelem_nvec_nsum_ntest_timing.txt
# where nelem is the number of elements in the vector, nvec is the nuber of
# vectors, nsum is the number of sums, ntest is the number of tests, and timing
# indicates if timing was enabled.
# -----------------------------------------------------------------------------

def main():

    import argparse
    import os, sys
    import shlex
    import glob

    import numpy as np
    import scipy.stats as st

    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mtick

    parser = argparse.ArgumentParser(
        description='Plot data from NVector performance tests')

    parser.add_argument('op', type=str,
                        help='Which NVector operation to plot')

    parser.add_argument('datadir', type=str,
                        help='Directory where test output files are located')

    parser.add_argument('--noplots', dest='noplots', action='store_true',
                        help='Turn on plots for time vs number of elements')

    parser.add_argument('--logx', dest='logx', action='store_true',
                        help='Generate plots for speedup with log scale for the x axis (number of elements')

    parser.add_argument('--fused', dest='fused', action='store_true',
                        help='Operation is a fused op')

    parser.add_argument('--show', dest='show', action='store_true',
                        help='Display plots rather than saving to file')

    parser.add_argument('--debug', dest='debug', action='store_true',
                        help='Turn on debugging output')

    # parse command line args
    args = parser.parse_args()

    if (args.debug):
        print(args)

    # check for test data directory
    if (not os.path.isdir(args.datadir)):
        print("ERROR:",args.datadir,"does not exist")
        sys.exit()

    # sort output files
    output_baseline = sorted(glob.glob(args.datadir+'/output*-old.log'))
    output_new      = sorted(glob.glob(args.datadir+'/output*-new.log'))
    output = output_baseline + output_new

    if (args.debug):
        print("output files")
        print(len(output))
        for i in range(len(output)):
            print(output[i])

    # figure out vector sizes, number of vectors, and number of sums
    nelem = []
    nvec  = []
    nsum  = []
    ntest = []

    # parse file names to get input parameters
    for f in output:

        split_fout = f.split("/")[-1]
        split_fout = split_fout.split("_")

        ne = int(split_fout[1])
        nv = int(split_fout[2])
        ns = int(split_fout[3])
        nt = int(split_fout[4])

        if (not ne in nelem):
            nelem.append(ne)

        if (not nv in nvec):
            nvec.append(nv)

        if (not ns in nsum):
            nsum.append(ns)

        if (not nt in ntest):
            ntest.append(nt)

    if (len(ntest) != 1):
        print("Warning: Unequal numbers of tests")

    nelem.sort()

    if (args.debug):
        print("nelem:",nelem, len(nelem))
        print("nvec: ",nvec,  len(nvec))
        print("nsum: ",nsum,  len(nsum))
        print("ntest:",ntest, len(ntest))

    # allocate numpy arrays for timing data
    avg_denom  = np.zeros([len(nvec), len(nelem)])
    sdev_denom = np.zeros([len(nvec), len(nelem)])
    avg_numer  = np.zeros([len(nvec), len(nelem)])
    sdev_numer = np.zeros([len(nvec), len(nelem)])
    avg_ratio = np.zeros([len(nvec), len(nelem)])

    # read 'baseline' files
    for f in output_baseline:
        if (args.debug):
            print("Reading:",f)
        # get test inputs from file name
        split_fout = f.split("/")[-1]
        split_fout = split_fout.split("_")
        ne = int(split_fout[1])
        nv = int(split_fout[2])
        ns = int(split_fout[3])
        with open(f) as fout:
            for line in fout:
                # split line into list
                split_line = shlex.split(line)
                # skip blank lines
                if (not split_line):
                    continue
                # tests finished, stop reading file
                if (split_line[0] == "Finished"):
                    break
                # check if the operation is the one we want and get data
                if (args.op == split_line[0]):
                    i = nvec.index(nv)
                    j = nelem.index(ne)
                    avg_numer[i][j]  = float(split_line[1])
                    sdev_numer[i][j] = float(split_line[2])

    # read output files
    for f in output_new:
        if (args.debug):
            print("Reading:",f)
        # get test inputs from file name
        split_fout = f.split("/")[-1]
        split_fout = split_fout.split("_")
        ne = int(split_fout[1])
        nv = int(split_fout[2])
        ns = int(split_fout[3])
        with open(f) as fout:
            for line in fout:
                # split line into list
                split_line = shlex.split(line)
                # skip blank lines
                if (not split_line):
                    continue
                # tests finished, stop reading file
                if (split_line[0] == "Finished"):
                    break
                # check if the operation is the one we want and get data
                if (args.op == split_line[0]):
                    i = nvec.index(nv)
                    j = nelem.index(ne)
                    avg_denom[i][j]  = float(split_line[1])
                    sdev_denom[i][j] = float(split_line[2])
                    avg_ratio[i][j] = avg_numer[i][j] / avg_denom[i][j]

    # --------------------------------------------------------------------------
    # Confidence Interval
    # --------------------------------------------------------------------------

    # allocate arrays for the upper and lower bounds of the confidence interval
    lower_denom   = np.zeros([len(nvec), len(nelem)])
    upper_denom   = np.zeros([len(nvec), len(nelem)])
    lower_numer = np.zeros([len(nvec), len(nelem)])
    upper_numer = np.zeros([len(nvec), len(nelem)])

    # critical value for 99% confidence interval
    if (ntest[0] < 30):
        # student's t distribution
        cv = st.t.interval(0.99, ntest[0]-1)[1]
    else:
        # normal distribution
        cv = st.norm.ppf(0.995)

    # confidence intervals
    cdev_denom  = cv * sdev_denom / np.sqrt(ntest[0])
    lower_denom = avg_denom - cdev_denom
    upper_denom = avg_denom + cdev_denom

    cdev_numer  = cv * sdev_numer / np.sqrt(ntest[0])
    lower_numer = avg_numer - cdev_numer
    upper_numer = avg_numer + cdev_numer

    # check if the new average times are within the baseline confidence interval
    denom_in = np.where(np.logical_and(avg_denom < upper_numer,
                                       avg_denom > lower_numer))

    # check if the baseline average times are within the new confidence interval
    numer_in = np.where(np.logical_and(avg_numer < upper_denom,
                                         avg_numer > lower_denom))

    # get which numbers of vectors and elements for new tests are in the
    # confidence interval of the baseline times
    df = np.zeros([len(denom_in[0])])
    vf = np.zeros([len(denom_in[0])])
    ef = np.zeros([len(denom_in[0])])

    for i in range(len(denom_in[0])):
        vf[i] = nvec[denom_in[0][i]]
        ef[i] = np.log2(nelem[denom_in[1][i]])
        df[i] = 1

    if (args.debug):
        print('vf:', vf)
        print('ef:', ef)

    # get which numbers of vectors and elements for baseline tests are in the
    # confidence interval of the new times
    du = np.zeros([len(numer_in[0])])
    vu = np.zeros([len(numer_in[0])])
    eu = np.zeros([len(numer_in[0])])

    for i in range(len(numer_in[0])):
        vu[i] = nvec[numer_in[0][i]]
        eu[i] = np.log2(nelem[numer_in[1][i]])
        du[i] = 1

    if (args.debug):
        print('vu:', vu)
        print('eu:', eu)

    # --------------------------------------------------------------------------
    # Output ratios
    # --------------------------------------------------------------------------

    # output the matrix of ratios in reverse row order to match the heatmap
    print(args.op)

    print("avg. new")
    for i in reversed(range(len(nvec))):
        print('%2d' % int(i+1), str(avg_denom[i]).replace('\n', ''))
    print()

    print("avg. baseline")
    for i in reversed(range(len(nvec))):
        print('%2d' % int(i+1), str(avg_numer[i]).replace('\n', ''))
    print()

    print("avg. ratio (speedup)")
    for i in reversed(range(len(nvec))):
        print('%2d' % int(i+1), str(avg_ratio[i]).replace('\n', ''))
    print()

    # --------------------------------------------------------------------------
    # Speedup v. Number of Elements Plots
    # --------------------------------------------------------------------------
    if (not args.noplots):

        colors = ['#000000','#a6cee3','#1f78b4','#b2df8a','#33a02c',
                  '#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6',
                  '#6a3d9a','#ffff99','#b15928']

        hatch = [ '/','\\','-','+','x','o','O','.','*']

        # --------------------------------------------------------------------------
        # Combined Number of Vectors Plots
        # --------------------------------------------------------------------------
        fig = plt.figure()
        ax = fig.add_subplot(111)

        if args.fused:
            indices = range(0,len(nvec))
        else:
            indices = range(len(nvec)-1,len(nvec))

        for i in indices:
            lab = 'num. vecs %d' % nvec[i]
            if (args.logx):
                ax.plot(nelem, avg_ratio[i],
                        color=colors[i], linestyle='-', label=lab)
                ax.set_xscale('log')
            else:
                ax.plot(nelem, avg_ratio[i],
                        color=colors[i], linestyle='-', label=lab)
                # # plot confidence interval
                # ax.fill_between(nelem, lower_denom[i], upper_denom[i],
                #                 color=colors[i], alpha=0.3)
                # ax.fill_between(nelem, lower_numer[i], upper_numer[i],
                #                 color=colors[i], hatch='.', alpha=0.3)

        ax.legend()
        ax.grid()

        plt.title('Average Speedup \n'+args.op)
        plt.xlabel('vector length')
        plt.ylabel('speedup (baseline/new)')

        if (args.show):
            plt.show()
        else:
            if (args.logx):
                fname=args.op+'-nvec-all-logx.pdf'
            else:
                fname=args.op+'-nvec-all.pdf'
                plt.ticklabel_format(axis='both',style='sci')
            plt.savefig(fname)
        plt.close()


# ===============================================================================

if __name__ == "__main__":
    main()

# EOF
