#!/usr/bin/env python
# -----------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
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

    parser.add_argument('--timevelem', dest='timevelem', action='store_true',
                        help='Turn on plots for time vs number of elements')

    parser.add_argument('--noheatmap', dest='heatmap', action='store_false',
                        help='Turn off heatmap plots')

    parser.add_argument('--loglog', dest='loglog', action='store_true',
                        help='Generate loglog plots for time vs number of elements')

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
    output = sorted(glob.glob(args.datadir+'/output*.txt'))

    # if (args.debug):
    #     print("output files")
    #     print(len(output))
    #     for i in range(len(output)):
    #         print(output[i])

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

    if (args.debug):
        print("nelem:",nelem, len(nelem))
        print("nvec: ",nvec,  len(nvec))
        print("nsum: ",nsum,  len(nsum))
        print("ntest:",ntest, len(ntest))

    # allocate numpy arrays for timing data
    avg_fused  = np.zeros([len(nvec), len(nelem)])
    sdev_fused = np.zeros([len(nvec), len(nelem)])

    avg_unfused  = np.zeros([len(nvec), len(nelem)])
    sdev_unfused = np.zeros([len(nvec), len(nelem)])

    avg_ratio = np.zeros([len(nvec), len(nelem)])

    # NVEC = np.zeros([len(nvec), len(nelem)])
    # NELM = np.zeros([len(nvec), len(nelem)])

    # read output files
    for f in output:

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

                    # NVEC[i][j] = nv
                    # NELM[i][j] = ne

                    avg_fused[i][j]  = float(split_line[1])
                    sdev_fused[i][j] = float(split_line[2])

                    avg_unfused[i][j]  = float(split_line[5])
                    sdev_unfused[i][j] = float(split_line[6])

                    avg_ratio[i][j] = avg_fused[i][j] / avg_unfused[i][j]

    if (args.debug):
        print(avg_fused)
        print(avg_unfused)
        print(avg_ratio)
        print(sdev_fused)
        print(sdev_unfused)

    # --------------------------------------------------------------------------
    # Confidence Interval
    # --------------------------------------------------------------------------

    # allocate arrays for the upper and lower bounds of the confidence interval
    lower_fused   = np.zeros([len(nvec), len(nelem)])
    upper_fused   = np.zeros([len(nvec), len(nelem)])
    lower_unfused = np.zeros([len(nvec), len(nelem)])
    upper_unfused = np.zeros([len(nvec), len(nelem)])

    # critical value for 99% confidence interval
    if (ntest[0] < 30):
        # student's t distribution
        cv = st.t.interval(0.99, ntest[0]-1)[1]
    else:
        # normal distribution
        cv = st.norm.ppf(0.995)

    # confidence intervals
    cdev_fused  = cv * sdev_fused / np.sqrt(ntest[0])
    lower_fused = avg_fused - cdev_fused
    upper_fused = avg_fused + cdev_fused

    cdev_unfused  = cv * sdev_unfused / np.sqrt(ntest[0])
    lower_unfused = avg_unfused - cdev_unfused
    upper_unfused = avg_unfused + cdev_unfused

    # check if the fused average times are within the unfused confidence interval
    fused_in = np.where(np.logical_and(avg_fused < upper_unfused,
                                       avg_fused > lower_unfused))

    # check if the unfused average times are within the fused confidence interval
    unfused_in = np.where(np.logical_and(avg_unfused < upper_fused,
                                         avg_unfused > lower_fused))

    # get which numbers of vectors and elements for fused tests are in the
    # confidence interval of the unfused times
    df = np.zeros([len(fused_in[0])])
    vf = np.zeros([len(fused_in[0])])
    ef = np.zeros([len(fused_in[0])])

    for i in range(len(fused_in[0])):
        vf[i] = nvec[fused_in[0][i]]
        ef[i] = np.log2(nelem[fused_in[1][i]])
        df[i] = 1

    if (args.debug):
        print(vf)
        print(ef)

    # get which numbers of vectors and elements for unfused tests are in the
    # confidence interval of the fused times
    du = np.zeros([len(unfused_in[0])])
    vu = np.zeros([len(unfused_in[0])])
    eu = np.zeros([len(unfused_in[0])])

    for i in range(len(unfused_in[0])):
        vu[i] = nvec[unfused_in[0][i]]
        eu[i] = np.log2(nelem[unfused_in[1][i]])
        du[i] = 1

    if (args.debug):
        print(vu)
        print(eu)

    # --------------------------------------------------------------------------
    # Output ratios
    # --------------------------------------------------------------------------

    # output the matrix of ratios in reverse row order to match the heatmap
    print(args.op)

    # print(avg_fused)
    # for i in reversed(range(len(nvec))):
    #     print('%2d' % int(i+1), str(avg_fused[i]).replace('\n', ''))
    # print

    # print(avg_unfused)
    # for i in reversed(range(len(nvec))):
    #     print('%2d' % int(i+1), str(avg_unfused[i]).replace('\n', ''))
    # print

    # print(NVEC)
    # print(NELM)
    # print(avg_ratio)
    for i in reversed(range(len(nvec))):
        print('%2d' % int(i+1), str(avg_ratio[i]).replace('\n', ''))
    print

    # --------------------------------------------------------------------------
    # Heat Map
    # --------------------------------------------------------------------------
    if (args.heatmap):

        x = np.arange(len(nelem)+1)-0.5 # x = log2(number of elements) = 0,1,2,...
        y = np.arange(len(nvec)+1)+1.5  # y = number of vectors = 2,3,4,...
        # y = np.arange(len(nvec)+1)+0.5  # y = number of vectors = 1,2,3,...
        X, Y = np.meshgrid(x, y)

        if (args.debug):
            print(x)
            print(y)

        # center the color bar around 1 (if possible)
        rmax = np.amax(avg_ratio)
        rmin = np.amin(avg_ratio)

        ext = 'neither'
        if (rmin > 1):
            cmap='Reds'
            norm = mpl.colors.Normalize(vmin=rmin, vmax=min(rmax,2))
            v = np.linspace(rmin, min(rmax,2), 10, endpoint=True)
            if (rmax > 2):
                ext = 'max'
        else:
            cmap='seismic'
            if (rmax-1 > 1):
                rrange = 1
                ext = 'max'
            else:
                rrange = max(abs(rmax-1),abs(rmin-1))

            v1 = np.linspace(1-rrange, 1, 5, endpoint=True)
            v2 = np.linspace(1, 1+rrange, 5, endpoint=True)
            v = np.append(v1,v2[1:])
            norm = mpl.colors.Normalize(vmin=1-rrange, vmax=1+rrange)

        # plot heatmap
        plt.pcolormesh(X, Y, avg_ratio, cmap=cmap, norm=norm)
        clb = plt.colorbar(ticks=v, extend=ext)
        clb.ax.set_title('Max = {0:.2f}\nMin = {1:.2f}'.format(rmax,rmin))

        # aff markers to indicate if the average time falls in a confidence interval
        plt.scatter(ef,vf,s=40,marker='^',c=df,label='fused')
        plt.scatter(eu,vu,s=40,marker='v',c=du,label='unfused')
        plt.legend(loc=9, bbox_to_anchor=(0.5, -0.1), ncol=2)

        # add legend for scatter plot
        art = []
        lgd = plt.legend(loc='lower right', bbox_to_anchor=(1.34, -0.17))
        art.append(lgd)

        # add labels and title
        plt.xticks(np.log2(nelem))
        plt.yticks(nvec)
        plt.xlabel('log2(num elements)')
        plt.ylabel('num vectors')
        plt.title('avg fused time / avg unfused time \n'+args.op)

        # display or save figure
        if (args.show):
            plt.show()
        else:
            plt.savefig(args.op+'-heatmap.pdf',
                        additional_artists=art,
                        bbox_inches="tight")
        plt.close()

    # --------------------------------------------------------------------------
    # Time vs Number of Elements Plots
    # --------------------------------------------------------------------------
    if (args.timevelem):

        colors = ['#000000','#a6cee3','#1f78b4','#b2df8a','#33a02c',
                  '#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6',
                  '#6a3d9a','#ffff99','#b15928']

        hatch = [ '/','\\','-','+','x','o','O','.','*']

        # --------------------------------------------------------------------------
        # Combined Number of Vectors Plots
        # --------------------------------------------------------------------------
        fig = plt.figure()
        ax = fig.add_subplot(111)

        for nv in nvec:

            i = nvec.index(nv)

            if (args.loglog):
                ax.loglog(nelem, avg_fused[i],
                          color=colors[i], linestyle='-', label=nv)
                ax.loglog(nelem, avg_unfused[i],
                          color=colors[i], linestyle='--', label=None)
            else:
                ax.plot(nelem, avg_fused[i],
                        color=colors[i], linestyle='-', label=nv)
                ax.plot(nelem, avg_unfused[i],
                        color=colors[i], linestyle='--', label=None)

            # plot confidence interval
            ax.fill_between(nelem, lower_fused[i], upper_fused[i],
                            color=colors[i], alpha=0.3)
            ax.fill_between(nelem, lower_unfused[i], upper_unfused[i],
                            color=colors[i], hatch='.', alpha=0.3)

        ax.legend()
        ax.grid()

        plt.title('Average Time Fused vs Unfused \n'+args.op)
        plt.xlabel('vector length')
        plt.ylabel('time (s)')

        if (args.show):
            plt.show()
        else:
            if (args.loglog):
                fname=args.op+'-nvec-all-loglog.pdf'
            else:
                fname=args.op+'-nvec-all.pdf'
                plt.ticklabel_format(axis='both',style='sci')
            plt.savefig(fname)
        plt.close()

        # --------------------------------------------------------------------------
        # Individual Number of Vectors Plots
        # --------------------------------------------------------------------------

        for nv in nvec:
            fig = plt.figure()
            ax  = fig.add_subplot(111)
            idx = nvec.index(nv)

            # plot run times
            if (args.loglog):
                ax.loglog(nelem, avg_fused[idx],
                          color='red', linestyle='-', label='Fused')
                ax.loglog(nelem, avg_unfused[idx],
                          color='blue', linestyle='--', label='Unfused')
            else:
                ax.plot(nelem, avg_fused[idx],
                        color='red', linestyle='-', label='Fused')
                ax.plot(nelem, avg_unfused[idx],
                        color='blue', linestyle='--', label='Unfused')

            # plot confidence intervals
            ax.fill_between(nelem, lower_fused[idx], upper_fused[idx],
                            color='red', alpha=0.2)
            ax.fill_between(nelem, lower_unfused[idx], upper_unfused[idx],
                            color='blue', hatch='.', alpha=0.2)

            ax.legend()
            ax.grid()

            plt.title('Average Time Fused vs Unfused with '+str(nv)+' vectors\n'+args.op)
            plt.xlabel('vector length')
            ax.set_ylabel('time (s)')

            if (args.show):
                plt.show()
            else:
                if (args.loglog):
                    fname=args.op+'-nvec-'+str(nv)+'-loglog.pdf'
                else:
                    fname=args.op+'-nvec-'+str(nv)+'.pdf'
                    plt.ticklabel_format(axis='both',style='sci')
                plt.savefig(fname)
            plt.close()

# ===============================================================================

if __name__ == "__main__":
    main()

# EOF

