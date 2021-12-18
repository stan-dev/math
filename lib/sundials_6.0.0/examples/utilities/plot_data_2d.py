#!/usr/bin/env python3
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
# Script to plot data from a 2-dimensional problem on a rectangular domain with
# uniform grid spacing.
#
# The input data file(s) must contain an M x N matrix of values where M is the
# number of output times and the N = 1 + V * X * Y columns contain the output
# data where V is the number of solution variables, X is the number of nodes in
# the x-direction, and Y is the number of nodes in the y-direction.
#
# The first column of each row is the output time for the solution data in the
# remaining columns. The output data must be ordered so the solution variables
# at a given location grouped together and the spatial node ordering varies
# first over x and then y. For example if V = 2, X = 2, and Y = 2 then the
# output file would contain the data
#
# t_0      u_0,0   v_0,0   u_1,0   v_1,0   u_0,1   v_0,1   u_1,1   v_1,1
# t_1      u_0,0   v_0,0   u_1,0   v_1,0   u_0,1   v_0,1   u_1,1   v_1,1
#  .         .       .       .       .       .       .       .       .
#  .         .       .       .       .       .       .       .       .
#  .         .       .       .       .       .       .       .       .
# t_{M-1}  u_0,0   v_0,0   u_1,0   v_1,0   u_0,1   v_0,1   u_1,1   v_1,1
#
# where u_i,j and v_i,j are the solution componetns at node (i,j). In general,
# the n-th solution component at node (i,j) is located at column index
# 1 + V * (X * j + i) + n with n = 0,...,V-1, i = 0,...,X-1, and j = 0,...,Y-1.
#
# Additionally data file(s) should contain the following header comment block
# describing the output data (note comment lines begin with #):
#
#   # nprocs  <number of mpi processes>
#   # nvar    <number of solution variables>
#   # nt      <number of output times>
#   # nx      <number of nodes in the x-direction>
#   # xl      <x-direction lower bound>
#   # xu      <x-direction upper bound>
#   # is      <subdomain global starting x-index>
#   # ie      <subdomain global ending x-index>
#   # ny      <number of nodes in the y-direction>
#   # yl      <y-direction lower bound>
#   # yu      <y-direction upper bound>
#   # js      <subdomain global starting y-index>
#   # je      <subdomain global ending y-index>
#
# With the exception of the subdomain indices, if any of the above values are
# not provided in the header the script will attempt to deduce the necessary
# values and will print a warning or halt with an error.
#
# The comment block may optionally contain additional entries used in creating
# plots:
#   # title    <plot title>
#   # varnames <variables names>
#
# Example usage:
#   plot_data_2d.py data_file.out
#   plot_data_2d.py data_file_proc_0.out data_file_proc_1.out
#   plot_data_2d.py data_file.out --plottype surface-ani
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# main routine
# -----------------------------------------------------------------------------
def main():

    import sys
    import argparse

    parser = argparse.ArgumentParser(description='''Plot 2D data files''')

    # List of input data files

    parser.add_argument('datafiles', type=str, nargs='+',
                        help='Data files to plot')

    # Plot type options

    group = parser.add_argument_group('Plot Options',
                                      '''Options to specify the type of plot to
                                      generate and what data to plot''')

    group.add_argument('--plottype', type=str,
                       choices=['surface', 'surface-ani',
                                'contour', 'contour-ani',
                                'slice', 'point'],
                       default='surface',
                       help='''Set the plot type''')

    group.add_argument('--plotvars', type=int, nargs='+',
                       help='''Variable indices to plot''')

    group.add_argument('--plottimes', type=int, nargs='+',
                       help='''Time indices to plot''')

    # Slice plot options

    group = parser.add_argument_group('Slice Plot Options',
                                      '''Options specific to the slice plot
                                      type''')

    group.add_argument('--slicetype', type=str, default='var',
                       choices=['var', 'time'],
                       help='''The slice plot type''')

    mxgroup = group.add_mutually_exclusive_group()

    mxgroup.add_argument('--yslice', type=int, default=-1,
                         help='''y index to plot''')

    mxgroup.add_argument('--xslice', type=int, default=-1,
                         help='''x index to plot''')

    # Point plot options

    group = parser.add_argument_group('Point Plot Options',
                                      '''Options specific to the point plot
                                      type''')

    group.add_argument('--point', type=int, nargs=2, default=[0, 0],
                       help='''x and y index to plot''')

    # Output options

    group = parser.add_argument_group('Output Options',
                                      '''Options for saving plots''')

    group.add_argument('--save', action='store_true',
                       help='''Save figure to file''')

    group.add_argument('--prefix', type=str,
                       help='''File name prefix for saving the figure''')

    group.add_argument('--merge', action='store_true',
                       help='''Merge PDF output files into a single file''')

    # Figure options

    group = parser.add_argument_group('Figure Options',
                                      '''Options to specify various figure
                                      properties''')

    group.add_argument('--labels', type=str, nargs='+',
                       help='''Data labels for the plot legend''')

    group.add_argument('--title', type=str,
                       help='''Plot title''')

    group.add_argument('--xlabel', type=str,
                       help='''x-axis label''')

    group.add_argument('--ylabel', type=str,
                       help='''y-axis label''')

    group.add_argument('--zlabel', type=str,
                       help='''z-axis label''')

    group.add_argument('--grid', action='store_true',
                       help='''Add grid to plot''')

    # Debugging options

    parser.add_argument('--debug', action='store_true',
                        help='Enable debugging')

    # parse command line args
    args = parser.parse_args()

    # create dictionary with header info
    info = read_header(args)

    # create matrix of subdomain info
    subdomains = read_subdomains(args, info)

    # read data
    time, xvals, yvals, zdata = read_data(args, info, subdomains)

    # setup plot info
    plot_settings(args, info, time, xvals, yvals, zdata)

    # Create plots
    if args.plottype == 'surface':
        plot_surface(args, info, time, xvals, yvals, zdata)

    if args.plottype == 'surface-ani':
        plot_surface_ani(args, info, time, xvals, yvals, zdata)

    if args.plottype == 'contour':
        plot_contour(args, info, time, xvals, yvals, zdata)

    if args.plottype == 'contour-ani':
        plot_contour_ani(args, info, time, xvals, yvals, zdata)

    if args.plottype == 'slice':

        # slice data
        if (args.yslice > -1) and (args.yslice < info['ny']):
            svals = xvals
            sdata = zdata[:, args.yslice, :, :]
            if args.xlabel:
                hlabel = args.xlabel
            else:
                hlabel = 'x'
            suffix = " at y = {:.4f}".format(yvals[args.yslice])
        elif (args.xslice > -1) and (args.xslice < info['nx']):
            svals = yvals
            sdata = zdata[:, :, args.xslice, :]
            if args.ylabel:
                hlabel = args.ylabel
            else:
                hlabel = 'y'
            suffix = " at x = {:.4f}".format(xvals[args.xslice])
        else:
            print("ERROR: invalid xslice or yslice option")
            sys.exit()

        if args.slicetype == 'var':
            plot_slice_vars(args, info, time, svals, sdata, hlabel, suffix)
        else:
            plot_slice_time(args, info, time, svals, sdata, hlabel, suffix)

    if args.plottype == 'point':

        # point data
        pdata = zdata[:, args.point[1], args.point[0], :]
        suffix = " at x = {:.4f}, y = {:.4f}".format(xvals[args.point[0]],
                                                     yvals[args.point[1]])

        plot_point(args, info, time, pdata, suffix)


# -----------------------------------------------------------------------------
# function to print info dictionary
# -----------------------------------------------------------------------------


def print_info(info):

    print("Info dictionary:")
    for k, v in info.items():
        print("  ", k, v)


# -----------------------------------------------------------------------------
# function to read data file header
# -----------------------------------------------------------------------------


def read_header(args):

    import sys
    import shlex
    import numpy as np

    # initialize dictionary of header info variables to None
    keys = ['title', 'varnames', 'nprocs', 'nvar', 'nt', 'nx', 'xl', 'xu',
            'ny', 'yl', 'yu']

    info = dict()
    for k in keys:
        info[k] = None

    # read the first input file and extract info from the header
    with open(args.datafiles[0]) as fn:

        # read the file line by line
        for line in fn:

            # skip empty lines
            if not line.strip():
                continue

            # exit after reading initial comment lines
            if "#" not in line:
                break

            # split line into list
            text = shlex.split(line)

            # plot title
            if "title" in line:
                info['title'] = " ".join(text[2:])
                continue

            # plot variable names
            if "vars" in line:
                info['varnames'] = text[2:]
                continue

            # total number of processes
            if "nprocs" in line:
                info['nprocs'] = int(text[2])
                continue

            # number of variables (at each spatial node)
            if "nvar" in line:
                info['nvar'] = int(text[2])
                continue

            # number of output times
            if "nt" in line:
                info['nt'] = int(text[2])
                continue

            # the global number of nodes in the x-direction, the x lower bound
            # (west) and the x upper bound (east)
            if "nx" in line:
                info['nx'] = int(text[2])
                continue
            if "xl" in line:
                info['xl'] = float(text[2])
                continue
            if "xu" in line:
                info['xu'] = float(text[2])
                continue

            # the global number of nodes in the y-direction, the y lower bound
            # (south) and the y upper bound (north)
            if "ny" in line:
                info['ny'] = int(text[2])
                continue
            if "yl" in line:
                info['yl'] = float(text[2])
                continue
            if "yu" in line:
                info['yu'] = float(text[2])
                continue

    # load data to deduce values and perform sanity checks
    data = np.loadtxt(args.datafiles[0], dtype=np.double)

    # try to fill in missing values
    if info['nvar'] is None:
        info['nvar'] = 1
        print("WARNING: nvar not provided. Using nvar = 1")

    if info['nt'] is None or info['nx'] is None or info['ny'] is None:

        # check if data exists
        if data.ndim != 2:
            print("ERROR: data file is not 2d")
            sys.exit()

        # number of output times
        if info['nt'] is None:
            info['nt'] = np.shape(data)[0]

        # number of spatial nodes
        if info['nx'] is None or info['ny'] is None:
            col = np.shape(data)[1] - 1  # exclude output times
            if info['nx'] is None and info['ny'] is not None:
                info['nx'] = col // (info['nvar'] * info['ny'])
            elif info['nx'] is not None and info['ny'] is None:
                info['ny'] = col // (info['nvar'] * info['nx'])
            else:
                info['nx'] = int(np.sqrt(col // info['nvar']))
                info['ny'] = info['nx']
                print("WARNING: nx and ny not provided. Using nx = ny =",
                      info['nx'])

    # sanity checks
    if info['nt'] != np.shape(data)[0]:
        print("ERROR: nt != nrows", info['nt'], np.shape(data)[0])
        sys.exit()

    if (info['nvar'] * info['nx'] * info['ny']) != (np.shape(data)[1] - 1):
        print("ERROR: nvar * nx * ny != ncols - 1")
        sys.exit()

    # check x-dimension lower and upper bounds
    if info['xl'] is None:
        print("WARNING: xl not provided, using xl = 0")
        info['xl'] = 0.0

    if info['xu'] is None:
        print("WARNING: xu not provided, using xu = 1")
        info['xu'] = 1.0

    # check y-dimension lower and upper bounds
    if info['yl'] is None:
        print("WARNING: yl not provided, using yl = 0")
        info['yl'] = 0.0

    if info['yu'] is None:
        print("WARNING: yu not provided, using yu = 1")
        info['yu'] = 1.0

    # check number of processes
    if info['nprocs'] is None:
        info['nprocs'] = len(args.datafiles)
        print("WARNING: nprocs not provided, using nprocs =", info['nprocs'])

    # check if all the expected input files were provided
    if len(args.datafiles) != info['nprocs']:
        print("ERROR: number of data files (", len(args.datafiles),
              ") does not match number of processes (", info['nprocs'], ")")
        sys.exit()

    if args.debug:
        print('title    = ', info['title'])
        print('varnames = ', info['varnames'])
        print('nprocs   = ', info['nprocs'])
        print('nvar     = ', info['nvar'])
        print('nt       = ', info['nt'])
        print('nx       = ', info['nx'])
        print('xl       = ', info['xl'])
        print('xu       = ', info['xu'])
        print('ny       = ', info['ny'])
        print('yl       = ', info['yl'])
        print('yu       = ', info['yu'])

    return info


# -----------------------------------------------------------------------------
# function to read data file subdomains
# -----------------------------------------------------------------------------


def read_subdomains(args, info):

    import sys
    import shlex
    import numpy as np

    # load subdomain information, store in table
    subdomains = np.zeros((info['nprocs'], 4), dtype=int)

    # get the spatial subdomain owned by each process
    if info['nprocs'] == 1:
        subdomains[0, 0] = 0
        subdomains[0, 1] = info['nx'] - 1
        subdomains[0, 2] = 0
        subdomains[0, 3] = info['ny'] - 1
    else:
        for idx, datafile in enumerate(args.datafiles):

            with open(datafile) as fn:

                # initialize found flags
                found_is = False
                found_ie = False
                found_js = False
                found_je = False

                # read the file line by line
                for line in fn:

                    # skip empty lines
                    if not line.strip():
                        continue

                    # exit after reading initial comment lines
                    if "#" not in line:
                        break

                    # split line into list
                    text = shlex.split(line)

                    # x-direction starting and ending index
                    if "is" in line:
                        subdomains[idx, 0] = int(text[2])
                        found_is = True
                        continue
                    if "ie" in line:
                        subdomains[idx, 1] = int(text[2])
                        found_ie = True
                        continue

                    # y-direction starting and ending index
                    if "js" in line:
                        subdomains[idx, 2] = int(text[2])
                        found_js = True
                        continue
                    if "je" in line:
                        subdomains[idx, 3] = int(text[2])
                        found_je = True
                        continue

                # check if subdomain indices were found
                if not (found_is and found_ie and found_js and found_je):
                    print("ERROR: could not find subdomain indices in",
                          datafile)
                    sys.exit()

    return subdomains


# -----------------------------------------------------------------------------
# read the data files
# -----------------------------------------------------------------------------


def read_data(args, info, subdomains):

    import numpy as np

    # initialize data arrays
    time = np.zeros(info['nt'])
    xvals = np.linspace(info['xl'], info['xu'], info['nx'])
    yvals = np.linspace(info['yl'], info['yu'], info['ny'])
    zdata = np.zeros((info['nt'], info['ny'], info['nx'], info['nvar']))

    # extract data
    for idx, datafile in enumerate(args.datafiles):

        if args.debug:
            print(datafile)

        # load data
        data = np.loadtxt(datafile, dtype=np.double)

        if args.debug:
            print(np.shape(data))

        if np.shape(data)[0] != info['nt']:
            print("WARNING: subdomain", str(idx), "has an incorrect number of"
                  "output times (", np.shape(data)[0], "vs", info['nt'], ")")
            info['nt'] = np.shape(data)[0]

        # x-subdomain indices
        istart = subdomains[idx, 0]
        iend = subdomains[idx, 1]

        # y-subdomain indices
        jstart = subdomains[idx, 2]
        jend = subdomains[idx, 3]

        # local length
        nxl = iend - istart + 1
        nyl = jend - jstart + 1

        if args.debug:
            print(istart, iend, nxl)
            print(jstart, jend, nyl)

        # reshape and save data
        time[:] = data[:, 0]
        for v in range(info['nvar']):
            for i in range(info['nt']):
                zdata[i, jstart:jend+1, istart:iend+1, v] = \
                    np.reshape(data[i, 1+v::info['nvar']], (nyl, nxl))

    return time, xvals, yvals, zdata


# -----------------------------------------------------------------------------
# setup info reused by different plots
# -----------------------------------------------------------------------------


def plot_settings(args, info, time, xvals, yvals, zdata):

    import numpy as np

    # determine extents of plots
    info['zmin'] = np.zeros(info['nvar'])
    info['zmax'] = np.zeros(info['nvar'])

    for v in range(info['nvar']):
        info['zmin'][v] = np.amin(zdata[:, :, :, v])
        info['zmax'][v] = np.amax(zdata[:, :, :, v])

    if args.debug:
        print("z max = ", info['zmax'])
        print("z min = ", info['zmin'])

    # which variables to plot
    if args.plotvars:
        info['pltvars'] = args.plotvars
    else:
        info['pltvars'] = range(info['nvar'])

    # which times to plot
    if args.plottimes:
        info['plttimes'] = args.plottimes
    else:
        info['plttimes'] = range(info['nt'])

    # x-axis label
    if args.xlabel:
        info['xlabel'] = args.xlabel
    else:
        info['xlabel'] = 'x'

    # y-axis label
    if args.ylabel:
        info['ylabel'] = args.ylabel
    else:
        info['ylabel'] = 'y'


# -----------------------------------------------------------------------------
# utility function to combine PDF files
# -----------------------------------------------------------------------------


def merge_pdf(mergefiles, fname):

    import os
    from PyPDF2 import PdfFileMerger

    merger = PdfFileMerger()

    for pdf in mergefiles:
        merger.append(pdf)

    merger.write(fname)
    merger.close()

    for pdf in mergefiles:
        os.remove(pdf)


# -----------------------------------------------------------------------------
# sufrace plot
# -----------------------------------------------------------------------------


def plot_surface(args, info, time, xvals, yvals, zdata):

    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib import cm

    # set x and y meshgrid objects
    X, Y = np.meshgrid(xvals, yvals)

    # generate plots
    for v in info['pltvars']:

        if args.merge:
            mergefiles = list()

        for t in info['plttimes']:

            # create figure and axes
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')

            ax.plot_surface(X, Y, zdata[t, :, :, v], rstride=1, cstride=1,
                            cmap=cm.jet, linewidth=0, antialiased=True,
                            shade=True)

            # set axis limits
            ax.set_xlim([info['xl'], info['xu']])
            ax.set_ylim([info['yl'], info['yu']])
            ax.set_zlim(info['zmin'][v], info['zmax'][v])

            # initial perspective
            ax.view_init(20, -120)

            # add axis labels
            plt.xlabel(info['xlabel'])
            plt.ylabel(info['ylabel'])

            # add z-axis label
            if args.zlabel:
                ax.set_zlabel(args.zlabel)
            elif info['varnames']:
                ax.set_zlabel(info['varnames'][v])
            else:
                ax.set_zlabel('z')

            # add title
            tstr = str(time[t])
            if args.title:
                title = args.title
            elif info['title']:
                title = info['title']
            else:
                title = 'Solution'
            plt.title(title + '\nt = ' + tstr)

            # add grid
            if args.grid:
                plt.grid()

            # save plot to file
            if args.save:
                if args.prefix:
                    fname = args.prefix + '_fig_surface_'
                else:
                    fname = 'fig_surface_'
                if info['varnames']:
                    fname += info['varnames'][v]
                else:
                    fname += 'var_' + repr(v).zfill(3)
                fname += '_t_' + repr(t).zfill(3) + '.pdf'
                plt.savefig(fname, bbox_inches='tight')
                if args.merge:
                    mergefiles.append(fname)
            else:
                plt.show()
            plt.close()

        if args.merge:
            if args.prefix:
                fname = args.prefix + '_fig_surface_'
            else:
                fname = 'fig_surface_'
            if info['varnames']:
                fname += info['varnames'][v]
            else:
                fname += 'var_' + repr(v).zfill(3)
            fname += '.pdf'
            merge_pdf(mergefiles, fname)


# -----------------------------------------------------------------------------
# animated sufrace plot
# -----------------------------------------------------------------------------


def plot_surface_ani(args, info, time, xvals, yvals, zdata):

    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation
    from matplotlib import cm

    def update_plot(frame_number, zarray, v, plot):
        plot[0].remove()
        plot[0] = ax.plot_surface(X, Y, zarray[frame_number, :, :, v],
                                  cmap=cm.jet)

        tstr = str(time[frame_number])
        if args.title:
            title = args.title
        elif info['title']:
            title = info['title']
        else:
            title = 'Solution'
        plt.title(title + '\nt = ' + tstr)

        return plot,

    # set x and y meshgrid objects
    X, Y = np.meshgrid(xvals, yvals)

    # generate plots
    for v in info['pltvars']:

        # create figure and axes
        fig = plt.figure()
        ax = plt.axes(projection='3d')

        plot = [ax.plot_surface(X, Y, zdata[0, :, :, v], rstride=1, cstride=1,
                                cmap=cm.jet, linewidth=0, antialiased=True,
                                shade=True)]

        # set axis limits
        ax.set_xlim([info['xl'], info['xu']])
        ax.set_ylim([info['yl'], info['yu']])
        ax.set_zlim([info['zmin'][v], info['zmax'][v]])

        # initial perspective
        ax.view_init(20, -120)

        # add x-axis label
        if args.xlabel:
            plt.xlabel(args.xlabel)
        else:
            ax.set_xlabel('x')

        # add y-axis label
        if args.ylabel:
            plt.ylabel(args.ylabel)
        else:
            ax.set_ylabel('y')

        # add z-axis label
        if args.zlabel:
            ax.set_zlabel(args.zlabel)
        elif info['varnames']:
            ax.set_zlabel(info['varnames'][v])
        else:
            ax.set_zlabel('z')

        # add grid
        if args.grid:
            plt.grid()

        fps = 2          # frame per sec
        frn = len(time)  # number of frames in the animation

        # create animation
        ani = animation.FuncAnimation(fig, update_plot, frn,
                                      fargs=(zdata, v, plot),
                                      interval=1000/fps)

        # save animation to file
        if args.save:
            if args.prefix:
                fname = args.prefix + '_ani_surface_'
            else:
                fname = 'ani_surface_'
            if info['varnames']:
                fname += info['varnames'][v]
            else:
                fname += 'var_' + repr(v).zfill(3)
            ani.save(fname + '.mp4', dpi=200, fps=fps)
        else:
            plt.show()
        plt.close()


# -----------------------------------------------------------------------------
# contour plot
# -----------------------------------------------------------------------------


def plot_contour(args, info, time, xvals, yvals, zdata):

    import numpy as np
    import matplotlib.pyplot as plt

    # set x and y meshgrid objects
    X, Y = np.meshgrid(xvals, yvals)

    # generate plots
    for v in info['pltvars']:

        levels = np.linspace(info['zmin'][v], info['zmax'][v], 100)
        ticks = np.linspace(info['zmin'][v], info['zmax'][v], 10)

        for t in info['plttimes']:

            # create figure and axes
            fig, ax = plt.subplots()

            cf = ax.contourf(X, Y, zdata[t, :, :, v], levels=levels,
                             cmap="coolwarm", extend="both")
            fig.colorbar(cf, ax=ax, fraction=0.046, pad=0.04, ticks=ticks)

            # set axis limits
            ax.set_xlim([info['xl'], info['xu']])
            ax.set_ylim([info['yl'], info['yu']])

            # add axis labels
            plt.xlabel(info['xlabel'])
            plt.ylabel(info['ylabel'])

            # add title
            tstr = str(time[t])
            if args.title:
                plt.title(args.title + ' at t = ' + tstr)
            elif info['title']:
                plt.title(info['title'] + ' at t = ' + tstr)
            else:
                plt.title('Solution at t = ' + tstr)

            # add grid
            if args.grid:
                plt.grid()

            # save plot to file
            if args.save:
                if args.prefix:
                    fname = args.prefix + '_fig_contour_'
                else:
                    fname = 'fig_contour_'
                if info['varnames']:
                    fname += info['varnames'][v]
                else:
                    fname += 'var_' + repr(v).zfill(3)
                fname += '_t_' + repr(t).zfill(3) + '.pdf'
                plt.savefig(fname, bbox_inches='tight')
            else:
                plt.show()
            plt.close()


# -----------------------------------------------------------------------------
# animated contour plot
# -----------------------------------------------------------------------------


def plot_contour_ani(args, info, time, xvals, yvals, zdata):

    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation

    def update_plot(frame_number, zarray, v, plot):
        plot[0] = ax.contourf(X, Y, zdata[frame_number, :, :, v],
                              levels=levels, cmap="coolwarm", extend="both")

        tstr = str(time[frame_number])
        if args.title:
            title = args.title
        elif info['title']:
            title = info['title']
        else:
            title = 'Solution'
        plt.title(title + '\nt = ' + tstr)

        return plot,

    # set x and y meshgrid objects
    X, Y = np.meshgrid(xvals, yvals)

    # generate plots
    for v in info['pltvars']:

        levels = np.linspace(info['zmin'][v], info['zmax'][v], 100)
        ticks = np.linspace(info['zmin'][v], info['zmax'][v], 10)

        # create figure and axes
        fig, ax = plt.subplots()

        plot = [ax.contourf(X, Y, zdata[0, :, :, v], levels=levels,
                            cmap="coolwarm", extend="both")]
        fig.colorbar(plot[0], ax=ax, fraction=0.046, pad=0.04, ticks=ticks)

        # set axis limits
        ax.set_xlim([info['xl'], info['xu']])
        ax.set_ylim([info['yl'], info['yu']])

        # add axis labels
        plt.xlabel(info['xlabel'])
        plt.ylabel(info['ylabel'])

        # add grid
        if args.grid:
            plt.grid()

        fps = 2          # frame per sec
        frn = len(time)  # number of frames in the animation

        # create animation
        ani = animation.FuncAnimation(fig, update_plot, frn,
                                      fargs=(zdata, v, plot),
                                      interval=1000/fps)

        # save animation to file
        if args.save:
            if args.prefix:
                fname = args.prefix + '_ani_contour_'
            else:
                fname = 'ani_contour_'
            if info['varnames']:
                fname += info['varnames'][v]
            else:
                fname += 'var_' + repr(v).zfill(3)
            ani.save(fname + '.mp4', dpi=200, fps=fps)
        else:
            plt.show()
        plt.close()


# -----------------------------------------------------------------------------
# 1d slice evolution with new figure for each variable
# -----------------------------------------------------------------------------


def plot_slice_vars(args, info, time, svals, sdata, hlabel, suffix):

    import numpy as np
    import matplotlib.pyplot as plt

    # determine extents of slice plot
    smin = np.zeros(info['nvar'])
    smax = np.zeros(info['nvar'])

    for v in range(info['nvar']):
        smin[v] = np.amin(sdata[:, :, v])
        smax[v] = np.amax(sdata[:, :, v])

    if args.debug:
        print(smin)
        print(smax)

    # set labels for the plot legend
    if args.labels:
        label = args.labels
    else:
        label = ["%.2f" % t for t in time]

    # create plot for each variable
    for v in info['pltvars']:

        # create figure and axes
        fig, ax = plt.subplots()

        # add each output time to the plot
        for t in info['plttimes']:
            ax.plot(svals, sdata[t, :, v], label=label[t])

        # set axis limits
        ax.set_xlim([svals[0], svals[-1]])
        ax.set_ylim([0.99 * smin[v], 1.01 * smax[v]])

        # add legend
        ax.legend(title="Times", bbox_to_anchor=(1.02, 1), loc="upper left")

        # add x-axis label
        ax.set_xlabel(hlabel)

        # add y-axis label
        if args.zlabel:
            ax.set_ylabel(args.zlabel)
        else:
            if info['varnames']:
                ax.set_ylabel(info['varnames'][v])
            else:
                ax.set_ylabel('variable ' + repr(v))

        # add title
        if args.title:
            plt.title(args.title + suffix)
        elif info['title']:
            plt.title(info['title'] + suffix)
        else:
            if info['varnames']:
                plt.title("Evolution of " + info['varnames'][v] + suffix)
            else:
                plt.title("Evolution of variable " + repr(v) + suffix)

        # add grid
        if args.grid:
            plt.grid()

        # save plot to file
        if args.save:
            if args.prefix:
                fname = args.prefix + '_fig_slice_'
            else:
                fname = 'fig_slice_'
            if info['varnames']:
                fname += info['varnames'][v]
            else:
                fname += 'var_' + repr(v).zfill(3)
            plt.savefig(fname + '.pdf', bbox_inches='tight')
        else:
            plt.show()
        plt.close()


# -----------------------------------------------------------------------------
# 1d slice evolution with new figure for each time
# -----------------------------------------------------------------------------


def plot_slice_time(args, info, time, svals, sdata, hlabel, suffix):

    import numpy as np
    import matplotlib.pyplot as plt

    # determine extents of slice plot
    smin = np.amin(sdata)
    smax = np.amax(sdata)

    if args.debug:
        print(smin)
        print(smax)

    # set labels for the plot legend
    if args.labels:
        label = args.labels
    elif info['varnames']:
        label = info['varnames']
    else:
        label = [None] * info['nvar']

    # create plot for each variable
    for t in info['plttimes']:

        # create figure and axes
        fig, ax = plt.subplots()

        # add each output time to the plot
        for v in info['pltvars']:
            ax.plot(svals, sdata[t, :, v], label=label[v])

        # set axis limits
        ax.set_xlim([svals[0], svals[-1]])
        ax.set_ylim([0.99 * smin, 1.01 * smax])

        # add legend
        ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left")

        # add x-axis label
        ax.set_xlabel(hlabel)

        # add y-axis label
        if args.zlabel:
            ax.set_ylabel(args.zlabel)

        # add title
        tstr = str(time[t])
        if args.title:
            plt.title(args.title + suffix + ' and t = ' + tstr)
        elif info['title']:
            plt.title(info['title'] + suffix + ' and t = ' + tstr)
        else:
            plt.title("Evolution" + suffix + ' and t = ' + tstr)

        # add grid
        if args.grid:
            plt.grid()

        # save plot to file
        if args.save:
            if args.prefix:
                fname = args.prefix + '_fig_slice_t_'
            else:
                fname = 'fig_slice_t_'
            fname += repr(t).zfill(3) + '.pdf'
            plt.savefig(fname, bbox_inches='tight')
        else:
            plt.show()
        plt.close()


# -----------------------------------------------------------------------------
# point evolution
# -----------------------------------------------------------------------------


def plot_point(args, info, time, pdata, suffix):

    import matplotlib.pyplot as plt

    # set labels for the plot legend
    if args.labels:
        label = args.labels
    elif info['varnames']:
        label = info['varnames']
    else:
        label = [None] * info['nvar']

    # create figure and axes
    fig, ax = plt.subplots()

    # create plot for each variable
    for v in info['pltvars']:
        ax.plot(time, pdata[:, v], label=label[v])

    # add legend
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left")

    # add x-axis label
    ax.set_xlabel("time")

    # add title
    if args.title:
        plt.title(args.title + suffix)
    elif info['title']:
        plt.title(info['title'] + suffix)
    else:
        plt.title("Evolution" + suffix)

    # add grid
    if args.grid:
        plt.grid()

    # save plot to file
    if args.save:
        if args.prefix:
            fname = args.prefix + '_fig_point'
        else:
            fname = 'fig_point'
        plt.savefig(fname + '.pdf', bbox_inches='tight')
    else:
        plt.show()
    plt.close()


# -----------------------------------------------------------------------------
# run the main routine
# -----------------------------------------------------------------------------


if __name__ == '__main__':
    import sys
    sys.exit(main())
