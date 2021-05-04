#!/bin/sh

#~ Copyright 2002-2019 Rene Rivera.
#~ Distributed under the Boost Software License, Version 1.0.
#~ (See accompanying file LICENSE_1_0.txt or copy at
#~ http://www.boost.org/LICENSE_1_0.txt)

# Reset the toolset.
B2_TOOLSET=
B2_OS=

# Run a command, and echo before doing so. Also checks the exit status and quits
# if there was an error.
echo_run ()
{
    echo "> $@"
    $@
    r=$?
    if test $r -ne 0 ; then
        exit $r
    fi
}

# Print an error message, and exit with a status of 1.
error_exit ()
{
    echo "
${@}

You can specify the toolset as the argument, i.e.:
    ./build.sh gcc

Toolsets supported by this script are:
    acc, clang, como, gcc, intel-darwin, intel-linux, kcc, kylix, mipspro,
    pathscale, pgi, qcc, sun, sunpro, tru64cxx, vacpp

For any toolset you can override the path to the compiler with the CXX
environment variable. You can also use additional flags for the compiler
with the CXXFLAGS environment variable.

A special toolset; cxx, is available which is used as a fallback when a more
specific toolset is not found and the cxx command is detected. The 'cxx'
toolset will use the CXX, CXXFLAGS, and LIBS environment variables, if present.

Similarly, the cross-cxx toolset is available for cross-compiling by using the
BUILD_CXX, BUILD_CXXFLAGS, and BUILD_LDFLAGS environment variables to compile
binaries that will be executed on the build system. This allows CXX etc. to be
set for cross-compilers to be propagated to subprocesses.
" 1>&2
    exit 1
}

# Check that a command is in the PATH.
test_path ()
{
    if `command -v command 1>/dev/null 2>/dev/null`; then
        command -v $1 1>/dev/null 2>/dev/null
    else
        hash $1 1>/dev/null 2>/dev/null
    fi
}

# Check that the OS name, as returned by "uname", is as given.
test_uname ()
{
    if test_path uname; then
        test `uname` = $*
    fi
}

# Check that the given command runs.
test_exec ()
{
    "$*" 1>/dev/null 2>/dev/null
}

# Check that the compiler can do C++11.
test_cxx11 ()
{
    if ! test $NO_CXX11_CHECK ; then
        case $1 in
            gcc) ( ${CXX:=g++} -x c++ -std=c++11 check_cxx11.cpp && rm -f a.out ) 1>/dev/null 2>/dev/null ;;
            intel-darwin) ( ${CXX:=icc} -xc++ check_cxx11.cpp && rm -f a.out ) 1>/dev/null 2>/dev/null ;;
            intel-linux) ( ${CXX:=icc} -xc++ check_cxx11.cpp && rm -f a.out ) 1>/dev/null 2>/dev/null ;;
            vacpp) ( ${CXX:=xlC_r} check_cxx11.cpp && rm -f a.out ) 1>/dev/null 2>/dev/null ;;
            xlcpp) ( ${CXX:=xlC_r} check_cxx11.cpp && rm -f a.out ) 1>/dev/null 2>/dev/null ;;
            como) ( ${CXX:=como} check_cxx11.cpp && rm -f a.out ) 1>/dev/null 2>/dev/null ;;
            kcc) ( ${CXX:=KCC} check_cxx11.cpp && rm -f a.out ) 1>/dev/null 2>/dev/null ;;
            kylix) ( ${CXX:=bc++} -tC -q check_cxx11.cpp && rm -f a.out ) 1>/dev/null 2>/dev/null ;;
            mipspro) ( ${CXX:=CC} -FE:template_in_elf_section -ptused check_cxx11.cpp && rm -f a.out ) 1>/dev/null 2>/dev/null ;;
            pathscale) ( ${CXX:=pathCC} check_cxx11.cpp && rm -f a.out ) 1>/dev/null 2>/dev/null ;;
            pgi) ( ${CXX:=pgc++} -std=c++11 check_cxx11.cpp && rm -f a.out ) 1>/dev/null 2>/dev/null ;;
            sun*) ( ${CXX:=CC} -std=c++11 check_cxx11.cpp && rm -f a.out ) 1>/dev/null 2>/dev/null ;;
            clang*) ( ${CXX:=clang++} -x c++ -std=c++11 check_cxx11.cpp && rm -f a.out ) 1>/dev/null 2>/dev/null ;;
            tru64cxx) ( ${CXX:=cc} check_cxx11.cpp && rm -f a.out ) 1>/dev/null 2>/dev/null ;;
            acc) ( ${CXX:=aCC} -AA check_cxx11.cpp && rm -f a.out ) 1>/dev/null 2>/dev/null ;;
            qcc) ( ${CXX:=QCC} check_cxx11.cpp && rm -f a.out ) 1>/dev/null 2>/dev/null ;;
            cxx) ( ${CXX:=cxx} check_cxx11.cpp && rm -f a.out ) 1>/dev/null 2>/dev/null ;;
            cross-cxx) ( ${CXX:=cxx} check_cxx11.cpp && rm -f a.out ) 1>/dev/null 2>/dev/null ;;
            *) test "0" = "1" ;;
        esac
    else
        test $NO_CXX11_CHECK
    fi
}

# Try and guess the toolset to bootstrap the build with...
guess_toolset ()
{
    if test_uname Darwin && test_cxx11 clang ; then B2_TOOLSET=clang
    elif test_uname IRIX && test_cxx11 mipspro ; then B2_TOOLSET=mipspro
    elif test_uname IRIX64 && test_cxx11 mipspro ; then B2_TOOLSET=mipspro
    elif test_uname OSF1 && test_cxx11 tru64cxx ; then B2_TOOLSET=tru64cxx
    elif test_uname QNX && test_path QCC && test_cxx11 qcc ; then B2_TOOLSET=qcc
    elif test_uname Linux && test_path xlC_r ; then
       if /usr/bin/lscpu | grep Byte | grep Little > /dev/null 2>&1 ; then
          # Little endian linux
          B2_TOOLSET=xlcpp
       else
          #Big endian linux
          B2_TOOLSET=vacpp
       fi
    elif test_uname AIX && test_path xlC_r && test_cxx11 vacpp ; then B2_TOOLSET=vacpp
    elif test_uname FreeBSD && test_path freebsd-version && test_path clang++ && test_cxx11 clang ; then B2_TOOLSET=clang
    elif test_path g++ && test_cxx11 gcc ; then B2_TOOLSET=gcc
    elif test_path clang++ && test_cxx11 clang ; then B2_TOOLSET=clang
    elif test_path icc && test_cxx11 intel-linux ; then B2_TOOLSET=intel-linux
    elif test -r /opt/intel/cc/9.0/bin/iccvars.sh && test_cxx11 intel-linux ; then
        B2_TOOLSET=intel-linux
        B2_TOOLSET_ROOT=/opt/intel/cc/9.0
    elif test -r /opt/intel_cc_80/bin/iccvars.sh && test_cxx11 intel-linux ; then
        B2_TOOLSET=intel-linux
        B2_TOOLSET_ROOT=/opt/intel_cc_80
    elif test -r /opt/intel/compiler70/ia32/bin/iccvars.sh && test_cxx11 intel-linux ; then
        B2_TOOLSET=intel-linux
        B2_TOOLSET_ROOT=/opt/intel/compiler70/ia32/
    elif test -r /opt/intel/compiler60/ia32/bin/iccvars.sh && test_cxx11 intel-linux ; then
        B2_TOOLSET=intel-linux
        B2_TOOLSET_ROOT=/opt/intel/compiler60/ia32/
    elif test -r /opt/intel/compiler50/ia32/bin/iccvars.sh && test_cxx11 intel-linux ; then
        B2_TOOLSET=intel-linux
        B2_TOOLSET_ROOT=/opt/intel/compiler50/ia32/
    elif test_path pgc++ && test_cxx11 pgi ; then B2_TOOLSET=pgi
    elif test_path pathCC && test_cxx11 pathscale ; then B2_TOOLSET=pathscale
    elif test_path como && test_cxx11 como ; then B2_TOOLSET=como
    elif test_path KCC && test_cxx11 kcc ; then B2_TOOLSET=kcc
    elif test_path bc++ && test_cxx11 kylix ; then B2_TOOLSET=kylix
    elif test_path aCC && test_cxx11 acc ; then B2_TOOLSET=acc
    elif test_uname HP-UX ; then B2_TOOLSET=acc
    elif test -r /opt/SUNWspro/bin/cc && test_cxx11 sunpro ; then
        B2_TOOLSET=sunpro
        B2_TOOLSET_ROOT=/opt/SUNWspro/
    # Test for some common compile command as the default fallback.
    elif test_path $CXX ; then B2_TOOLSET=cxx
    elif test_path cxx ; then
        B2_TOOLSET=cxx
        CXX=cxx
    elif test_path cpp ; then
        B2_TOOLSET=cxx
        CXX=cpp
    elif test_path CC ; then
        B2_TOOLSET=cxx
        CXX=CC
    fi
    if test "$B2_TOOLSET" = "" ; then
        error_exit "Could not find a suitable toolset."
    fi
}

check_debug_build ()
{
    while test $# -gt 0
    do
        case "$1" in
            --debug) return 0 ;;
        esac
        shift
    done
    return 1
}

# The one option we support in the invocation
# is the name of the toolset to force building
# with.
case "$1" in
    --guess-toolset) NO_CXX11_CHECK=1 ; guess_toolset ; echo "$B2_TOOLSET" ; exit 1 ;;
    -*) guess_toolset ;;
    ?*) B2_TOOLSET=$1 ; shift ;;
    *) guess_toolset ;;
esac

# We need a C++11 compiler. Check here and given some feedback about it.
if ! test_cxx11 $B2_TOOLSET ; then
    error_exit "
A C++11 capable compiler is required for building the B2 engine.
Toolset '$B2_TOOLSET' does not appear to support C++11.

** Note, the C++11 capable compiler is _only_ required for building the B2
** engine. The B2 build system allows for using any C++ level and any other
** supported language and resource in your projects.
"
fi

case $B2_TOOLSET in

    gcc)
        CXX=${CXX:=g++}
        CXX_VERSION_OPT=${CXX_VERSION_OPT:=--version}
        # Check whether it's MinGW GCC, which has Windows headers and none of POSIX ones.
        machine=$(${CXX} -dumpmachine 2>/dev/null)
        if test $? -ne 0 ; then
            echo "B2_TOOLSET is gcc, but the 'gcc' command cannot be executed."
            echo "Make sure 'gcc' is in PATH, or use a different toolset."
            exit 1
        fi
        case $machine in
        *mingw*)
        # MinGW insists that its bin directory be in PATH.
        if test -r ${B2_TOOLSET_ROOT}bin/gcc ; then
            export PATH=${B2_TOOLSET_ROOT}bin:$PATH
        fi
        B2_CXX="${CXX} -x c++ -std=c++11"
        B2_CXXFLAGS_RELEASE="-O2 -s"
        B2_CXXFLAGS_DEBUG="-O0 -g"
        B2_OS="NT"
        ;;

        *cygwin*)
        B2_CXX="${CXX} -x c++ -std=gnu++11"
        B2_CXXFLAGS_RELEASE="-O2 -s"
        B2_CXXFLAGS_DEBUG="-O0 -g"
        ;;

        *)
        B2_CXX="${CXX} -x c++ -std=c++11"
        B2_CXXFLAGS_RELEASE="-O2 -s"
        B2_CXXFLAGS_DEBUG="-O0 -g"
        esac
    ;;

    intel-darwin)
        CXX=${CXX:=icc}
        CXX_VERSION_OPT=${CXX_VERSION_OPT:=--version}
        B2_CXX="${CXX} -xc++"
        B2_CXXFLAGS_RELEASE="-O3 -s"
        B2_CXXFLAGS_DEBUG="-O0 -g -p"
    ;;

    intel-linux)
        CXX=${CXX:=icc}
        CXX_VERSION_OPT=${CXX_VERSION_OPT:=--version}
        test_path ${CXX} >/dev/null 2>&1
        if test $? ; then
            echo "Found ${CXX} in environment"
            B2_TOOLSET_ROOT=`echo ${CXX}| sed -e 's/bin.*\/icc//'`
            # probably the most widespread
            ARCH=intel64
        else
            echo "No intel compiler in current path"
            echo "Look in a few old place for legacy reason"
            if test -r /opt/intel/cc/9.0/bin/iccvars.sh ; then
                B2_TOOLSET_ROOT=/opt/intel/cc/9.0/
            elif test -r /opt/intel_cc_80/bin/iccvars.sh ; then
                B2_TOOLSET_ROOT=/opt/intel_cc_80/
            elif test -r /opt/intel/compiler70/ia32/bin/iccvars.sh ; then
                B2_TOOLSET_ROOT=/opt/intel/compiler70/ia32/
            elif test -r /opt/intel/compiler60/ia32/bin/iccvars.sh ; then
                B2_TOOLSET_ROOT=/opt/intel/compiler60/ia32/
            elif test -r /opt/intel/compiler50/ia32/bin/iccvars.sh ; then
                B2_TOOLSET_ROOT=/opt/intel/compiler50/ia32/
            fi
        fi
        if test -r ${B2_TOOLSET_ROOT}bin/iccvars.sh ; then
            # iccvars does not change LD_RUN_PATH. We adjust LD_RUN_PATH here in
            # order not to have to rely on ld.so.conf knowing the icc library
            # directory. We do this before running iccvars.sh in order to allow a
            # user to add modifications to LD_RUN_PATH in iccvars.sh.
            if test -z "${LD_RUN_PATH}"; then
                LD_RUN_PATH="${B2_TOOLSET_ROOT}lib"
            else
                LD_RUN_PATH="${B2_TOOLSET_ROOT}lib:${LD_RUN_PATH}"
            fi
            export LD_RUN_PATH
            . ${B2_TOOLSET_ROOT}bin/iccvars.sh $ARCH
        fi
        B2_CXX="${CXX} -xc++"
        B2_CXXFLAGS_RELEASE="-O3 -s"
        B2_CXXFLAGS_DEBUG="-O0 -g -p"
    ;;

    vacpp)
        CXX=${CXX:=xlC_r}
        CXX_VERSION_OPT=${CXX_VERSION_OPT:=-qversion}
        B2_CXX="${CXX}"
        B2_CXXFLAGS_RELEASE="-O3 -s -qstrict -qinline"
        B2_CXXFLAGS_DEBUG="-g -qNOOPTimize -qnoinline -pg"
    ;;

    xlcpp)
        CXX=${CXX:=xlC_r}
        CXX_VERSION_OPT=${CXX_VERSION_OPT:=-qversion}
        B2_CXX="${CXX}"
        B2_CXXFLAGS_RELEASE="-s -O3 -qstrict -qinline"
        B2_CXXFLAGS_DEBUG="-g -qNOOPTimize -qnoinline -pg"
    ;;

    como)
        CXX=${CXX:=como}
        CXX_VERSION_OPT=${CXX_VERSION_OPT:=--version}
        B2_CXX="${CXX}"
        B2_CXXFLAGS_RELEASE="-O3 --inlining"
        B2_CXXFLAGS_DEBUG="-O0 -g --no_inlining --long_long"
    ;;

    kcc)
        CXX=${CXX:=KCC}
        CXX_VERSION_OPT=${CXX_VERSION_OPT:=--version}
        B2_CXX="KCC"
        B2_CXXFLAGS_RELEASE="+K2 -s"
        B2_CXXFLAGS_DEBUG="+K0 -g"
    ;;

    kylix)
        CXX=${CXX:=bc++}
        CXX_VERSION_OPT=${CXX_VERSION_OPT:=--version}
        B2_CXX="bc++ -tC -q"
        B2_CXXFLAGS_RELEASE="-O2 -vi -w-inl -s"
        B2_CXXFLAGS_DEBUG="-Od -v -vi-"
    ;;

    mipspro)
        CXX=${CXX:=CC}
        CXX_VERSION_OPT=${CXX_VERSION_OPT:=--version}
        B2_CXX="${CXX} -FE:template_in_elf_section -ptused"
        B2_CXXFLAGS_RELEASE="-Ofast -g0 \"-INLINE:none\" -s"
        B2_CXXFLAGS_DEBUG="-O0 -INLINE -g"
    ;;

    pathscale)
        CXX=${CXX:=pathCC}
        CXX_VERSION_OPT=${CXX_VERSION_OPT:=--version}
        B2_CXX="${CXX}"
        B2_CXXFLAGS_RELEASE="-O3 -inline -s"
        B2_CXXFLAGS_DEBUG="-O0 -noinline -ggdb"
    ;;

    pgi)
        CXX=${CXX:=pgc++}
        CXX_VERSION_OPT=${CXX_VERSION_OPT:=--version}
        B2_CXX="${CXX} -std=c++11"
        B2_CXXFLAGS_RELEASE="-fast -s"
        B2_CXXFLAGS_DEBUG="-O0 -gopt"
    ;;

    sun*)
        CXX=${CXX:=CC}
        CXX_VERSION_OPT=${CXX_VERSION_OPT:=-V}
        if test -z "${B2_TOOLSET_ROOT}" -a -r /opt/SUNWspro/bin/CC ; then
            B2_TOOLSET_ROOT=/opt/SUNWspro/
        fi
        if test -r "${B2_TOOLSET_ROOT}/bin/CC" ; then
            PATH=${B2_TOOLSET_ROOT}bin:${PATH}
            export PATH
        fi
        B2_CXX="${CXX} -std=c++11"
        B2_CXXFLAGS_RELEASE="-xO4 -s"
        B2_CXXFLAGS_DEBUG="-g"
    ;;

    clang*)
        CXX=${CXX:=clang++}
        CXX_VERSION_OPT=${CXX_VERSION_OPT:=--version}
        B2_CXX="${CXX} -x c++ -std=c++11"
        B2_TOOLSET=clang
        B2_CXXFLAGS_RELEASE="-O3 -s"
        B2_CXXFLAGS_DEBUG="-O0 -fno-inline -g"
    ;;

    tru64cxx)
        CXX_VERSION_OPT=${CXX_VERSION_OPT:=--version}
        B2_CXX="cc"
        B2_CXXFLAGS_RELEASE="-O5 -inline speed -s"
        B2_CXXFLAGS_DEBUG="-O0 -pg -g"
    ;;

    acc)
        CXX=${CXX:=aCC}
        CXX_VERSION_OPT=${CXX_VERSION_OPT:=--version}
        B2_CXX="${CXX} -AA"
        B2_CXXFLAGS_RELEASE="-O3 -s"
        B2_CXXFLAGS_DEBUG="+d -g"
    ;;

    qcc)
        CXX=${CXX:=QCC}
        CXX_VERSION_OPT=${CXX_VERSION_OPT:=--version}
        B2_CXX="${CXX}"
        B2_CXXFLAGS_RELEASE="-O3 -Wc,-finline-functions"
        B2_CXXFLAGS_DEBUG="O0 -Wc,-fno-inline -gstabs+"
    ;;

    cxx)
        CXX=${CXX:=cxx}
        CXX_VERSION_OPT=${CXX_VERSION_OPT:=--version}
        B2_CXX="${CXX}"
    ;;

    cross-cxx)
        CXX=${BUILD_CXX:=cxx}
        CXXFLAGS=${BUILD_CXXFLAGS}
        CXX_VERSION_OPT=${CXX_VERSION_OPT:=--version}
        B2_CXX="${CXX}"
    ;;

    *)
        error_exit "Unknown toolset: $B2_TOOLSET"
    ;;
esac

echo "
###
###
### Using '$B2_TOOLSET' toolset.
###
###
"
echo_run ${CXX} ${CXX_VERSION_OPT}
echo "
###
###
"
B2_SOURCES="\
 builtins.cpp \
 class.cpp \
 command.cpp \
 compile.cpp \
 constants.cpp \
 cwd.cpp \
 debug.cpp \
 debugger.cpp \
 execcmd.cpp \
 filesys.cpp \
 frames.cpp \
 function.cpp \
 glob.cpp\
 hash.cpp \
 hcache.cpp \
 hdrmacro.cpp \
 headers.cpp \
 jam.cpp \
 jamgram.cpp \
 lists.cpp \
 make.cpp \
 make1.cpp \
 md5.cpp \
 mem.cpp \
 modules.cpp \
 native.cpp \
 object.cpp \
 option.cpp \
 output.cpp \
 parse.cpp \
 pathsys.cpp \
 regexp.cpp \
 rules.cpp \
 scan.cpp \
 search.cpp \
 jam_strings.cpp \
 startup.cpp \
 subst.cpp \
 sysinfo.cpp \
 timestamp.cpp \
 variable.cpp \
 w32_getreg.cpp \
 modules/order.cpp \
 modules/path.cpp \
 modules/property-set.cpp \
 modules/regex.cpp \
 modules/sequence.cpp \
 modules/set.cpp \
 "
case $B2_OS in
    NT)
    B2_SOURCES="${B2_SOURCES} execnt.cpp filent.cpp pathnt.cpp"
    ;;

    *)
    B2_SOURCES="${B2_SOURCES} execunix.cpp fileunix.cpp pathunix.cpp"
    ;;
esac

if check_debug_build "$@" ; then B2_CXXFLAGS="${B2_CXXFLAGS_DEBUG}"
else B2_CXXFLAGS="${B2_CXXFLAGS_RELEASE} -DNDEBUG"
fi
echo_run ${B2_CXX} ${CXXFLAGS} ${B2_CXXFLAGS} ${B2_SOURCES} -o b2
echo_run cp b2 bjam
