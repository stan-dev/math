/*  Copyright 2019 Rene Rivera
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or http://www.boost.org/LICENSE_1_0.txt)
 */

#include "sysinfo.h"
#include "jam.h"
#include "output.h"

#include <thread>

#if defined(OS_MACOSX)
#include <sys/types.h>
#include <sys/sysctl.h>
#endif

#if !defined(OS_NT)
#include <unistd.h>
#endif

#if defined(OS_LINUX)
// Need to define this in case it's not as that's the only way to get the
// sched_* APIs.
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <sched.h>
#endif


b2::system_info::system_info()
{
}

namespace
{
    unsigned int macosx_physicalcpu()
    {
        #if defined(OS_MACOSX)
        int out_hw_ncpu = 0;
        size_t len_hw_ncpu = sizeof(out_hw_ncpu);
        int result = ::sysctlbyname(
            "hw.physicalcpu", &out_hw_ncpu, &len_hw_ncpu, nullptr, 0);
        if (result == 0) return out_hw_ncpu;
        #endif
        return 0;
    }

    unsigned int macosx_logicalcpu()
    {
        #if defined(OS_MACOSX)
        int out_hw_ncpu = 0;
        size_t len_hw_ncpu = sizeof(out_hw_ncpu);
        int result = ::sysctlbyname(
            "hw.logicalcpu", &out_hw_ncpu, &len_hw_ncpu, nullptr, 0);
        if (result == 0) return out_hw_ncpu;
        #endif
        return 0;
    }

    unsigned int sched_affinity_cpu_count()
    {
        #if defined(CPU_COUNT_S)
        ::cpu_set_t cpu_set;
        if (::sched_getaffinity(0, sizeof(cpu_set_t), &cpu_set) == 0)
        {
            return CPU_COUNT_S(sizeof(cpu_set_t), &cpu_set);
        }
        #endif
        return 0;
    }

    unsigned int sysconf_nprocs_configured()
    {
        #if defined(_SC_NPROCESSORS_ONLN)
        return ::sysconf(_SC_NPROCESSORS_CONF);
        #else
        return 0;
        #endif
    }

    unsigned int sysconf_nprocs_online()
    {
        #if defined(_SC_NPROCESSORS_ONLN)
        return ::sysconf(_SC_NPROCESSORS_ONLN);
        #else
        return 0;
        #endif
    }

    unsigned int std_thread_hardware_concurrency()
    {
        return std::thread::hardware_concurrency();
    }
}

unsigned int b2::system_info::cpu_core_count()
{
    if (cpu_core_count_ == 0)
    {
        cpu_thread_count_ = macosx_physicalcpu();
    }
    if (cpu_thread_count_ == 0)
    {
        cpu_thread_count_ = sysconf_nprocs_configured();
    }
    if (cpu_core_count_ <= 0)
    {
        cpu_core_count_ = 1;
    }
    return cpu_core_count_;
}

unsigned int b2::system_info::cpu_thread_count()
{
    if (cpu_thread_count_ == 0)
    {
        cpu_thread_count_ = macosx_logicalcpu();
    }
    if (cpu_thread_count_ == 0)
    {
        cpu_thread_count_ = sched_affinity_cpu_count();
    }
    if (cpu_thread_count_ == 0)
    {
        cpu_thread_count_ = sysconf_nprocs_online();
    }
    if (cpu_thread_count_ == 0)
    {
        cpu_thread_count_ = std_thread_hardware_concurrency();
    }
    if (cpu_thread_count_ == 0)
    {
        cpu_thread_count_ = cpu_core_count();
    }
    return cpu_thread_count_;
}
