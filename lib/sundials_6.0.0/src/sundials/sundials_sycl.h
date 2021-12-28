/* ---------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * ---------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * ---------------------------------------------------------------------------
 * This header files defines internal utility functions and macros for working
 * with SYCL.
 * ---------------------------------------------------------------------------*/

#include <CL/sycl.hpp>
#include <sundials/sundials_types.h>

#ifndef _SUNDIALS_SYCL_H
#define _SUNDIALS_SYCL_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* Get the maximum work group size (block size) for a queue */
#define SYCL_BLOCKDIM(q)  (q->get_device().get_info<::sycl::info::device::max_work_group_size>())

/* Grid (work group) stride loop */
#define GRID_STRIDE_XLOOP(item, iter, max)              \
  for (sunindextype iter = item.get_global_id(0);       \
       iter < max;                                      \
       iter += item.get_global_range(0))

/* Sycl parallel for loop */
#define SYCL_FOR(q, total, block, item, loop)           \
  q->submit([&](::sycl::handler& h) {                   \
      h.parallel_for(::sycl::nd_range<1>{total,block},  \
                     [=](::sycl::nd_item<1> item)       \
                     { loop }); });

/* Sycl parallel for loop with stream for ouput */
#define SYCL_FOR_DEBUG(q, total, block, item, loop)     \
  q->submit([&](::sycl::handler& h) {                   \
      ::sycl::stream out(1024, 256, h);                 \
      h.parallel_for(::sycl::nd_range<1>{total,block},  \
                     [=](::sycl::nd_item<1> item)       \
                     { loop }); });

/* Sycl parallel for loop with reduction */
#define SYCL_FOR_REDUCE(q, total, block, item, rvar, rop, loop) \
  q->submit([&](::sycl::handler& h) {                           \
      h.parallel_for(::sycl::nd_range<1>{total,block},          \
                     ::sycl::reduction(rvar, rop),              \
                     [=](::sycl::nd_item<1> item, auto& rvar)   \
                     { loop }); });

/* Sycl parallel for loop with reduction and stream for ouput */
#define SYCL_FOR_REDUCE_DEBUG(q, total, block, item, rvar, rop, loop)   \
  q->submit([&](::sycl::handler& h) {                                   \
      ::sycl::stream out(1024, 256, h);                                 \
      h.parallel_for(::sycl::nd_range<1>{total,block},                  \
                     ::sycl::reduction(rvar, rop),                      \
                     [=](::sycl::nd_item<1> item, auto& rvar)           \
                     { loop }); });

#ifdef __cplusplus  /* wrapper to enable C++ usage */
}
#endif

#endif /* _SUNDIALS_SYCL_H */
