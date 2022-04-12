/*
 * -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 */

/* NOTE: SUNDIALS_HOST_DEVICE and SUNDIALS_DEVICE_INLINE must be defined
         before including this file */

#include <limits>

namespace sundials
{
namespace reductions
{
namespace impl
{

template<typename Arg1, typename Arg2, typename Result>
struct BinaryOperator
{
  using first_arg_type = Arg1;
  using second_arg_type = Arg2;
  using result_arg_type = Result;
};

template<typename Ret, typename Arg1 = Ret, typename Arg2 = Arg1>
struct plus : public BinaryOperator<Arg1, Arg2, Ret>
{
  SUNDIALS_HOST_DEVICE constexpr Ret operator()(const Arg1& lhs,
                                                const Arg2& rhs) const
  {
    return Ret{lhs} + rhs;
  }

  static SUNDIALS_HOST_DEVICE SUNDIALS_DEVICE_INLINE constexpr Ret identity()
  {
    return Ret{0};
  }
};

template<typename Ret, typename Arg1 = Ret, typename Arg2 = Arg1>
struct maximum : public BinaryOperator<Arg1, Arg2, Ret>
{
  SUNDIALS_HOST_DEVICE constexpr Ret operator()(const Arg1& lhs,
                                                const Arg2& rhs) const
  {
    return (lhs >= rhs) ? lhs : rhs;
  }

  static SUNDIALS_HOST_DEVICE SUNDIALS_DEVICE_INLINE constexpr Ret identity()
  {
    return std::numeric_limits<Ret>::min();
  }
};

template<typename Ret, typename Arg1 = Ret, typename Arg2 = Arg1>
struct minimum : public BinaryOperator<Arg1, Arg2, Ret>
{
  SUNDIALS_HOST_DEVICE constexpr Ret operator()(const Arg1& lhs,
                                                const Arg2& rhs) const
  {
    return (rhs < lhs) ? rhs : lhs;
  }

  static SUNDIALS_HOST_DEVICE SUNDIALS_DEVICE_INLINE constexpr Ret identity()
  {
    return std::numeric_limits<Ret>::max();
  }
};

} // impl
} // reductions
} // sundials
