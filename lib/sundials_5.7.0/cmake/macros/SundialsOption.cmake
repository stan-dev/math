# ---------------------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# ---------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2021, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------------------
# Provides the macro:
#
#   SUNDIALS_OPTION(<variable> <type> <docstring> <default value>
#                   [SHOW_IF dependencies]
#                   [DEPENDS_ON dependencies]
#                   [DEPNDS_ON_THROW_ERROR])
#
# Within CMake creates a cache variable <variable> and sets it to the value
# <default value> if <variable> is not yet defined. <type> may be any of the
# types valid for CMake's set command (FILEPATH, PATH, STRING, BOOL, INTERNAL).
# <docstring> is a description of the <variable>.
#
# The SHOW_IF option can be used to provide variables which must evaluate
# to true in order to make <variable> a CACHE variable. If not all dependencies
# are met <variable> is made an INTERNAL variable.
#
# The DEPENDS_ON option can be used to provide variable which must evaluate
# to true in order to set <variable>. If the DEPENDS_ON dependencies do not
# all evaluate to true, then <variable> is set to "" and a warning is printed.
# The DEPENDS_ON_THROW_ERROR option will change the warning to be an error.
#
# The OPTIONS macro can be used to provide a list of valid values for the
# variable.
#
# The ADVANCED option can be used to make <variable> an advanced CMake option.
# ---------------------------------------------------------------------------


macro(sundials_option NAME TYPE DOCSTR DEFAULT_VALUE)
  set(options DEPENDS_ON_THROW_ERROR ADVANCED)   # macro options
  set(multiValueArgs OPTIONS SHOW_IF DEPENDS_ON) # macro keyword inputs followed by multiple values

  # parse inputs and create variables sundials_option_<keyword>
  cmake_parse_arguments(sundials_option "${options}" "${oneValueArgs}"
    "${multiValueArgs}" ${ARGN} )

  # check if dependencies for this option have been met
  set(all_show_if_dependencies_met TRUE)
  if(sundials_option_SHOW_IF)
      foreach(_dependency ${sundials_option_SHOW_IF})
        if(NOT ${_dependency})
          set(all_show_if_dependencies_met FALSE)
        endif()
    endforeach()
  endif()

  # check if dependencies for this option have been met
  set(all_depends_on_dependencies_met TRUE)
  if(sundials_option_DEPENDS_ON)
      foreach(_dependency ${sundials_option_DEPENDS_ON})
        if(NOT ${_dependency})
          set(all_depends_on_dependencies_met FALSE)
          list(APPEND depends_on_dependencies_not_met "${_dependency},")
        endif()
    endforeach()
  endif()

  if(all_show_if_dependencies_met AND all_depends_on_dependencies_met)

    if(DEFINED ${NAME})
      # if the variable is already defined, keep its value
      set(${NAME} ${${NAME}} CACHE ${TYPE} ${DOCSTR} FORCE)
    else()
      # if the variable is not already defined, use the default value
      set(${NAME} ${DEFAULT_VALUE} CACHE ${TYPE} ${DOCSTR})
    endif()

    # make the option advanced if necessary
    if(sundials_option_ADVANCED)
      mark_as_advanced(FORCE ${NAME})
    endif()

  elseif(NOT all_depends_on_dependencies_met)

    # don't set the normal variable unless all dependencies were met
    if(DEFINED ${NAME})
      if(${NAME})
        string(CONCAT _warn_msg_string
          "The variable ${NAME} was set to ${${NAME}} "
          "but not all of its dependencies "
          "(${depends_on_dependencies_not_met}) evaluate to TRUE."
        )
        if(sundials_option_DEPENDS_ON_THROW_ERROR)
          print_error("${_warn_msg_string}" "Setting ${NAME} to \"\".")
        else()
          print_warning("${_warn_msg_string}" "Setting ${NAME} to \"\".")
        endif()
      endif()
      set(${NAME} "" CACHE INTERNAL ${DOCSTR} FORCE)
    endif()

  elseif(NOT all_show_if_dependencies_met)

    # hide the normal variable if not all dependencies were met
    if(DEFINED ${NAME})
      set(${NAME} ${${NAME}} CACHE INTERNAL ${DOCSTR} FORCE)
    else()
      set(${NAME} ${DEFAULT_VALUE} CACHE INTERNAL ${DOCSTR})
    endif()

  endif()

  if(sundials_option_OPTIONS)
    if(NOT (${NAME} IN_LIST sundials_option_OPTIONS))
      list(JOIN sundials_option_OPTIONS ", " _options_msg)
      print_error("Value of ${NAME} must be one of ${_options_msg}")
    endif()
    get_property(is_in_cache CACHE ${NAME} PROPERTY TYPE)
    if(is_in_cache)
      set_property(CACHE ${NAME} PROPERTY STRINGS ${sundials_option_OPTIONS})
    endif()
    unset(is_in_cache)
  endif()

  unset(all_show_if_dependencies_met)
  unset(all_depends_on_dependencies_met)
  unset(depends_on_dependencies_not_met)

endmacro()
