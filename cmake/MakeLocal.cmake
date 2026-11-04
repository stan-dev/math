option(STAN_USE_MAKE_LOCAL "Read supported settings from make/local" ON)
set(STAN_MAKE_LOCAL "${CMAKE_SOURCE_DIR}/make/local" CACHE FILEPATH "Local Make settings file")

macro(stan_cache_default name value type description)
  if(NOT DEFINED ${name})
    set(${name} "${value}" CACHE ${type} "${description}")
  endif()
endmacro()

if(STAN_USE_MAKE_LOCAL)
  file(GLOB local_settings CONFIGURE_DEPENDS "${STAN_MAKE_LOCAL}")
  if(local_settings)
    find_package(Python3 REQUIRED COMPONENTS Interpreter)
    find_program(STAN_SETTINGS_MAKE NAMES gmake make mingw32-make)
    if(NOT STAN_SETTINGS_MAKE)
      message(FATAL_ERROR "Reading make/local requires GNU Make. Use STAN_USE_MAKE_LOCAL=OFF to ignore make/local.")
    endif()
    execute_process(
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_LIST_DIR}/make_local.py"
        --make "${STAN_SETTINGS_MAKE}" --source "${CMAKE_SOURCE_DIR}"
        --local "${STAN_MAKE_LOCAL}" --output "${CMAKE_BINARY_DIR}/make-local.cmake"
      RESULT_VARIABLE settings_status
      ERROR_VARIABLE settings_error)
    if(NOT settings_status EQUAL 0)
      message(FATAL_ERROR "${settings_error}")
    endif()
    include("${CMAKE_BINARY_DIR}/make-local.cmake")
    set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
      "${CMAKE_CURRENT_LIST_DIR}/make_local.py")
    message(STATUS "Loaded Make settings from ${STAN_MAKE_LOCAL}")
  endif()
endif()

function(stan_apply_make_local_flags)
  foreach(language C CXX)
    separate_arguments(flags NATIVE_COMMAND "${STAN_MAKE_LOCAL_${language}_FLAGS}")
    foreach(flag IN LISTS flags)
      add_compile_options("$<$<COMPILE_LANGUAGE:${language}>:${flag}>")
    endforeach()
  endforeach()
  separate_arguments(flags NATIVE_COMMAND "${STAN_MAKE_LOCAL_CPP_FLAGS}")
  add_compile_options(${flags})
  separate_arguments(flags NATIVE_COMMAND "${STAN_MAKE_LOCAL_LINK_FLAGS}")
  add_link_options(${flags})
  separate_arguments(libraries NATIVE_COMMAND "${STAN_MAKE_LOCAL_LIBRARIES_FLAGS}")
  link_libraries(${libraries})
endfunction()
