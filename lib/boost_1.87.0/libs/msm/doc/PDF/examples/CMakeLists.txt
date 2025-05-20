cmake_minimum_required(VERSION 3.17)

# project-specific includes
include_directories(.)
include_directories(../../../include)


IF (WIN32)
include_directories(C:/dev/boost_1_84_0)
set(LIBS_BOOST C:/dev/boost_1_84_0/stage/lib)
ELSE()
include_directories(/home/bob/.vs/boost_1_84_0)
set(LIBS_BOOST /home/bob/.vs/boost_1_84_0/stage/lib)
ENDIF()


# set the project name
project(msm_examples)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if (MSVC)
  add_compile_options(/bigobj)
else ()

endif ()


# add the executable
add_executable( msm_examples
                #back
                ../../../include/boost/msm/back/state_machine.hpp
                ../../../include/boost/msm/back/metafunctions.hpp
                ../../../include/boost/msm/back/dispatch_table.hpp
                #back11
                ../../../include/boost/msm/back11/state_machine.hpp
                ../../../include/boost/msm/back11/metafunctions.hpp
                ../../../include/boost/msm/back11/dispatch_table.hpp
                #front
                ../../../include/boost/msm/front/state_machine_def.hpp
                ../../../include/boost/msm/front/states.hpp
                ../../../include/boost/msm/front/functor_row.hpp
                ../../../include/boost/msm/front/euml/euml_typeof.hpp
                ../../../include/boost/msm/front/euml/state_grammar.hpp
                ../../../include/boost/msm/front/euml/common.hpp
                #sources
                ${PROJECT_SOURCE_DIR}/History.cpp
              )

target_include_directories(msm_examples PUBLIC ${INCLUDE_ASYNCHRONOUS})


IF (WIN32)
target_link_directories(msm_examples PUBLIC ${LIBS_BOOST})
ELSE()
target_link_libraries(msm_examples rt)
target_link_libraries(msm_examples stdc++fs)
target_link_libraries(msm_examples pthread)
target_link_libraries(msm_examples stdc++)
target_link_libraries(msm_examples m)

target_link_libraries(msm_examples ${LIBS_BOOST}/libboost_unit_test_framework.a)
target_link_libraries(msm_examples ${LIBS_BOOST}/libboost_serialization.a)
ENDIF()

