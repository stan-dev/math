find_package(Doxygen)

if (DOXYGEN_FOUND)

configure_file(Doxyfile.in Doxyfile)
add_custom_target(
  docs
  COMMAND
  ${DOXYGEN_EXECUTABLE}
  WORKING_DIRECTORY
  ${CMAKE_CURRENT_BINARY_DIR}
  VERBATIM)

endif (DOXYGEN_FOUND)
