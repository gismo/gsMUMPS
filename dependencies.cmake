# dependencies.cmake for gsMUMPS module
# This file is automatically included by the main CMakeLists.txt
# when gsMUMPS is in GISMO_OPTIONAL

# Copy FindMUMPS.cmake to main cmake directory
if(EXISTS "${PROJECT_SOURCE_DIR}/optional/gsMUMPS/cmake/FindMUMPS.cmake")
  execute_process(COMMAND ${CMAKE_COMMAND} -E copy
    "${PROJECT_SOURCE_DIR}/optional/gsMUMPS/cmake/FindMUMPS.cmake"
    "${PROJECT_SOURCE_DIR}/cmake/")
endif()

# Find MUMPS
find_package(MUMPS REQUIRED)

if(MUMPS_FOUND)
  message(STATUS "MUMPS found: ${MUMPS_INCLUDE_DIR}")
  message(STATUS "MUMPS libraries: ${MUMPS_LIBRARIES}")

  # Add MUMPS include directories to global includes
  set(GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} ${MUMPS_INCLUDE_DIR}
    CACHE INTERNAL "${PROJECT_NAME} include directories")

  # Add MUMPS libraries to global linker
  set(gismo_LINKER ${gismo_LINKER} ${MUMPS_LIBRARIES}
    CACHE INTERNAL "${PROJECT_NAME} extra linker objects")

else()
  message(WARNING "MUMPS not found. gsMUMPS module will be disabled.")
endif()