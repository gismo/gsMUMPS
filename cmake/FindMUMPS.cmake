## ############################
## Module to find MUMPS package
## ############################

set(MUMPS_FOUND "NO")

message("MUMPS_DIR: ${MUMPS_DIR}")

if(MUMPS_DIR)
  # MUMPS_DIR is specified
  set(MUMPS_INCLUDE_DIR ${MUMPS_DIR}/include)
  set(MUMPS_LIBRARY_DIR ${MUMPS_DIR}/lib)
elseif(DEFINED ENV{MUMPS_DIR})
  # MUMPS_DIR is specified
  set(MUMPS_INCLUDE_DIR $ENV{MUMPS_DIR}/include)
  set(MUMPS_LIBRARY_DIR $ENV{MUMPS_DIR}/lib)
elseif(DEFINED ENV{MUMPS_LIBRARY_DIR})
  # MUMPS_INCLUDE_DIR & LIBRARY_DIR is specified
  set(MUMPS_INCLUDE_DIR $ENV{MUMPS_INCLUDE_DIR})
  set(MUMPS_LIBRARY_DIR $ENV{MUMPS_LIBRARY_DIR})
else()
  # Otherwise look for standard or specific places
  find_path (MUMPS_INCLUDE_DIR
    NAMES zmumps_struc.h
    PATHS
    /opt/MUMPS/include
    /usr/MUMPS/include
    /usr/local/opt/mumps/include
    /usr/include
    /usr/include/MUMPS
    ~/include)

  find_path(MUMPS_LIBRARY_DIR
    NAMES libzmumps.so
    PATHS
    /opt/MUMPS/lib
    /usr/MUMPS/lib
    /usr/local/opt/mumps/lib
    /usr/lib
    /usr/lib64/openmpi/lib
    /usr/lib/x86_64-linux-gnu/
    ~/lib)
endif()

if(MUMPS_INCLUDE_DIR AND MUMPS_LIBRARY_DIR)
  set(MUMPS_FOUND YES)

  find_library(MUMPS_COMMON_LIBRARY
    NAMES mumps_common
    PATHS ${MUMPS_LIBRARY_DIR}
    NO_DEFAULT_PATH)

  find_library(MUMPS_Z_LIBRARY
    NAMES zmumps
    PATHS ${MUMPS_LIBRARY_DIR}
    NO_DEFAULT_PATH)

  find_library(MUMPS_C_LIBRARY
    NAMES cmumps
    PATHS ${MUMPS_LIBRARY_DIR}
    NO_DEFAULT_PATH)

  find_library(MUMPS_D_LIBRARY
    NAMES dmumps
    PATHS ${MUMPS_LIBRARY_DIR}
    NO_DEFAULT_PATH)

  find_library(MUMPS_PORD_LIBRARY
    NAMES pord
    PATHS ${MUMPS_LIBRARY_DIR}
    NO_DEFAULT_PATH)

  set(MUMPS_LIBRARIES ${MUMPS_D_LIBRARY} ${MUMPS_Z_LIBRARY} ${MUMPS_COMMON_LIBRARY})

else()
  if(MUMPS_FIND_REQUIRED)
    message(FATAL_ERROR "MUMPS not found, please set MUMPS_DIR to your MUMPS install directory")
  endif()
endif()

# # FROM: https://github.com/NGSolve/ngsolve/blob/master/cmake/cmake_modules/FindMUMPS.cmake
# find_path (MUMPS_DIR include/zmumps_c.h HINTS ENV MUMPS_DIR DOC "Mumps Directory")
# if(EXISTS ${MUMPS_DIR}/include/zmumps_c.h)
#     set(MUMPS_FOUND YES)
#     set(MUMPS_INCLUDES ${MUMPS_DIR})
#     find_path (MUMPS_INCLUDE_DIR mumps_compat.h HINTS "${MUMPS_DIR}" PATH_SUFFIXES include NO_DEFAULT_PATH)
#     list(APPEND MUMPS_INCLUDES ${MUMPS_INCLUDE_DIR})
#     find_library(LIB_MUMPS_COMMON mumps_common PATHS ${MUMPS_DIR}/lib)
#     find_library(LIB_MUMPS_D dmumps PATHS ${MUMPS_DIR}/lib)
#     find_library(LIB_MUMPS_Z zmumps PATHS ${MUMPS_DIR}/lib)
#     find_library(LIB_PORD pord PATHS ${MUMPS_DIR}/lib)
#     # find_library(LIB_PARMETIS parmetis HINTS ${PARMETIS_DIR}/lib REQUIRED)
#     # find_library(LIB_METIS metis HINTS ${PARMETIS_DIR}/lib REQUIRED)

#     # if (NOT USE_MKL)
#     #     find_library(LIB_SCALAPACK scalapack HINTS ${SCALAPACK_DIR}/lib REQUIRED)
#     # endif()

#     set(MUMPS_LIBRARIES ${LIB_MUMPS_D} ${LIB_MUMPS_Z} ${LIB_MUMPS_COMMON})

#     if (LIB_PORD)
#        list(APPEND MUMPS_LIBRARIES ${LIB_PORD})
#     endif()

# endif(EXISTS ${MUMPS_DIR}/include/zmumps_c.h)
# include(FindPackageHandleStandardArgs)
# find_package_handle_standard_args(MUMPS DEFAULT_MSG MUMPS_LIBRARIES MUMPS_INCLUDES)
