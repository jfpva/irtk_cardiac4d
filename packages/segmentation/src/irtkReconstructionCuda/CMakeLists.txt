project(irtkReconstructionCuda)

cmake_minimum_required(VERSION 2.6)
#cmake_policy(SET CMP0010 OLD)

#set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)
#set(LIBRARY_OUTPUT_PATH  ${CMAKE_SOURCE_DIR}/lib)

#if IRTK_SOURCE_DIR is set, we are building within the IRTK source tree
if(IRTK_SOURCE_DIR)
SET(BUILD_OUT_OF_IRTK_TREE FALSE)
message("building within IRTK source tree")
else(IRTK_SOURCE_DIR)
SET(BUILD_OUT_OF_IRTK_TREE TRUE)
message("building out of IRTK source tree")
endif(IRTK_SOURCE_DIR)

if (WIN32)
  # for msvc 8 deprecated warning removal
  add_definitions(-D_CRT_SECURE_NO_DEPRECATE)
  add_definitions(-D_CRT_NONSTDC_NO_DEPRECATE)
  add_definitions(-D_CRT_SECURE_NO_DEPRECATE_GLOBALS)
  add_definitions(-DWIN32 -D_WIN32)
  add_definitions(-DNOMINMAX)

  #set(CMAKE_DEBUG_POSTFIX "D")

  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} /NODEFAULTLIB:LIBCMT.lib;libcmtd.lib")
endif(WIN32)

if (WIN32)
  add_definitions(-DGL_GLEXT_PROTOTYPES -bigobj)
endif()
# CUDA
option(CUDA_BUILD_EMULATION "enable emulation mode" OFF)
find_package(CUDA REQUIRED)
# FindCUDA.cmake of CMake 2.8.7 does not look for NPP
if (NOT DEFINED CUDA_npp_LIBRARY AND COMMAND find_cuda_helper_libs AND NOT CUDA_VERSION VERSION_LESS "4.0")
  find_cuda_helper_libs(npp)
endif ()
include_directories(${CUDA_INCLUDE_DIRS})
# hide variable as SDK now (usually) part of CUDA Toolkit
mark_as_advanced(FORCE CUDA_SDK_ROOT_DIR)
# additional directories where to look for SDK which are
# not yet added in FindCUDA.cmake of CMake
list(APPEND CUDA_SDK_SEARCH_PATH "${CUDA_TOOLKIT_ROOT_DIR}/samples") 
# look for required SDK header files as suggested in FindCUDA.cmake
find_path(CUDA_HELPER_INCLUDE_DIR
  helper_cuda.h
  PATHS         ${CUDA_SDK_SEARCH_PATH}
  PATH_SUFFIXES "common/inc"
  DOC           "Location of helper_cuda.h of the CUDA SDK."
  NO_DEFAULT_PATH
)
find_path(CUDA_HELPER_INCLUDE_DIR helper_cuda.h DOC "Location of helper_cuda.h of the CUDA SDK.")
mark_as_advanced(FORCE CUDA_HELPER_INCLUDE_DIR)
# either add found path to include search path or show user the
# appropriate CMake variables to tell us where to find the headers
if (CUDA_HELPER_INCLUDE_DIR)
  include_directories(${CUDA_HELPER_INCLUDE_DIR})
elseif (NOT CUDA_SDK_ROOT_DIR)
  mark_as_advanced(CLEAR CUDA_SDK_ROOT_DIR)
else ()
  mark_as_advanced(CLEAR CUDA_HELPER_INCLUDE_DIR)
endif ()

if( ${CUDA_VERSION} VERSION_LESS "5")
  message(INFO "Building with CUDA < 5.0")
else()
  message(STATUS "Building with CUDA >= 5.0")
  # we build for all useful compute capabilities (C.P > 2.0)
  set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS};-gencode;arch=compute_30,code=sm_30;-use_fast_math;")

endif()


if(BUILD_OUT_OF_IRTK_TREE)

# we do not need to search for IRTK for in source tree build
find_package(IRTK REQUIRED)
message("${IRTK_LIBRARIES_DIR}")
include_directories(${IRTK_INCLUDE_DIRS})

#set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_PATH} ${IRTK_DIR}/../config)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

FIND_PACKAGE(TBB REQUIRED)
IF (TBB_FOUND)
    # Attention: DO NOT define TBB_DEPRECATED by default or before including the
    #            other TBB header files, in particular parallel_for. The deprecated
    #            behavior of parallel_for is to not choose the chunk size (grainsize)
    #            automatically!
    #
    # http://software.intel.com/sites/products/documentation/doclib/tbb_sa/help/tbb_userguide/Automatic_Chunking.htm
    ADD_DEFINITIONS(-DHAS_TBB)
    INCLUDE_DIRECTORIES(${TBB_INCLUDE_DIRS})
    LINK_DIRECTORIES(${TBB_LIBRARY_DIRS})
    LINK_LIBRARIES(${TBB_LIBRARIES})
ENDIF (TBB_FOUND)

FIND_PACKAGE(OpenCV REQUIRED)
include_directories(${OpenCV_DIR}/include ${OpenCV_DIR}/include/opencv)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

else(BUILD_OUT_OF_IRTK_TREE)

include_directories(../../../common++/include ../../../image++/include ../../../geometry++/include ../../../contrib++/include)

SET(IRTK_LIBRARIES contrib++ image++ geometry++ common++) 
IF (USE_VXL)
  SET(IRTK_LIBRARIES irtkvnl irtkvcl irtkvnl_algo irtknetlib ${IRTK_LIBRARIES})
ELSE (USE_VXL)
  SET(IRTK_LIBRARIES recipes ${IRTK_LIBRARIES})
ENDIF (USE_VXL)

endif(BUILD_OUT_OF_IRTK_TREE)
link_directories(${IRTK_LIBRARIES_DIR})

include_directories(../../include)

#todo needs to go in extra lib //
cuda_add_library(reconstruction_cuda_lib STATIC reconstruction_cuda.cu reconstruction_cuda.cuh gaussfilter.cu SKconvolution.cu sum.cu gaussFilterConvolution.cuh)

if(NOT BUILD_OUT_OF_IRTK_TREE)
INSTALL_TARGETS(/lib reconstruction_cuda_lib)
endif(NOT BUILD_OUT_OF_IRTK_TREE)

INCLUDE_DIRECTORIES(${TBB_INSTALL_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR} ${IRTK_DIR}/../external/gco-v3.0 ${IRTK_DIR}/../packages/registration/include ${IRTK_DIR}/../packages/segmentation/include)

cuda_add_executable(reconstructionCuda reconstruction_cuda.cpp irtkReconstructionCuda.cc ../../include/irtkReconstructionCuda.h)
target_link_libraries(reconstructionCuda ${IRTK_LIBRARIES} ${TBB_LIBRARIES} ${OpenCV_LIBS} ${VTK_LIBRARIES} reconstruction_cuda_lib ${CUDA_cufft_LIBRARY})
