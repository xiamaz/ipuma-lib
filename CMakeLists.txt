cmake_minimum_required(VERSION 3.10)
project(ipuma-lib LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" )
include(AddGPLibrary) 
include(CMakeDefinitions.txt)

set(POPLIBS_ENABLED_IPU_ARCH_NAMES "ipu2")
# if(NOT DEFINED POPLIBS_ENABLED_IPU_ARCH_NAMES)
#     set(POPLIBS_ENABLED_IPU_ARCH_NAMES "cpu;ipu1;ipu2" CACHE STRING
#         "The architectures to build the poplibs library objects for")
# endif()
message("Enabled architectures: ${POPLIBS_ENABLED_IPU_ARCH_NAMES}")

# set(DEFAULT_TEST_VARIANTS "cpu;IpuModel2;")
set(DEFAULT_TEST_VARIANTS "cpu;")

set(POPC_EXECUTABLE "popc")

# Set up popc flags.
set(POPC_FLAGS -I ${CMAKE_SOURCE_DIR}/include)

if(CMAKE_BUILD_TYPE STREQUAL Debug)
  list(APPEND POPC_FLAGS -g)
elseif(CMAKE_BUILD_TYPE STREQUAL RelWithDebInfo)
  list(APPEND POPC_FLAGS -O3 -g)
else()
  list(APPEND POPC_FLAGS -O3)
endif()

# If requested, optimise the debug builds too.
if (OPTIMISE_IPU_DEBUG AND CMAKE_BUILD_TYPE STREQUAL Debug)
  list(APPEND POPC_FLAGS -O3)
endif()



include_directories("extern/plog/include")

add_subdirectory("extern/include/ssw")
list (APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/extern/seqan3/build_system")

add_subdirectory(src)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  set(PACKAGE_TESTS ON)
  set(PACKAGE_MAIN ON)
  set(PACKAGE_CPU ON)
endif()

set(IPUMA_BIN_INCLUDES "${PROJECT_SOURCE_DIR}/src" "${PROJECT_SOURCE_DIR}/extern/cxxopts/include" "${PROJECT_SOURCE_DIR}/extern/include")

if(PACKAGE_MAIN)
  add_subdirectory(ipusw)
endif()

if(PACKAGE_CPU)
  add_subdirectory(cpusw)
endif()

if(PACKAGE_TESTS)
  enable_testing()
  add_subdirectory(test)
  # add_subdirectory(test/core)
endif()

#   - POPC_EXECUTABLE
#   - POPC_FLAGS