cmake_minimum_required(VERSION 3.10)
project(ipuma-lib LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" )
include(AddGPLibrary) 
include(CMakeDefinitions.txt)

set(POPLIBS_ENABLED_IPU_ARCH_NAMES "cpu;ipu2")
# if(NOT DEFINED POPLIBS_ENABLED_IPU_ARCH_NAMES)
#     set(POPLIBS_ENABLED_IPU_ARCH_NAMES "cpu;ipu1;ipu2" CACHE STRING
#         "The architectures to build the poplibs library objects for")
# endif()
message("Enabled architectures: ${POPLIBS_ENABLED_IPU_ARCH_NAMES}")

set(DEFAULT_TEST_VARIANTS "cpu;IpuModel2;")

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


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

enable_testing()

include_directories("plog/include")
add_subdirectory(googletest)

add_subdirectory(src)
add_subdirectory(test)

#   - POPC_EXECUTABLE
#   - POPC_FLAGS