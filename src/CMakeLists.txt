add_gp_library(
  NAME algoipu
  ASM_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/codelets/codelet_batch_asm.S
    ${CMAKE_CURRENT_SOURCE_DIR}/codelets/codelet_multi_asm.S
  CPP_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/codelets/helper_utils.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/codelets/codelet_batch_affine.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/codelets/codelet_batch_asm.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/codelets/codelet_multi_asm.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/codelets/codelet_multi_affine.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/codelets/codelet_xdrop.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/codelets/codelet_greedy_xdrop.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/codelets/codelet_multi_band_xdrop.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/codelets/codelet_multi_xdrop.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/codelets/codelet_seed_extend_multi_xdrop.cpp
  HEADERS
)

add_subdirectory(swatlib)

add_library(IPUMA_LIBRARY "ipu_batch_affine.cpp" "ipu_base.cpp" "driver.cpp" "partition.cpp" "types.cpp" "ipu_config.cpp")

find_package(OpenMP REQUIRED)
if(OPENMP_FOUND OR OpenMP_CXX_FOUND)
    # Set THREADED if OpenMP is found
    message("Threaded")
    target_compile_definitions(IPUMA_LIBRARY PUBLIC THREADED)
    if(TARGET OpenMP::OpenMP_CXX)
        target_link_libraries(IPUMA_LIBRARY PUBLIC OpenMP::OpenMP_CXX)
    else()
        target_compile_options(IPUMA_LIBRARY PUBLIC "${OpenMP_CXX_FLAGS}")
        target_link_libraries(IPUMA_LIBRARY PUBLIC "${OpenMP_CXX_FLAGS}")
    endif()
endif()
target_include_directories(IPUMA_LIBRARY PUBLIC "${PROJECT_SOURCE_DIR}/extern/include")
target_link_libraries(IPUMA_LIBRARY PUBLIC poplar popops poputil swatlib)
add_dependencies(IPUMA_LIBRARY  algoipu.gp_codelets)

include(FetchContent)
if (NOT channel_POPULATED)
    FetchContent_Declare(channel URL https://github.com/andreiavrammsd/cpp-channel/archive/refs/tags/v0.7.2.tar.gz)
    FetchContent_Populate(channel)
    target_include_directories(IPUMA_LIBRARY PUBLIC ${channel_SOURCE_DIR}/include)
    # OR
    # add_subdirectory(${channel_SOURCE_DIR}/)
endif()