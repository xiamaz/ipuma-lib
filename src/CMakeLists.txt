add_gp_library(
  NAME algoipu
  ASM_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/codelets/codelet_batch_asm.S
    ${CMAKE_CURRENT_SOURCE_DIR}/codelets/codelet_multi_asm.S
    ${CMAKE_CURRENT_SOURCE_DIR}/codelets/codelet_striped_asm.S
    ${CMAKE_CURRENT_SOURCE_DIR}/codelets/codelet_multi_striped_asm.S
  CPP_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/codelets/helper_utils.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/codelets/codelet_batch_affine.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/codelets/codelet_batch_asm.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/codelets/codelet_multi_asm.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/codelets/codelet_multi_affine.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/codelets/codelet_multi_striped.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/codelets/codelet_multi_striped_asm.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/codelets/codelet_striped_asm.cpp
  HEADERS
)

add_subdirectory(swatlib)

add_library(IPUMA_LIBRARY "ipu_batch_affine.cpp" "ipu_base.cpp" "driver.cpp" "partition.cpp" "types.cpp" "ipu_config.cpp")
target_link_libraries(IPUMA_LIBRARY poplar popops poputil swatlib)
add_dependencies(IPUMA_LIBRARY  algoipu.gp_codelets)