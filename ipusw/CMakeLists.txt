add_executable(ipusw main.cpp)

target_include_directories(ipusw PRIVATE "${PROJECT_SOURCE_DIR}/src" "${PROJECT_SOURCE_DIR}/extern/cxxopts/include" "${PROJECT_SOURCE_DIR}/extern/include")
target_link_libraries(ipusw IPUMA_LIBRARY)

set_target_properties(ipusw PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")