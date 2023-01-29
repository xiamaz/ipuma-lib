add_executable(xdrop main.cpp)

find_package (Threads)

target_include_directories(xdrop PRIVATE "${PROJECT_SOURCE_DIR}/src" "${PROJECT_SOURCE_DIR}/extern/cxxopts/include" "${PROJECT_SOURCE_DIR}/extern/include")
target_link_libraries(xdrop PUBLIC swatlib IPUMA_LIBRARY)

set_target_properties(xdrop PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")