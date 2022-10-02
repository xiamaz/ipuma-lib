add_executable(xdrop main.cpp)

# find_package (Threads)

target_include_directories(xdrop PRIVATE "${PROJECT_SOURCE_DIR}/src" "${PROJECT_SOURCE_DIR}/extern/cxxopts/include" "${PROJECT_SOURCE_DIR}/extern/include")
# target_link_libraries(xdrop PUBLIC IPUMA_LIBRARY) # ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(xdrop PUBLIC swatlib) # ${CMAKE_THREAD_LIBS_INIT})

set_target_properties(xdrop PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")


list (APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../seqan3/build_system")
add_executable(seqdrop seqxdrop.cpp)

find_package (seqan3 3.0 REQUIRED)
target_link_libraries (seqdrop seqan3::seqan3)