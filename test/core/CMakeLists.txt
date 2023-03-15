add_subdirectory("${PROJECT_SOURCE_DIR}/extern/googletest" "googletest")
set(BINARY ipuma-core)

set(TEST_SOURCES main.cpp elba_integration.cpp)

set(SOURCES ${TEST_SOURCES})

add_executable(${BINARY} ${TEST_SOURCES} ${MHM2_TARGET_OBJECTS})

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DSEQAN_ENABLE_DEBUG=0")
# set(CMAKE_CXX_FLAGS_DEBUG "-O1 -DSEQAN_ENABLE_DEBUG=0 -g -Wall -Wpedantic -fsanitize=address")


target_include_directories(ipuma-core PRIVATE
    "${CMAKE_SOURCE_DIR}/src"
    "${CMAKE_SOURCE_DIR}/extern/include"
    "${PROJECT_SOURCE_DIR}/extern/seqan/include"
)

add_test(NAME ${BINARY} COMMAND ${BINARY})

set_target_properties( ${BINARY}
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/test"
)

target_link_libraries(
  ${BINARY} PUBLIC gtest IPUMA_LIBRARY SSW_LIBRARY)
install(TARGETS ${BINARY} DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} -V)
add_dependencies(check ${BINARY})