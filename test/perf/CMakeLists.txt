set(BINARY ipuma-perf)

set(TEST_SOURCES main.cpp ipuma-perf.cpp)

set(SOURCES ${TEST_SOURCES})

add_executable(${BINARY} ${TEST_SOURCES} ${MHM2_TARGET_OBJECTS})

target_include_directories(ipuma-perf PRIVATE "${CMAKE_SOURCE_DIR}/src" "${CMAKE_SOURCE_DIR}/extern/include")

# add_test(NAME ${BINARY} COMMAND ${BINARY})

set_target_properties( ${BINARY}
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/test"
)

target_link_libraries(
  ${BINARY} PUBLIC gtest IPUMA_LIBRARY)
install(TARGETS ${BINARY} DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)