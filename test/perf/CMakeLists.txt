set(BINARY ipuma-perf)

file(
  GLOB TEST_SOURCES
  LIST_DIRECTORIES false
  *.hpp *.cpp)

set(SOURCES ${TEST_SOURCES})

include_directories("${CMAKE_SOURCE_DIR}/src")
include_directories("${CMAKE_SOURCE_DIR}/src/ssw")

add_executable(${BINARY} ${TEST_SOURCES} ${MHM2_TARGET_OBJECTS})

# add_test(NAME ${BINARY} COMMAND ${BINARY})

message("MHM2_LINK_LIBRARIES: ${MHM2_LINK_LIBRARIES}")

set_target_properties( ${BINARY}
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/test"
)

target_link_libraries(
  ${BINARY} PUBLIC ${MHM2_LINK_LIBRARIES} ${UPCXX_UTILS_LIBRARIES}
                   ${UPCXX_LIBRARIES} gtest)
install(TARGETS ${BINARY} DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)