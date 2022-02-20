add_subdirectory("${PROJECT_SOURCE_DIR}/googletest" "googletest")

if(ENABLE_POPLAR)
  add_definitions(-DENABLE_IPUS)
  message(STATUS "Building tests for IPUs with poplar")
endif()

set(BINARY ipuma-test)

set(TEST_SOURCES main.cpp simple-test.cpp parity-test.cpp partition-test.cpp dev-test.cpp mn_cmp-test.cpp)

# set(SOURCES ${TEST_SOURCES})

include_directories("${CMAKE_SOURCE_DIR}/src")
include_directories("ssw")

add_subdirectory(ssw)

add_executable(${BINARY} ${TEST_SOURCES})

add_test(NAME ${BINARY} COMMAND ${BINARY})

target_link_libraries(
  ${BINARY} PUBLIC gtest SSW_LIBRARY IPUMA_LIBRARY)
install(TARGETS ${BINARY} DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

# enable "make check"
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} -V)
add_dependencies(check ${BINARY})

add_subdirectory(perf)