add_executable(cpusw main.cpp)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
target_link_libraries(cpusw PRIVATE Threads::Threads)

target_compile_features(cpusw PRIVATE cxx_std_17)
target_include_directories(cpusw PRIVATE ${IPUMA_BIN_INCLUDES} ${SEQAN3_INCLUDE_DIRS})
target_link_libraries (cpusw PRIVATE SSW_LIBRARY IPUMA_LIBRARY)
set_target_properties(cpusw PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")