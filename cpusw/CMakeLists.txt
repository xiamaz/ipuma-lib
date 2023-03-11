add_executable(cpusw main.cpp)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
target_link_libraries(cpusw PRIVATE Threads::Threads)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(cpusw PUBLIC OpenMP::OpenMP_CXX)
endif()

add_library( genometools SHARED IMPORTED )
set_target_properties( genometools PROPERTIES IMPORTED_LOCATION ${PROJECT_SOURCE_DIR}/extern/genometools/lib/libgenometools.a )

target_compile_definitions(cpusw PUBLIC WITHOUT_CAIRO)
target_compile_features(cpusw PRIVATE cxx_std_17)
target_include_directories(cpusw PRIVATE
    ${IPUMA_BIN_INCLUDES}
    "${PROJECT_SOURCE_DIR}/extern/seqan/include"
    "${PROJECT_SOURCE_DIR}/extern/genometools/src"
    "${PROJECT_SOURCE_DIR}/extern/genometools/obj"
)
target_link_libraries (cpusw PRIVATE SSW_LIBRARY IPUMA_LIBRARY genometools)
set_target_properties(cpusw PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")