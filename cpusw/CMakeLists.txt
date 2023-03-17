add_executable(cpusw main.cpp)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
target_link_libraries(cpusw PRIVATE Threads::Threads)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -march=native -mavx2")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -mavx2")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -DSEQAN_ENABLE_DEBUG=0 -g -Wall -Wpedantic -fsanitize=address")

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(cpusw PUBLIC OpenMP::OpenMP_CXX)
endif()

add_library(ksw2 ${PROJECT_SOURCE_DIR}/extern/ksw2/ksw2_extz2_sse.c)
target_include_directories(ksw2 PUBLIC "${PROJECT_SOURCE_DIR}/extern/ksw2")

add_library( genometools SHARED IMPORTED )
set_target_properties(
    genometools 
    PROPERTIES
    IMPORTED_LOCATION ${PROJECT_SOURCE_DIR}/extern/genometools/lib/libgenometools.a
)

add_library(libgaba SHARED IMPORTED )
set_target_properties(
    libgaba 
    PROPERTIES
    IMPORTED_LOCATION ${PROJECT_SOURCE_DIR}/extern/libgaba/libgaba.a
)

target_compile_definitions(cpusw PUBLIC WITHOUT_CAIRO)
target_compile_features(cpusw PRIVATE cxx_std_17)
target_include_directories(cpusw PRIVATE
    ${IPUMA_BIN_INCLUDES}
    "${PROJECT_SOURCE_DIR}/extern/seqan/include"
    "${PROJECT_SOURCE_DIR}/extern/genometools/src"
    "${PROJECT_SOURCE_DIR}/extern/genometools/obj"
    "${PROJECT_SOURCE_DIR}/extern/libgaba"
    "${PROJECT_SOURCE_DIR}/extern/ksw2"
)
target_link_libraries (cpusw PRIVATE SSW_LIBRARY IPUMA_LIBRARY genometools libgaba ksw2)
set_target_properties(cpusw PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
