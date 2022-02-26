add_library(swatlib swatlib.h
matrix.cpp matrix.h
similarity.cpp similarity.h
encoding.cpp encoding.h
fasta.cpp fasta.h
timing.hpp
vector.hpp
)

target_include_directories(swatlib PUBLIC "${PROJECT_SOURCE_DIR}/extern/include")