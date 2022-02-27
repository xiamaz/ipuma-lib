# all c++ must be at least c++14

add_library(SSW_LIBRARY ssw.cpp ssw_core.cpp)
target_compile_features(SSW_LIBRARY PRIVATE cxx_std_14)
set_target_properties(SSW_LIBRARY PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS OFF CXX_EXTENSIONS OFF)
