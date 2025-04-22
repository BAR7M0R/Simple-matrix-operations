add_library(Matrix_tools STATIC Src/Matrix_tools.cpp)
target_compile_features(Matrix_tools PUBLIC cxx_std_20)
target_include_directories(Matrix_tools PUBLIC Inc)
target_link_libraries(Matrix_tools PUBLIC Matrix)

