add_library(Matrix STATIC Src/Matrix.cpp)
target_compile_features(Matrix PUBLIC cxx_std_20)
target_include_directories(Matrix PUBLIC Inc)

