cmake_minimum_required(VERSION 3.5)

project(Simple_matrix_operations
        VERSION 0.5.0
        LANGUAGES CXX)

enable_testing()

set(CMAKE_PREFIX_PATH "C:/Users/Dell/.vcpkg-clion/vcpkg/installed/x64-mingw-dynamic")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

find_package(Matplot++ CONFIG REQUIRED)
find_package(Catch2 3 REQUIRED)

add_subdirectory(Libs/Matrix)
add_subdirectory(Libs/Matrix_tools)
add_executable(my_app main.cpp)
target_compile_features(my_app PRIVATE cxx_std_20)
target_link_libraries(my_app PRIVATE Matrix)

