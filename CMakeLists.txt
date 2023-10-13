cmake_minimum_required(VERSION 3.4 FATAL_ERROR)
set(CMAKE_CXX_STANDARD 11)
project(faradio)

# 输出编译信息
message("Start compile!")
message("##############")
aux_source_directory(./src DIRSRCS)
aux_source_directory(./src/include DIRSRCS)
message(${DIRSRCS})

# 添加include头文件
include_directories(./src/include ./src)

add_subdirectory(pybind11)
#pybind11_add_module(faradio src/wraper.cpp src/include/Detector.cpp src/include/Field.cpp src/include/Tracer.cpp)
pybind11_add_module(faradio ${DIRSRCS})

# EXAMPLE_VERSION_INFO is defined by setup.py and passed into the C++ code as a
# define (VERSION_INFO) here.
target_compile_definitions(faradio
                           PUBLIC VERSION_INFO=${EXAMPLE_VERSION_INFO})
