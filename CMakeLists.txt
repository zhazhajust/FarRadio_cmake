cmake_minimum_required(VERSION 3.4 FATAL_ERROR)
set(CMAKE_CXX_STANDARD 11)
project(faradio)

# 输出编译信息
message("Start compile!")
message("##############")
aux_source_directory(./src DIRSRCS)
aux_source_directory(./src/include/Detector DIRSRCS)
aux_source_directory(./src/include/Field DIRSRCS)
aux_source_directory(./src/include/Tracer DIRSRCS)
aux_source_directory(./src/include/FaradioMPI DIRSRCS)
message(${DIRSRCS})

# OPENMP
add_definitions(-DOMP)

# DEBUG
if(CMAKE_BUILD_TYPE AND (CMAKE_BUILD_TYPE STREQUAL "Debug"))
    add_definitions(-DDEBUG)
else(CMAKE_BUILD_TYPE AND (CMAKE_BUILD_TYPE STREQUAL "Release"))
    add_definitions(-DEIGEN_NO_DEBUG)
endif(CMAKE_BUILD_TYPE)

# add mpi
find_package(MPI REQUIRED)

# 添加编译选项
#set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
#set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})

message(${MPI_LIBRARIES})

# 添加include头文件
include_directories(./src/include ./pybind11/include ${MPI_INCLUDE_PATH})

# 添加子目录
add_subdirectory(pybind11)
#pybind11_add_module(faradio src/wraper.cpp src/include/Detector.cpp src/include/Field.cpp src/include/Tracer.cpp)
pybind11_add_module(faradio ${DIRSRCS})

# Link the MPI library
target_link_libraries(faradio PRIVATE ${MPI_LIBRARIES})

# EXAMPLE_VERSION_INFO is defined by setup.py and passed into the C++ code as a
# define (VERSION_INFO) here.
target_compile_definitions(faradio
                           PUBLIC VERSION_INFO=${EXAMPLE_VERSION_INFO})