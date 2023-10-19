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
if(MPI_FOUND)    
    include_directories(${MPI_INCLUDE_PATH})
    link_libraries(${MPI_LIBRARIES})
    message(STATUS "MPI library found")
    message(${MPI_LIBRARIES})
else(MPI_FOUND)
    add_definitions(-DNONMPI)
    message(FATAL_ERROR MPI library not found”)
endif(MPI_FOUND)

# 添加include头文件
include_directories(./src/include ./pybind11/include)

# 添加子目录
add_subdirectory(pybind11)
pybind11_add_module(faradio ${DIRSRCS})

# EXAMPLE_VERSION_INFO is defined by setup.py and passed into the C++ code as a
# define (VERSION_INFO) here.
target_compile_definitions(faradio
                           PUBLIC VERSION_INFO=${EXAMPLE_VERSION_INFO})