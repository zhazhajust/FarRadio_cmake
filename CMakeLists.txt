cmake_minimum_required(VERSION 3.4 FATAL_ERROR)
set(CMAKE_CXX_STANDARD 11)
project(faradio)

# Set Compiler
#set(CMAKE_C_COMPILER /opt/homebrew/bin/gcc-13)
#set(CMAKE_CXX_COMPILER /opt/homebrew/bin/g++-13)

# 输出编译信息
message("Start compile!")
message("##############")
aux_source_directory(./src DIRSRCS)
aux_source_directory(./src/include/Detector DIRSRCS)
aux_source_directory(./src/include/Field DIRSRCS)
aux_source_directory(./src/include/Tracer DIRSRCS)

# OPENMP
find_package(OpenMP)
if(OpenMP_FOUND)
    message(STATUS "found openmp")
    include_directories(${OpenMP_INCLUDE_PATH})
    link_libraries(${OpenMP_LIBRARIES})
    if(APPLE)
        add_compile_options(-Xpreprocessor -fopenmp)
        link_libraries(-lomp)
    else()
        add_compile_options(-lgomp)
    endif()
else()
    message("openmp not found!")
endif()

# DEBUG
if(CMAKE_BUILD_TYPE AND (CMAKE_BUILD_TYPE STREQUAL "Debug"))
    add_definitions(-DDEBUG)
else(CMAKE_BUILD_TYPE AND (CMAKE_BUILD_TYPE STREQUAL "Release"))
    add_definitions(-DEIGEN_NO_DEBUG)
endif(CMAKE_BUILD_TYPE)

# add mpi
find_package(MPI)
if(MPI_FOUND)    
    include_directories(${MPI_INCLUDE_PATH})
    link_libraries(${MPI_LIBRARIES})
    aux_source_directory(./src/include/FaradioMPI DIRSRCS)
    message(STATUS "MPI library found")
    message(${MPI_LIBRARIES})
else(MPI_FOUND)
    add_definitions(-DNONMPI)
    message("MPI library not found")
endif(MPI_FOUND)

# list source dir
message(${DIRSRCS})

# 添加include头文件
include_directories(./src/include ./pybind11/include)

# 添加子目录
add_subdirectory(pybind11)
pybind11_add_module(faradio ${DIRSRCS})

# EXAMPLE_VERSION_INFO is defined by setup.py and passed into the C++ code as a
# define (VERSION_INFO) here.
target_compile_definitions(faradio
                           PUBLIC VERSION_INFO=${EXAMPLE_VERSION_INFO})