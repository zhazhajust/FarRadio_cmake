#ifndef FIELD_H
#define FIELD_H

#include <memory>
#include <pybind11/pybind11.h>
#include <iostream>

namespace py = pybind11;

class Field
{
public:
    Field();
    Field(int dim);
    virtual ~Field();
    
    inline double& __attribute__((always_inline)) operator()(int i)
    {
        //DEBUGEXEC( if( i>=dims_[0] || j>=dims_[1] || k >= dims_[2] ) ERROR( name << "Out of limits "<< i << " " << j << " " << k ) );
        //return data[i][j][k];
        
        //unsigned int ni = dims_[0];
        //unsigned int nj = dims_[1];
        //unsigned int nk = dims_[2];

        return this->get_data()[i];
    };

//protected:

    //double* get_data(){
    //    return this->_data;
    //};

    double* get_data(){
        return this->_data.get();
    };

    int global_dim;
    //double* _data;
    std::shared_ptr<double> _data;
};

class Field3D : public Field
{
public:
    Field3D();
    Field3D(int size_0, int size_1, int size_2);
    ~Field3D();
    inline double& __attribute__((always_inline)) operator()(int i, int j, int k)
    {
        //DEBUGEXEC( if( i>=dims_[0] || j>=dims_[1] || k >= dims_[2] ) ERROR( name << "Out of limits "<< i << " " << j << " " << k ) );
        //return data[i][j][k];
        
        int ni = dims_[0];
        int nj = dims_[1];
        int nk = dims_[2];
        //std::cout << "idx: " << k + j * nk + i * nk * nj << std::endl;
        return this->get_data()[k + j * nk + i * nk * nj];
    };

    inline void setValue(int i, int j, int k, double value){
        int ni = dims_[0];
        int nj = dims_[1];
        int nk = dims_[2];
        //std::cout << "setValue: " << k + j * nk + i * nk * nj << " to be " << value << std::endl;
        
        this->get_data()[k + j * nk + i * nk * nj] = value;
        return;
    };

    inline int __attribute__((always_inline)) get_dim(int i)
    {
        if (i > 2){
            return 0;
        };
        return this->dims_[i];
    };

    /*
    //to python memoryview
    inline py::memoryview to_memoryview(){
        int dim[3];
        dim[0] = this->get_dim(0);
        dim[1] = this->get_dim(1);
        dim[2] = this->get_dim(2);
        return py::memoryview::from_buffer(
            this->get_data(),               // buffer pointer
            {dim[0], dim[1], dim[2]},                 // shape (rows, cols)
            {sizeof(double) * dim[1] * dim[2], sizeof(double) * dim[2], sizeof(double)} // strides
        );
    };

    //to python memoryview
    inline py::memoryview to_memoryview_1d(){
        return py::memoryview::from_memory(
            this->get_data(),               // buffer pointer
            sizeof(double) * this->get_dim(0) * this->get_dim(1) * this->get_dim(2) // strides
        );
    };
    */
   
    //double ***data;

//protected:
    int dims_[3];
};

#endif
