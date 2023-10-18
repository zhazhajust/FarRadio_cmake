#ifndef FIELD_H
#define FIELD_H

#include <memory>
#include <iostream>
#include <pybind11/pybind11.h>
#include "utils.hpp"

namespace py = pybind11;

class Field
{
public:
    Field();
    Field(int dim);
    virtual ~Field();
    
    inline double& __attribute__((always_inline)) operator()(int i)
    {
        return this->get_data()[i];
    };

    inline double& __attribute__((always_inline)) data(int i){
        return this->get_data()[i];
    }

    inline int __attribute__((always_inline)) size(){
        return this->global_dim;
    }

    double* get_data(){
        return this->_data.get();
    };
protected:
    int global_dim;
    std::shared_ptr<double> _data;
};

class Field3D : public Field
{
public:
    Field3D();
    Field3D(int size_0, int size_1, int size_2);
    ~Field3D();

    inline double& __attribute__((always_inline)) data(int i, int j, int k){
        return this->get_data()[k + j * this->dims_[2] + i * this->dims_[1] * this->dims_[2]];
    }

    inline double& __attribute__((always_inline)) operator()(int i, int j, int k)
    {
        //DEBUGEXEC( if( i>=dims_[0] || j>=dims_[1] || k >= dims_[2] ) ERROR( name << "Out of limits "<< i << " " << j << " " << k ) );
        int ni = dims_[0];
        int nj = dims_[1];
        int nk = dims_[2];
        return this->get_data()[k + j * nk + i * nk * nj];
    };

    inline void setValue(int i, int j, int k, double value){
        int ni = dims_[0];
        int nj = dims_[1];
        int nk = dims_[2];
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

protected:
    int dims_[3];
};

#endif
