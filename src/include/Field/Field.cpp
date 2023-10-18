#ifndef FIELD_CPP
#define FIELD_CPP
#include <iostream>
#include "Field.hpp"

using namespace std;

Field::Field(int global_dim): global_dim(global_dim)
{
    this->_data = std::shared_ptr<double> (new double[global_dim], std::default_delete<double[]>());
    for(int i = 0; i < global_dim; i++){
        this->data(i) = 0.0;
    };
    Address(" Field constructor called");
};

Field::~Field()
{
    Address(" Field destructor called");
};

Field3D::Field3D(int size_0, int size_1, int size_2): 
    Field(size_0*size_1*size_2), dims_{size_0, size_1, size_2}
{
    Address(" Field3D constructor called");
};

Field3D::~Field3D()
{
    Address(" Field3D destructor called");
};

#endif