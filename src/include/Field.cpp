#ifndef FIELD_CPP
#define FIELD_CPP
#include <iostream>
#include "Field.hpp"

using namespace std;

Field::Field(int global_dim): global_dim(global_dim)
{
    //this->_data = new double[global_dim];

    //std::shared_ptr<double> p(new double[global_dim], [](double* p) {delete[] p;});
    //std::shared_ptr<double> p(new double[global_dim], std::default_delete<double[]>());
    //this->_data = p;

    this->_data = std::shared_ptr<double> (new double[global_dim], std::default_delete<double[]>());

    //this->_data = std::make_shared<double>(global_dim);

    for(int i = 0; i < global_dim; i++){
        this->get_data()[i] = 0.0;
    };
    cout << "Address:" << hex << long(this) << " Field constructor called" << endl;
};

Field::~Field()
{
    //delete[] this->_data;
    cout << "Address:" << hex << long(this) << " Field destructor called" << endl;
};

Field3D::Field3D(int size_0, int size_1, int size_2): 
    Field(size_0*size_1*size_2), dims_{size_0, size_1, size_2}
{
    //dims_[0] = size_0;
    //dims_[1] = size_1;
    //dims_[2] = size_2;
    //global_dim = size_0 * size_1 * size_2;
    /*
    //initialize 3D array
    data = new double **[size_0];
    for(int i = 0; i < size_1; i++){
        data[i] = new double *[size_1];
        for(int j = 0; j < size_2; j++){
            data[i][j] = new double[size_2];
        }
    };
    */
    cout << "Address:" << hex << long(this) << " Field3D constructor called" << endl;
    //vector<vector<vector<double>>> e4(this->nf[2],this->nf[3],3);
};

Field3D::~Field3D()
{
    /*
    delete this->data;
    */
    cout << "Address:" << hex << long(this) << " Field3D destructor called" << endl;
};

#endif