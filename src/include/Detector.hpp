#ifndef DETECTOR_H
#define DETECTOR_H
#include <iostream>
#include <vector>
#include "Field.hpp"
#include "Vector.hpp"
#include "Tracer.hpp"
#include "pybind11/eigen.h"
#include "Eigen/Dense"
#include "Eigen/Core"
#include <pybind11/pybind11.h>
#include <unsupported/Eigen/CXX11/Tensor>
#include <memory>

namespace py = pybind11;
using namespace std;

typedef Eigen::Matrix<double, 3, 1> Vec3d;
//typedef Eigen::MatrixX3d Vec3dArr;
typedef Eigen::Matrix<double, Eigen::Dynamic, 3> Vec3dArr;
typedef Eigen::Tensor<double, 3> Tensor3D;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixXd;

class SpheDetector
{
public:
    SpheDetector(vector<double> dmin, vector<double> dmax, 
        vector<int> nf);

    ~SpheDetector();

    void init_screen();

    /*
    Field3D& get_emf(){
        return *(this->emf);
    }
    Field3D& get_screen_potisions(){
        return *(this->screen_potisions);
    }
    */
    
    std::shared_ptr<Field3D> get_emf(){
        return this->emf;
    }
    std::shared_ptr<Field3D> get_screen_potisions(){
        return this->screen_potisions;
    }
    

    /*
    inline py::memoryview get_emf(){
        return this->emf->to_memoryview();
    }

    inline py::memoryview get_screen_potisions(){
        return this->screen_potisions->to_memoryview();
    }
    */

    /*
    const Tensor3D get_emf(){
        return this->emf;
    }

    const Tensor3D get_screen_potisions(){
        return this->screen_potisions;
    }
    */
   
    const MatrixXd get_screen_x(){
        return screen_x;
    }
    
    vector<double> get_dmax(){
        return this->dmax;
    }
    vector<double> get_dmin(){
        return this->dmin;
    }
    vector<int> get_nf(){
        return this->nf;
    }

    //void cmp_emf(Particle& tracer);

    void cmp_emf(Eigen::Ref<const Vec3dArr> position_arr, 
        Eigen::Ref<const Vec3dArr> position_prev_arr, Eigen::Ref<const Vec3dArr> beta_arr, 
        Eigen::Ref<const Vec3dArr> beta_prev_arr, double time, double charge, int particle_nums, double dt);

protected:
    //Field3D* emf;
    //Field3D* screen_potisions;

    std::shared_ptr<Field3D> emf;
    std::shared_ptr<Field3D> screen_potisions;

    //Tensor3D emf;
    //Tensor3D screen_potisions;

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> screen_x;

    vector<double> dmax;
    vector<double> dmin;
    vector<int> nf;
    double d1;
    double d2;
    double d3;
};

#endif