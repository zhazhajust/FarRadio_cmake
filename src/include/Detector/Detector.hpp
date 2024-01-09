#ifndef DETECTOR_H
#define DETECTOR_H
#include <vector>
#include <memory>
#include <iostream>
#include <pybind11/pybind11.h>
#include <unsupported/Eigen/CXX11/Tensor>
#include "utils.hpp"
#include "Field/Field.hpp"
#include "Tracer/Tracer.hpp"
#include "Eigen/Core"
#include "Eigen/Dense"
#include "pybind11/eigen.h"

#ifndef NONMPI
    #include "FaradioMPI/FaradioMPI.hpp"
#endif

namespace py = pybind11;
using namespace std;

typedef Eigen::Tensor<double, 3> Tensor3D;
typedef Eigen::Matrix<double, 3, 1> Vec3d;
typedef Eigen::Matrix<double, 3, 1> Vec1d;
typedef Eigen::Matrix<double, Eigen::Dynamic, 3> Vec3dArr;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vec1dArr;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixXd;

class SpheDetector
{
public:
    SpheDetector(vector<double> dmin, vector<double> dmax, 
        vector<int> nf);

    ~SpheDetector();

    void init_screen();
    
    std::shared_ptr<Field3D> get_emf(){
        return this->emf;
    }
    std::shared_ptr<Field3D> get_screen_potisions(){
        return this->screen_potisions;
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
#ifndef NONMPI
    std::shared_ptr<FaradioMPI> get_mpi(){
        return this->faradio_mpi;
    }
#endif
    void check_boundary(int l, int j, int k);
    void deposite_potential(double time_ret, double time_ret_prev, const Vec3d& far_field, int j, int k);
    void cmp_emf_single_particle(const Vec3d& position, const Vec3d& position_prev, 
    const Vec3d& beta, const Vec3d& beta_prev, double time, double charge, double dt);
    void cmp_emf(Eigen::Ref<const Vec3dArr> position_arr, 
        Eigen::Ref<const Vec3dArr> position_prev_arr, Eigen::Ref<const Vec3dArr> beta_arr, 
        Eigen::Ref<const Vec3dArr> beta_prev_arr, double time, Eigen::Ref<const Vec1dArr> charge, double dt);
    void reduce();

    void set_approx(bool if_approx){
        this->if_approx = if_approx;
    }

protected:
    bool if_approx;
    std::shared_ptr<Field3D> emf;
    std::shared_ptr<Field3D> screen_potisions;
#ifndef NONMPI
    std::shared_ptr<FaradioMPI> faradio_mpi;
#endif
    vector<double> time_det;
    vector<double> dmax;
    vector<double> dmin;
    vector<int> nf;
    double d1;
    double d2;
    double d3;
    double R;
};

#endif