#include <vector>
#include <iostream>
#include "Vector.hpp"
//#include "Field.hpp"
#include "Detector.hpp"
//#include "Tracer.hpp"
//#include "input_data.h"

using namespace std;

#include <math.h>

#define PI acos(-1)

#define Q 1.60217662e-19 // electron charge
#define C 299792458 // speed of light

#define calcu_far_field(n, beta_prev, beta_dot, R) \
    n.cross((n - beta_prev).cross(beta_dot)) / pow(1 - beta_prev.dot(n), 3) / R

#define ALL Eigen::placeholders::all

/*
R_vec cmp_field(R_vec n, R_vec beta_prev, R_vec beta_dot, double R)
{
    return n % ((n - beta_prev) % beta_dot) / pow(1 - beta_prev * n, 3) / R;
};
*/

SpheDetector::SpheDetector(vector<double> dmin, vector<double> dmax, 
        vector<int> nf): dmin(dmin), dmax(dmax), nf(nf)
{
    // delta of time, theta, phi
    d1 = (dmax[0] - dmin[0])/nf[0];
    d2 = (dmax[1] - dmin[1])/nf[1];
    d3 = (dmax[2] - dmin[2])/nf[2];

    // grid of time, theta, phi
    //emf = new Field3D(nf[0], nf[1], nf[2]);
    //std::shared_ptr<Field3D> emf (new Field3D(nf[0], nf[1], nf[2]));
    this->emf = std::make_shared<Field3D>(nf[0], nf[1], nf[2]);

    //emf = Tensor3D(nf[0], nf[1], nf[2]);
    //emf.setZero();

    // grid of screen positions: theta, phi, {x, y, z}
    //screen_potisions = new Field3D(nf[1], nf[2], 3);
    this->screen_potisions = std::make_shared<Field3D>(nf[1], nf[2], 3);

    //screen_potisions = Tensor3D(nf[1], nf[2], 3);
    //screen_potisions.setZero();

    //screen_x = MatrixXd(nf[1], nf[2]);

    screen_x = Eigen::MatrixXd::Zero(nf[1], nf[2]);

    // write buffer
    //write_buffer = new double [nf[0]*nf[1]*nf[2]];

    // init screen positions
    this->init_screen();

    cout << "Address:" << hex << long(this) << " SpheDetector constructor called" << endl;
    //vector<vector<vector<double>>> e4(this->nf[2],this->nf[3],3);
};

SpheDetector::~SpheDetector()
{
    //delete this->write_buffer;

    //delete this->emf;
    //delete this->screen_potisions;

    cout << "Address:" << hex << long(this) << " SpheDetector destructor called" << endl;
};


//void SpheDetector::cmp_emf(vector<R_vec> position, vector<R_vec> position_prev,
//        vector<R_vec> beta, vector<R_vec> beta_prev, double time, 
//        vector<double> charge, int particle_nums, double dt)

void SpheDetector::init_screen(){
    for(int j = 0; j < this->nf[1]; j++){
        for(int k = 0; k < this->nf[2]; k++){
 
            double theta = this->dmin[1]+j*this->d2;
            double phi = this->dmin[2]+k*this->d3;

            (*this->screen_potisions)(j, k, 0) = sin(theta)*cos(phi);
            (*this->screen_potisions)(j, k, 1) = sin(theta)*sin(phi);
            (*this->screen_potisions)(j, k, 2) = cos(theta);

            //this->screen_potisions->setValue(j, k, 0, sin(theta)*cos(phi));
            //this->screen_potisions->setValue(j, k, 1, sin(theta)*sin(phi));
            //this->screen_potisions->setValue(j, k, 2, cos(theta));

            //this->screen_potisions(j, k, 0) = sin(theta)*cos(phi);
            //this->screen_potisions(j, k, 1) = sin(theta)*sin(phi);
            //this->screen_potisions(j, k, 2) = cos(theta);

            this->screen_x(j, k) = sin(theta)*cos(phi);

            //auto temp = this->screen_potisions->get_data();
            //temp[0 + k * this->nf[2] + j * this->nf[2] * this->nf[1]] = sin(theta)*cos(phi);
            //temp[1 + k * this->nf[2] + j * this->nf[2] * this->nf[1]] = sin(theta)*sin(phi);
            //temp[2 + k * this->nf[2] + j * this->nf[2] * this->nf[1]] = cos(theta);

        }
    };

    /*
    for(int j = 0; j < this->nf[1]; j++){
        for(int k = 0; k < this->nf[2]; k++){
            cout << (*this->screen_potisions)(j, k, 0) << " ";
        }
    };
    */
}

/*
void SpheDetector::cmp_emf_single_particle(Eigen::Ref<const Vec3d> position, 
    Eigen::Ref<const Vec3d> position_prev, Eigen::Ref<const Vec3d> beta, 
    Eigen::Ref<const Vec3d> beta_prev, double time, double charge, int particle_nums, double dt){

    double time_prev = time - dt;
    const Vec3d beta_dot = (beta - beta_prev) / dt;

    for(int j = 0; j < this->nf[1]; j++){
        for(int k = 0; j < this->nf[2]; j++){
            Vec3d n((*this->screen_potisions)(j, k, 0), 
                    (*this->screen_potisions)(j, k, 1), 
                    (*this->screen_potisions)(j, k, 2));

            Vec3d far_field = Q / C * calcu_far_field(n, beta_prev[i], beta_dot[i], R);

            // interpolate linearly
            double time_ret = time - n*position[i];
            double time_ret_prev = time_prev - n*position_prev[i];

            unsigned int idx_time = (time_ret - this->dmin[0]) / this->d1 + 1;
            unsigned int idx_time_prev = (time_ret_prev - this->dmin[0]) / this->d1 + 1;

            for(int l = idx_time_prev; l < idx_time; l++){
                double temp = (l * this->d1 - time_ret_prev) / (time_ret - time_ret_prev);
                (*this->emf)(idx_time, j, k) += temp * far_field[1];
            }
            //
        }
    }

}
*/

/*
void SpheDetector::interp_1d(double time_ret, double time_ret_prev, Vec3d far_field, double& data_ref){
    unsigned int idx_time = (time_ret - this->dmin[0]) / this->d1 + 1;
    unsigned int idx_time_prev = (time_ret_prev - this->dmin[0]) / this->d1 + 1;

    for(int l = idx_time_prev; l < idx_time; l++){
        double temp = (l * this->d1 - time_ret_prev) / (time_ret - time_ret_prev);
        data_ref += temp * far_field[1];
    }
}
*/

#define INTERP1D(time_ret, time_ret_prev, far_field) \
    int idx_time = (time_ret - this->dmin[0]) / this->d1 + 1; \
    int idx_time_prev = (time_ret_prev - this->dmin[0]) / this->d1 + 1; \
    for(int l = idx_time_prev; l < idx_time; l++){ \
        double temp = (l * this->d1 - time_ret_prev) / (time_ret - time_ret_prev); \
        this->emf->operator()(idx_time, j, k) += temp * far_field[1]; \
    }
        //(*this->emf)(idx_time, j, k) += temp * far_field[1]; \

        
void SpheDetector::cmp_emf(Eigen::Ref<const Vec3dArr> position_arr, 
    Eigen::Ref<const Vec3dArr> position_prev_arr, Eigen::Ref<const Vec3dArr> beta_arr, 
    Eigen::Ref<const Vec3dArr> beta_prev_arr, double time, double charge, int particle_nums, double dt)
{
    double time_prev = time - dt;
    double R = this->dmin[4];
    Vec3dArr beta_dot_arr(position_arr.rows(), 3);
    beta_dot_arr = (beta_arr - beta_prev_arr) / dt;

    for(int i = 0; i < particle_nums; i++){
        //theta
        for(int j = 0; j < this->nf[1]; j++){
            //phi
            for(int k = 0; j < this->nf[2]; j++){
                //Vec3d n((*this->screen_potisions)(j, k, 0), 
                //        (*this->screen_potisions)(j, k, 1), 
                //        (*this->screen_potisions)(j, k, 2));

                Vec3d n(this->screen_potisions->operator()(j, k, 0), this->screen_potisions->operator()(j, k, 1), 
                    this->screen_potisions->operator()(j, k, 2));

                Vec3d far_field = Q / C * calcu_far_field(n, Vec3d (beta_prev_arr(i, ALL)), 
                    Vec3d (beta_dot_arr(i, ALL)), R);

                // interpolate linearly
                double time_ret = time - n.dot(position_arr(i, ALL));
                double time_ret_prev = time_prev - n.dot(position_prev_arr(i, ALL));
                INTERP1D(time_ret, time_ret_prev, far_field);
            }
        }
    }
}

/*
void SpheDetector::cmp_emf(Particle& tracer)
{
    double time = tracer.get_time();
    double dt = tracer.get_dt();
    vector<R_vec> position = tracer.get_position();
    vector<R_vec> position_prev = tracer.get_position_prev();
    vector<R_vec> beta = tracer.get_beta();
    vector<R_vec> beta_prev = tracer.get_beta_prev();
    double charge = tracer.get_charge();
    int particle_nums = position.size();

    vector<int> iters(3, 0);

    double R = this->dmin[4];
    double time_prev = time - dt;
    vector<R_vec> beta_dot(particle_nums);
    
//#ifdef OMP
//#pragma omp parallel for
//#endif
    for(int i = 0; i < particle_nums; i++){
        beta_dot[i] = (beta[i] - beta_prev[i]) / dt;
    }

//#ifdef OMP
//#pragma omp parallel for
//#endif
    for(int i = 0; i < particle_nums; i++){
        for(int j = 0; j < this->nf[1]; j++){
            for(int k = 0; j < this->nf[2]; j++){
                double theta = this->dmin[1]+k*this->d2;
                double phi = this->dmin[2]+j*this->d3;
                (*this->screen_potisions)(k+1, j+1, 0) = sin(theta)*cos(phi);
                (*this->screen_potisions)(k+1, j+1, 1) = sin(theta)*sin(phi);
                (*this->screen_potisions)(k+1, j+1, 2) = cos(theta);
            }
        }
    };

//#ifdef OMP
//#pragma omp parallel for
//#endif
    for(int i = 0; i < particle_nums; i++){
        for(int j = 0; j < this->nf[1]; j++){
            for(int k = 0; j < this->nf[2]; j++){
                R_vec n((*this->screen_potisions)(j, k, 0), 
                        (*this->screen_potisions)(j, k, 1), 
                        (*this->screen_potisions)(j, k, 2));

                R_vec far_field = Q / C * far_field(n, beta_prev[i], beta_dot[i], R);

                // interpolate linearly
                double time_ret = time - n*position[i];
                double time_ret_prev = time_prev - n*position_prev[i];

                unsigned int idx_time = (time_ret - this->dmin[0]) / this->d1 + 1;
                unsigned int idx_time_prev = (time_ret_prev - this->dmin[0]) / this->d1 + 1;
                for(int l = idx_time_prev; l < idx_time; l++){
                    double temp = (l * this->d1 - time_ret_prev) / (time_ret - time_ret_prev);
                    //this->emf->data[j][k][0] += far_field[0];
                    (*this->emf)(idx_time, j, k) += temp * far_field[1];
                    //this->emf->data[j][k][2] += far_field[2];
                }
            }
        }
    };
};
*/