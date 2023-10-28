#include <math.h>
#include <vector>
#include <iostream>
#include "Detector.hpp"

#ifdef _OPENMP
    #include <omp.h>
#endif

using namespace std;

#define PI acos(-1)
#define Q 1.60217662e-19 // electron charge
#define C 299792458 // speed of light
#define ALL Eigen::placeholders::all
#define NORMAL(x) (sqrt(x(0)*x(0) + x(1)*x(1) + x(2)*x(2)))
#define Lienard_Wiechert(n, beta_prev, beta_dot, R) \
    (n.cross((n - beta_prev).cross(beta_dot)) / pow(1 - beta_prev.dot(n), 3) / R)

Vec3d calcu_field(const Vec3d& n, const Vec3d& beta_prev, const Vec3d& beta_dot, double R){
    return n.cross((n - beta_prev).cross(beta_dot)) / pow(1 - beta_prev.dot(n), 3) / R;
}

SpheDetector::SpheDetector(vector<double> dmin, vector<double> dmax, 
    vector<int> nf): dmin(dmin), dmax(dmax), nf(nf), if_approx(true)
{
    // delta of time, theta, phi
    d1 = (dmax[0] - dmin[0])/nf[0];
    d2 = (dmax[1] - dmin[1])/nf[1];
    d3 = (dmax[2] - dmin[2])/nf[2];
    this->R = dmin[3];
    this->emf = std::make_shared<Field3D>(nf[0], nf[1], nf[2]);
    this->screen_potisions = std::make_shared<Field3D>(nf[1], nf[2], 3);

    // init screen positions
    this->init_screen();
    this->time_det = vector<double>(this->nf[0]);
    for (int k = 0; k < this->nf[0]; k++) {
        double tpos = this->dmin[0] + k * this->d1;
        //double tpos = this->dmin[3]+ this->dmin[0] + k * this->d1;
        this->time_det[k] = tpos;
    };

#ifdef _OPENMP
    cout << "OpenMP is used" << endl;
#pragma omp parallel
    {
        int threadid = omp_get_thread_num();
        if(threadid == 0){
            int nthreads = omp_get_num_threads();
            cout << "nthreads: " << nthreads << endl;
        }
    }
#endif

#ifndef NONMPI
    this->faradio_mpi = FaradioMPI::getFaradioMPI();
#endif
    Address(" SpheDetector constructor called");
};

SpheDetector::~SpheDetector()
{
    Address(" SpheDetector destructor called");
};

void SpheDetector::init_screen(){
    for(int j = 0; j < this->nf[1]; j++){
        for(int k = 0; k < this->nf[2]; k++){
 
            double theta = this->dmin[1]+j*this->d2;
            double phi = this->dmin[2]+k*this->d3;

            this->screen_potisions->data(j, k, 0) = sin(theta)*cos(phi);
            this->screen_potisions->data(j, k, 1) = sin(theta)*sin(phi);
            this->screen_potisions->data(j, k, 2) = cos(theta);

        }
    };
}

void SpheDetector::check_boundary(int l, int j, int k){
    if(l < 0 || l > this->emf->get_dim(0) || j < 0 || j > this->emf->get_dim(1) 
        || k < 0 || k > this->emf->get_dim(2)) return;
}

void SpheDetector::deposite_potential(double time_ret, double time_ret_prev, const Vec3d& far_field, int j, int k){
    int idx_time = (time_ret - this->dmin[0]) / this->d1;
    int idx_time_prev = (time_ret_prev - this->dmin[0]) / this->d1;
    double temp_partial, next_temp;
    for(int l = idx_time_prev; l < idx_time; l++){
        next_temp = (l + 1) * this->d1 + this->dmin[0];
        temp_partial = next_temp - time_ret_prev;
#ifdef _OPENMP
    #pragma omp critical
#endif
        this->emf->data(l, j, k) += temp_partial * far_field[1];
        time_ret_prev = next_temp;
    }
    temp_partial = time_ret - time_ret_prev;
#ifdef _OPENMP
    #pragma omp critical
#endif
    this->emf->data(idx_time, j, k) += temp_partial * far_field[1];
}

#define CHECK_BOUNDARY(l, j, k) \
    if(l < 0 || l > this->emf->get_dim(0) || j < 0 || j > this->emf->get_dim(1) \
        || k < 0 || k > this->emf->get_dim(2)) continue;

void SpheDetector::cmp_emf_single_particle(const Vec3d& position, const Vec3d& position_prev, 
    const Vec3d& beta, const Vec3d& beta_prev, double time, double charge, double dt){
    
    double time_prev = time - dt;
    const Vec3d beta_dot = (beta - beta_prev) / dt;

    //this->R = this->dmin[3]+(this->nf[0]-1)*this->d3;

    //theta
    for(int j = 0; j < this->nf[1]; j++){
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(j, position, position_prev, \
    beta, beta_prev, beta_dot, time, time_prev, charge, dt) schedule(dynamic, 1)
#endif
        //phi
        for(int k = 0; k < this->nf[2]; k++){
//#ifdef APROX
            Vec3d n;
            double time_ret, time_ret_prev;
            double R;
            if(if_approx){
                n = Vec3d(this->screen_potisions->data(j, k, 0), 
                    this->screen_potisions->data(j, k, 1), 
                    this->screen_potisions->data(j, k, 2));
                // ret time points
                time_ret = time - n.dot(position);
                time_ret_prev = time_prev - n.dot(position_prev);
                R = this->dmin[3];
            }
//#else
            else{
                n = Vec3d(this->screen_potisions->data(j, k, 0) * this->dmin[3] - position(0),
                    this->screen_potisions->data(j, k, 1) * this->dmin[3] - position(1),
                    this->screen_potisions->data(j, k, 2)* this->dmin[3] - position(2));
                time_ret = time + n.norm() - this->dmin[3];
                n = Vec3d(this->screen_potisions->data(j, k, 0) * this->dmin[3] - position_prev(0),
                    this->screen_potisions->data(j, k, 1) * this->dmin[3] - position_prev(1),
                    this->screen_potisions->data(j, k, 2) * this->dmin[3] - position_prev(2));
                R = n.norm();
                n.normalize();
                time_ret_prev = time_prev + R - this->dmin[3];
            }
//#endif
            double temp_left = this->time_det[0];
            double temp_right = this->time_det[this->nf[0] - 1];

            if ((time_ret_prev < temp_left) || (time_ret < temp_left) || 
            (time_ret > temp_right) || (time_ret_prev > temp_right)) continue;

            int idx_time = ((time_ret - this->time_det[0]) / this->d1);
            int idx_time_prev = ((time_ret_prev - this->time_det[0]) / this->d1);
            Vec3d far_field = charge * Lienard_Wiechert(n, beta_prev, beta_dot, R);

            // Interpolation for it
            for (int it = idx_time_prev; it < idx_time; it++) {
                double next_temp = (it + 1) * this->d1 + this->time_det[0]; 
                double temp = next_temp - time_ret_prev;
                //CHECK_BOUNDARY(it, j, k);
                this->emf->data(it, j, k) += temp * far_field[1];
                time_ret_prev = next_temp;
            }
            // Interpolation for tit
            double temp = time_ret - time_ret_prev;
            //CHECK_BOUNDARY(idx_time, j, k);
            this->emf->data(idx_time, j, k) += temp * far_field[1];

        }
    }
}

void SpheDetector::cmp_emf(Eigen::Ref<const Vec3dArr> position_arr, 
    Eigen::Ref<const Vec3dArr> position_prev_arr, Eigen::Ref<const Vec3dArr> beta_arr, 
    Eigen::Ref<const Vec3dArr> beta_prev_arr, double time, double charge, double dt)
{

    for(int i = 0; i < position_arr.rows(); i++){

        const Vec3d& beta=beta_arr.row(i);
        const Vec3d& beta_prev=beta_prev_arr.row(i);
        const Vec3d& position=position_arr.row(i);
        const Vec3d& position_prev=position_prev_arr.row(i);
        
        this->cmp_emf_single_particle(position, position_prev, beta, beta_prev,
            time, charge, dt);

#ifdef _OPENMP
        //int threadid = omp_get_thread_num();
        //cout << "nthreads: " << threadid << endl;
#endif

    }
    return;
}

void SpheDetector::reduce(){
#ifndef NONMPI
    if(faradio_mpi->rank() == 0){
        faradio_mpi->reduce(MPI_IN_PLACE, this->emf->get_data(), this->emf->size(), 0);
    }else{
        faradio_mpi->reduce(this->emf->get_data(), nullptr, this->emf->size(), 0);
    }
#endif
}
