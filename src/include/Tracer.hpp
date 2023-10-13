#ifndef PARTICLE_H
#define PARTICLE_H
#include <vector>
#include <iostream>
#include "Vector.hpp"

using namespace std;

class Particle
{
public:
    Particle(vector<vector<double>> positions, vector<vector<double>> betas,
        vector<vector<double>> position_prev, vector<vector<double>> beta_prev,  
        double time, double dt, double charge, int iter_id);
    ~Particle();

    vector<R_vec>& get_position(){return _position;};
    vector<R_vec>& get_beta(){return _beta;};
    double get_time(){return _time;};
    vector<R_vec>& get_position_prev(){return _position_prev;};
    vector<R_vec>& get_beta_prev(){return _beta_prev;};
    double get_charge(){return _charge;};
    double get_dt(){return _dt;};
    int get_iter_id(){return _iter_id;};

protected:
    vector<R_vec> _position;
    vector<R_vec> _position_prev;
    vector<R_vec> _beta;
    vector<R_vec> _beta_prev;

    double _time;
    double _dt;
    double _charge;
    int _iter_id;
};

#endif