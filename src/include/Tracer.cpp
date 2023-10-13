#include "Tracer.hpp"

Particle::Particle(vector<vector<double>> positions, vector<vector<double>> betas,
        vector<vector<double>> position_prev, vector<vector<double>> beta_prev,  
        double time, double dt, double charge, int iter_id):
     _time(time), _dt(dt), _charge(charge), _iter_id(iter_id)
{
    int particle_nums = positions.size();
    for(int i = 0; i < particle_nums; i++){
        _position.push_back(R_vec(positions[i][0], positions[i][1], positions[i][2]));
        _beta.push_back(R_vec(betas[i][0], betas[i][1], betas[i][2]));
        _beta_prev.push_back(R_vec(beta_prev[i][0], beta_prev[i][1], beta_prev[i][2]));
        _position_prev.push_back(R_vec(position_prev[i][0], position_prev[i][1], position_prev[i][2]));
    }
    cout << "Address:" << hex << long(this) << " Particle constructor called" << endl;
};

Particle::~Particle()
{
    cout << "Address:" << hex << long(this) << " Particle destructor called" << endl;
};