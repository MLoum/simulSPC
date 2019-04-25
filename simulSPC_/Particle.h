//
// Created by MatthieuL on 07/11/2018.
//

#ifndef SIMULSPC_PARTICLE_H
#define SIMULSPC_PARTICLE_H

#include "Experiment.h"
#include <random>

class Experiment;

class Particle
{
private :
    double x_, y_, z_; // in nm

    double r_hydro_;      // hydrodynamic radius in nm
    double Dx_, Dy_, Dz_;	// diffusion coefficients in nm^2 .ns-1

    double abs_cross_section_; // in m^-2
    double quantum_yield_;

    Experiment *exp_;

    std::default_random_engine generator_;

    std::normal_distribution<double> rand_walk_gaussian_;


public:
    Particle();

    Particle(double x_);

    Particle(Experiment &exp);
    ~Particle();
    void move();
    void light_matter_interraction();
    void inter_particle_interraction();

    //void set_exp(Experiment *exp) {exp_ = exp;}
    void initParticleParam(Experiment *exp);
};

#endif //SIMULSPC_PARTICLE_H
