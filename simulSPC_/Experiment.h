//
// Created by MatthieuL on 07/11/2018.
//

#ifndef SIMULSPC_EXPERIMENT_H
#define SIMULSPC_EXPERIMENT_H

#include "OpticalSetup.h"
#include "Solvent.h"
#include "Particle.h"
#include "PSF.h"
#include "MDF.h"


class OpticalSetup;
class Solvent;
class Particle;
class PSF;
class MDF;

//typedef long long __int64;
//typedef unsigned long long __uint64;

class Experiment{

public:
    double time_step_; // in ns
    double space_step_; //in nm
    Solvent solvent_;
    OpticalSetup opticalSetup_;

    int nb_of_particle_;

    Particle *particle_list;
    long long *photon_list;

    void tick();

    Experiment();
    ~Experiment();


};


#endif //SIMULSPC_EXPERIMENT_H
