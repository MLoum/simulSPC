//
// Created by MatthieuL on 07/11/2018.
//

#include "Experiment.h"

void Experiment::tick() {

    for(int i=0; i<nb_of_particle_; i++)
    {
        particle_list[i].move();
//        particle_list[i].inter_particle_interraction();
        particle_list[i].light_matter_interraction();
        // Photon probability -> branching
    }
}

Experiment::~Experiment() {
delete(particle_list);
delete(photon_list);
}

Experiment::Experiment()  {

    //TODO read from a txt file.
    // Dummy defulat values
    nb_of_particle_ = 10;
    time_step_ = 50;
    space_step_ = 100;

    solvent_.box_size_radial_ = 1E5;
    solvent_.box_size_axial_ = 5E5;

    particle_list = new Particle[nb_of_particle_];
    for(int i=0; i<nb_of_particle_; i++)
        particle_list[i].initParticleParam(this);

    photon_list = new __int64[1000000];

    // create PSF, CEF and MDR
    
}
