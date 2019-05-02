//
// Created by MatthieuL on 07/11/2018.
//

#include "Experiment.h"
#include "MDF_gaussian.h"
#include <time.h>

void Experiment::tick() {

    for(int i=0; i<nb_of_particle_; i++)
    {
        particle_list[i].move();
//        particle_list[i].inter_particle_interraction();
        particle_list[i].light_matter_interraction();
        // Photon probability -> branching
		//TODO 
    }

	//photon_vector_.push_back(photon);
}

Experiment::~Experiment() {
delete[] particle_list;
delete[] photon_list;
}

void Experiment::convertPhoton_vectorToList()
{
	if (photon_list != NULL)
		delete[] photon_list;
	photon_list = new __int64[photon_vector_.size()];

	for (int i = 0; i < photon_vector_.size(); i++)
		photon_list[i] = photon_vector_[i].compress_data();
}

Experiment::Experiment()  {

	// choosing and random number generator... at "random" for now e.g. gsl_rng_default gsl_rng_taus	 gsl_rng_mt19937
	rngGenerator_ = gsl_rng_alloc(gsl_rng_taus);
	// rng seed
	gsl_rng_set(rngGenerator_, time(NULL));

    //TODO read from a txt file.
    // Dummy default values
    nb_of_particle_ = 10;
    time_step_ = 50;
    space_step_ = 100;

    solvent_.box_size_radial_ = 1E5;
    solvent_.box_size_axial_ = 5E5;

	opticalSetup_.laser_exc_.wl = 405;
	opticalSetup_.laser_exc_.intensity = 1;
	opticalSetup_.objective_.NA = 0.95;

    particle_list = new Particle[nb_of_particle_];
    for(int i=0; i<nb_of_particle_; i++)
        particle_list[i].initParticleParam(this);

    photon_list = new __int64[1000000];

    // create PSF, CEF and MDF
	// TODO read from file the type of MDF.
	mdf_ = new MDF_gaussian(this);
    
}
