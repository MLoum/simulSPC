//
// Created by MatthieuL on 07/11/2018.
//

#include "Particle.h"
#include <stdio.h>
#include <gsl/gsl_randist.h>
//#include <random>

Particle::Particle(Experiment &exp) {

    exp_ = &exp;

    //TODO from File
    Dx_ = 5E-3;
    Dy_ = 5E-3;
    Dz_ = 5E-3;

    x_ = gsl_rng_uniform(exp_->rngGenerator_);
    y_ = gsl_rng_uniform(exp_->rngGenerator_);
    z_ = gsl_rng_uniform(exp_->rngGenerator_);
}

void Particle::initParticleParam(Experiment *exp)
{
    exp_ = exp;

    //TODO from File and gaussian distribution.
    Dx_ = 5E-3;
    Dy_ = 5E-3;
    Dz_ = 5E-3;

	abs_cross_section_ = 1E-16;
	quantum_yield_ = 0.95;

	//TODO precomputed value of abs_cross_section_/(h*nu)*quantum_yield_

    x_ = gsl_rng_uniform(exp_->rngGenerator_) * exp_->solvent_.box_size_radial_;
    y_ = gsl_rng_uniform(exp_->rngGenerator_) * exp_->solvent_.box_size_radial_;
    z_ = gsl_rng_uniform(exp_->rngGenerator_) * exp_->solvent_.box_size_radial_;
}

void Particle::move() {

    //TODO Dx * exp_->time_step_ change unit to spare a multiplciation
    x_ += Dx_ * gsl_ran_gaussian(exp_->rngGenerator_, 1) * exp_->time_step_;
    y_ += Dy_ * gsl_ran_gaussian(exp_->rngGenerator_, 1) * exp_->time_step_;
    z_ += Dz_ * gsl_ran_gaussian(exp_->rngGenerator_, 1) * exp_->time_step_;

    // printf("%f %f %f\n", x_, y_, z_);

    // Periodic boundary condition
    if(x_ < 0)
    {x_ += exp_->solvent_.box_size_radial_;}
    else if(x_ >= exp_->solvent_.box_size_radial_)
    {x_ -= exp_->solvent_.box_size_radial_;}

    if(y_ < 0)
    {y_ += exp_->solvent_.box_size_radial_;}
    else if(y_ >= exp_->solvent_.box_size_radial_)
    {y_ -= exp_->solvent_.box_size_radial_;}

    if(z_ < 0)
    {z_ += exp_->solvent_.box_size_axial_;}
    else if(z_ >= exp_->solvent_.box_size_axial_)
    {z_ -= exp_->solvent_.box_size_axial_;}

}

Particle::~Particle() {

}

void Particle::light_matter_interraction() {
	
	// 1) Based on particle position, angle, absorption cross-section and so on, calculate 
	//via the Molecular Detection Efficiency Function (MDF), the probability to detect a photon from the molecule.	
	double photon_probability = exp_->mdf_->get_detection_efficiency(this);

	// 2) Knowing the mean probability of counting a photon, use poissonian statistics to see if it is the case.
	if (gsl_ran_poisson(exp_->rngGenerator_, photon_probability) > 0)
	{
		// A photon was detected
		// Inform the exp structure to
		int dummy_photon = 1;

		//TODO microtime.
	}
		



	
}

void Particle::inter_particle_interraction(){

}


Particle::Particle()
{

}
