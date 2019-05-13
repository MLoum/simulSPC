//
// Created by MatthieuL on 07/11/2018.
//

#include "Particle.h"
#include <stdio.h>
#include <gsl/gsl_randist.h>
#define PI 3.1415159265
#define kb 1.38E-23

Particle::Particle(Experiment &exp) {

    exp_ = &exp;

    x_ = gsl_rng_uniform(exp_->rngGenerator_);
    y_ = gsl_rng_uniform(exp_->rngGenerator_);
    z_ = gsl_rng_uniform(exp_->rngGenerator_);
}

void Particle::initParticleParam(Experiment *exp, int ligne_start, int ligne_end)
{
    exp_ = exp;

	double r_hydro_max, r_hydro_sigma;

	r_hydro_max = exp_->init_parameter("r_hydro=", ligne_start, ligne_end);
	r_hydro_sigma = exp_->init_parameter("r_hydro=", ligne_start, ligne_end);

	r_hydro_ = gsl_ran_gaussian(exp_->rngGenerator_, r_hydro_sigma * r_hydro_max / 100) + r_hydro_max;

	abs_cross_section_ = exp_->init_parameter("cross_section=", ligne_start, ligne_end);
	quantum_yield_ = exp_->init_parameter("quantum_yield=",  ligne_start, ligne_end);

	Dx_ = 1E9 * (kb * exp_->solvent_.temperature_) / (6 * PI * exp_->solvent_.viscosity_ * r_hydro_ * 1E-9);

	//TODO precomputed value of abs_cross_section_/(h*nu)*quantum_yield_

    x_ = gsl_rng_uniform(exp_->rngGenerator_) * exp_->solvent_.box_size_radial_;
    y_ = gsl_rng_uniform(exp_->rngGenerator_) * exp_->solvent_.box_size_radial_;
    z_ = gsl_rng_uniform(exp_->rngGenerator_) * exp_->solvent_.box_size_axial_;
}

void Particle::brownianMotionTest()
{
	//TODO Copy the first version
}

void Particle::move() {

    //TODO Dx * exp_->time_step_ change unit to spare a multiplciation
    x_ += Dx_ * gsl_ran_gaussian(exp_->rngGenerator_, 1) * exp_->time_step_;
    y_ += Dy_ * gsl_ran_gaussian(exp_->rngGenerator_, 1) * exp_->time_step_;
    z_ += Dz_ * gsl_ran_gaussian(exp_->rngGenerator_, 1) * exp_->time_step_;

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

bool Particle::light_matter_interraction() {
	
	// 1) Based on particle position, angle, absorption cross-section and so on, calculate 
	//via the Molecular Detection Efficiency Function (MDF), the probability to detect a photon from the molecule.	
	double photon_probability = exp_->mdf_->get_detection_efficiency(this);

	// 2) Knowing the mean probability of counting a photon, use poissonian statistics to see if it is the case.
	if (gsl_ran_poisson(exp_->rngGenerator_, photon_probability) > 0)
	{
		// A photon was detected
		//TODO microtime.
		return true;
	}
	else
		return false;	
}

void Particle::inter_particle_interraction(){

}


Particle::Particle()
{

}

/* Debug function that should disappear at some point*/
void Particle::display_particle_position()
{
	cout << x_ << "; " << y_ << "; " << z_ << "; " << endl;
}
