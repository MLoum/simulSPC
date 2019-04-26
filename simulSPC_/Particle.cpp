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

    //TODO from File
    Dx_ = 5E-3;
    Dy_ = 5E-3;
    Dz_ = 5E-3;

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
	
	// 1) Get PSF at particle position
	//exp_->psf_->get_E_field(*this);

	// 2) calculate the mean number of absorbed and then reemited photon during one tick of the experiment

	// 3) Calculate the probability is collected by the microscope via the CEF
	double photon_probability = 0.5;


	// 4) Knowing the mean probability of counting a photon, use poissonian statistics to see if it is the case.
	if (gsl_ran_poisson(exp_->rngGenerator_, photon_probability) > 0)
		int dummy_photon = 1;

	/*
	if (GL_diffAngulaire)
	{
		CoeffAbsorption = pow(sin(p->Theta)*cos(p->Phi), 4); //Absorption à 2 photons linéaire laser polarisé en x
		CoeffEmission = pow(sin(p->Theta), 2);;								// collection sans polariseur
		coeffDipoleScalaireChampElec = CoeffAbsorption * CoeffEmission; //2photons
	}
	else
		coeffDipoleScalaireChampElec = 1;

	if (OpticalSetup->ProfilExcitation[r][(int)(p->z)] == 0)
	{
		*NbreMoyenPhotonEmis = 0;	return;
	}

	up->IntensiteLaserenW * OpticalSetup->ProfilExcitation[r][(int)(p->z)];

	if (OpticalSetup->TypeExcitation == MONO_PHOTONIQUE)	
		*NbreMoyenPhotonEmis = ProfilIntensiteA_endroitParticule / OpticalSetup->EnergiePhoton * p->SectionEfficace * p->RendementQuantique * Param->pasDeTemps * 1e-9;
	*/

	
}

void Particle::inter_particle_interraction(){

}


Particle::Particle()
{

}
