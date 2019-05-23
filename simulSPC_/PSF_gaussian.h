#pragma once
#include "PSF.h"
#include "Experiment.h"

/**
 * \brief      Simple Gaussian-Lorentzian PSF
 *
 * \details		This PSF doesn't take in account the vectorial nature of the electrical field. Stated differently, the polarisation is uniform. The PSF is precalculated using the fact that there is a  symetry of revolution of the PSF around the optical axis.
 * \sa	PSF.h
 */
class PSF_gaussian :
	public PSF
{
public:
	double **psf_;	
	
	//E_field get_E_field(Particle &p);

	PSF_gaussian();
	PSF_gaussian(Experiment *exp);
	virtual double get_probability_absorption(Particle *p);
	virtual void explore_psf();
	~PSF_gaussian();

};

