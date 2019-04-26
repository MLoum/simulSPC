#pragma once
#include "PSF.h"
#include "Experiment.h"

/**
 * \brief      Simple Gaussian-Gaussian PSF
 *
 * \details		This PSF doesn't take in account the vectorial nature of the electrical field. The PSF is precalculated using the fact that there is a  symetry of revolution of the PSF around the optical axis.
 * \sa	PSF.h
 */
class PSF_gaussian :
	public PSF
{
public:
	double **psf_;
	
	E_field get_E_field(Particle &p);

	PSF_gaussian();
	PSF_gaussian(OpticalSetup OpticalSetup_);
	~PSF_gaussian();

};

