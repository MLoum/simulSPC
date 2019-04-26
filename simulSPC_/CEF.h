#pragma once

#include "OpticalSetup.h"
#include "Experiment.h"
#include "Particle.h"

class Experiment;
class OpticalSetup;

/**
 * \brief      Abstract class for representing the Collection Efficiency Function (CEF)
 *
 * \details		A Collection Efficiency Function (CEF) gives the probability of collecting a photon based on the position and the angular orientation and the emission pattern (isotropic, dipole, ...) of a particle
 * \sa
 */
class CEF
{
protected:
	OpticalSetup *OpticalSetup_;
	Experiment *exp_;

public:
	/**
	 * \brief    Pure virtual method
	 * \param    Particle &p  Particle with which to calculate the CEF
	 * \return   A double with the probability ([0,1]) of detecting a photon at the position and with the orientation of the p particle.
	 */
	virtual double getCollectionProbability(Particle &p) = 0;
	CEF();
	~CEF();
};

