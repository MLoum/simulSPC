//
// Created by MatthieuL on 07/11/2018.
//

#ifndef SIMULSPC_PSF_H
#define SIMULSPC_PSF_H


#include "OpticalSetup.h"
#include "Experiment.h"
#include "Particle.h"
#include "E_field.h"

class Experiment;
class OpticalSetup;
class Particle;
/**
 * \brief      Abstract class for representing the Point Spread Function (PSF)
 *
 * \details		A Point Spread Function (PSF) gives the value of the electrical field of a laser focused by an high numerical apperture microscope objective.
 * \sa
 */
class PSF
{
protected:
    OpticalSetup *OpticalSetup_;
	Experiment *exp_;

public:
	PSF();
    PSF(OpticalSetup OpticalSetup_);
    ~PSF();

	/**
	 * \brief    Pure virtual method
	 * \param    Particle &p  Particle with which to calculate the PSF
	 * \return   The laser electrical field at the position of the particle
	 */
	virtual E_field get_E_field(Particle &p) = 0;

};





#endif //SIMULSPC_PSF_H
