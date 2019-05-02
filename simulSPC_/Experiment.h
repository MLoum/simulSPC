//
// Created by MatthieuL on 07/11/2018.
//

#ifndef SIMULSPC_EXPERIMENT_H
#define SIMULSPC_EXPERIMENT_H

/**
*  \file Experiment.h
 * \version   1.0
 * \date       2018 11 07
 * \brief      Main structure of the programm
 *
 * \details    This class gather all the instance of the parameters (time_step_, space_step_, OpticalSetup, solvent, PSF, CEF, ...) and the list of the particle particle_list. 
 *	The main function is "tick" who simulate the experience for one tick of the main clock.
 * Some members are for convenience (like in python)
 * \copyright Moltech-Anjou UMR6200 Samson Team
 */

class MDF;
#include "OpticalSetup.h"
#include "Solvent.h"
#include "Particle.h"
#include "MDF.h"
#include <gsl/gsl_randist.h>


class OpticalSetup;
class Solvent;
class Particle;
//class PSF;
//class CEF;

//class MDF_gaussian;

//typedef long long __int64;
//typedef unsigned long long __uint64;



/**
 * \brief      Main structure of the programm
 *
 * \details    This class gather all the instance of the parameters (time_step_, space_step_, OpticalSetup, solvent, PSF, CEF, ...) and the list of the particle particle_list.
 *	The main function is "tick" who simulate the experience for one tick of the main clock.
 * Some members are for convenience (like in python)
 * \sa
 */
class Experiment{

public:
	friend class Particle;
    double time_step_; //!< Time step of the simulation in \b nanosecond, typical value 50ns
    double space_step_; //!< Space step in \b nanometer for the PSF and CEF if they are pre-computed, typical value 50nm (lambda/10)
    Solvent solvent_; //!< See the corresponding class
    OpticalSetup opticalSetup_; //!< See the corresponding class
	//PSF *psf_; //!< See the corresponding class
	//CEF *cef_; //!< See the corresponding class
	MDF *mdf_; //!< See the corresponding class

	gsl_rng *rngGenerator_; //!< GSL algorithm for random number generation.

    int nb_of_particle_;

    Particle *particle_list; //!< Array of particle
    long long *photon_list; //!< Array of the time arrivals (in clock tick unit) of the collected photon. This is the main output of the software.

public:
	/**
	 * \brief       Main method of the Experiment class. It computes one tick of the experiment (particle movement, particle interaction, ligth and matter interraction, ...)
	 * \details    Called in the int main void of the software.
	 */
    void tick();

    Experiment();
    ~Experiment();

};


#endif //SIMULSPC_EXPERIMENT_H
