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
#include <vector>
#include <fstream>
#include <string>
#include <iostream>
#include <algorithm>
#include "Solvent.h"
#include "Particle.h"
#include "MDF.h"
#include "Photon.h"
#include <gsl/gsl_randist.h>

using namespace std;

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

	__int64 macroClock_; //!< Clock for the experiment. At each tick it increase by one. If a photon is detected its arrival time is equal to macroClock_;

    Particle *particle_list_; //!< Array of particle
	std::vector<Photon> photon_vector_;   //!< Vector of the generated photon (time arrivals (in clock tick unit), microtime, channel). This is the main output of the software.
    long long *photon_array_; 

	vector<string> iniFilevector_;

public:
	/**
	 * \brief       Main method of the Experiment class. It computes one tick of the experiment (particle movement, particle interaction, ligth and matter interraction, ...)
	 * \details    Called in the int main void of the software.
	 */
    void tick();

    Experiment();
    ~Experiment();
	void convertPhoton_vectorToList();
	void write_photon_vector();

	void read_ini_file();
	double init_parameter(std::string parameter, char start, char end);
};


#endif //SIMULSPC_EXPERIMENT_H
