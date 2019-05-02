//
// Created by MatthieuL on 07/11/2018.
//

#ifndef SIMULSPC_PARTICLE_H
#define SIMULSPC_PARTICLE_H

#include "Experiment.h"
#include <random>

class Experiment;

/**
 * \brief      Class containing all the information on a particle (position, angular, quantum yield, ...)
 *
 * \details    
 * \sa
 */
class Particle
{
public :
    double x_, y_, z_; //!< position in nanometer
	double phi1_, phi2_, phi_3; 

    
    double Dx_, Dy_, Dz_;	//!<  translationnal diffusion coefficients in nm^2 .ns-1
	double r_hydro_;      //!<  hydrodynamic radius in nm

    double abs_cross_section_; //!<  in cm^-2
    double quantum_yield_;

    Experiment *exp_;	


public:
    Particle();

    Particle(double x_);

    Particle(Experiment &exp);
    ~Particle();
    void move();

	/**
	 * \brief     Called at each tick of the main clock for each particle in order to know this particle has emitted a photon that is detected by the microscope
	 * \details   There are four steps. 
				First get the electrical field at the position of the particle (via the calculated PSF).
				Second  calculate the mean number of absorbed and then reemited photon during one tick of the experiment 
				Third, Based on the particle position and its angular orientation and its emission pattern, calculate (via the CEF) the probability that an eventual emitted photon is effectively collected by the microscope 
				Fourth, Knowing the mean probability of counting a photon, use poissonian statistics to see if it is the case.

	 *
	 */
    void light_matter_interraction();
	/**
	 * \brief    Called at each tick of the main clock for each particle. Simulate the interractions between close enough particles (like aggregation or FRET) during one tick the main clock
	 * \details   Empty for the time being

	 */
    void inter_particle_interraction();

    //void set_exp(Experiment *exp) {exp_ = exp;}
	/**
	 * \brief    Method called just after the instanciation of an array of particle and acting as a constructor.
	*	\sa Experiment
	 */
    void initParticleParam(Experiment *exp);
};

#endif //SIMULSPC_PARTICLE_H
