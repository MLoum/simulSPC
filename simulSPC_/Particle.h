//
// Created by MatthieuL on 07/11/2018.
//

#ifndef SIMULSPC_PARTICLE_H
#define SIMULSPC_PARTICLE_H

#include "Experiment.h"
#include "math.h"
//random -> still needed ?
#include <random>
#include <iostream>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>

class Experiment;

struct Monitor_variable;


/**
 * \brief      Class containing all the information on a particle (position, angular, quantum yield, ...)
 *
 * \details
 * \sa
 */
class Particle
{
public:
	double x_, y_, z_; //!< position in nanometer
	//TODO rotation matrix from the referential of the particle  to the referential of the laboratory.

	double Dx_, Dy_, Dz_;	//!<  rotational diffusion coefficients in radian^2.ns-1
	double Dphi1_, Dphi2_, Dphi3_;	//!<  translationnal diffusion coefficients in nm^2 .ns-1
	double r_hydro_;      //!<  hydrodynamic radius in nm
	double mean_delta_x_squared_, mean_delta_y_squared_, mean_delta_z_squared_;
	double mean_delta_phi1_squared_, mean_delta_phi2_squared_, mean_delta_phi3_squared_;

	double abs_cross_section_; //!<  in cm^-2
	double quantum_yield_;
	double brightness_; //!<  precomputed value of abs_cross_section_*quantum_yield_

	bool is_point_;
	bool is_isotropic_;

	Experiment *exp_;

	gsl_matrix* mat_transf_ = gsl_matrix_calloc(3, 3); //!<  Matrix to transform a vector described in the referential of the particle to the laboratory referential
	
	gsl_matrix* mat_rot1 = gsl_matrix_calloc(3, 3);
	gsl_matrix* mat_rot2 = gsl_matrix_calloc(3, 3);
	gsl_matrix* mat_rot3 = gsl_matrix_calloc(3, 3);
	gsl_matrix* mat_inter = gsl_matrix_calloc(3, 3);
	

	gsl_vector* abs_dipole_vector = gsl_vector_calloc(3);
	gsl_vector* abs_dipole_vector_lab = gsl_vector_calloc(3);
	gsl_vector* inter_vector = gsl_vector_calloc(3);

	std::vector<Monitor_variable> monitor_variable_vector;


public:
	Particle();

	Particle(double x_);

	Particle(Experiment &exp);
	~Particle();
	/**
 * \brief     Called at each tick of the main clock for each particle in order to translation and rotation
 * \details   	
 *
 */
	void move();

	/**
	 * \brief     Called at each tick of the main clock for each particle in order to know this particle has emitted a photon that is detected by the microscope
	 * \details   There are four steps.
				First get the electrical field at the position of the particle (via the calculated PSF).
				Second  calculate the mean number of absorbed and then reemited photon during one tick of the experiment
				Third, Based on the particle position and its angular orientation and its emission pattern, calculate (via the CEF) the probability that an eventual emitted photon is effectively collected by the microscope
				Fourth, Knowing the mean probability of counting a photon, use poissonian statistics to see if it is the case.

		\return Return TRUE if a photon emitted by this particle was detected
	 *
	 */
	bool light_matter_interraction();
	/**
	 * \brief    Called at each tick of the main clock for each particle. Simulate the interractions between close enough particles (like aggregation or FRET) during one tick the main clock
	 * \details   Empty for the time being

	 */
	void inter_particle_interraction();

	/**
	 * \brief    Conveniency function for multithreaded calculation.

	 */
	bool move_and_light_matter_interraction() {
		move(); 
		return  light_matter_interraction();
	}

	//void set_exp(Experiment *exp) {exp_ = exp;}
	/**
	 * \brief    Method called just aftr the instanciation of an array of particle and acting as a constructor.
				It initialize a particle based on the information gathered from the ini file.
	*	\sa Experiment
	 */
	void initParticleParam(Experiment* exp, int ligne_start, int ligne_end);

	/**
 * \brief    Return the absorption dipole in the ->referential of the lab<-
			
 */
	gsl_vector* get_rotated_abs_dipole_vector();

	//void set_exp(Experiment *exp) {exp_ = exp;}
	/**
	 * \brief    Export in a text file the translational and angular position of a particle in order to test its adequacy with a brownian motion
	*	\sa Experiment
	 */
	void brownianMotionTest();

	
	void display_particle_position();

	void monitor_particle(int mode=1);

};

struct Monitor_variable
{
	//double x_, y_, z_;
	double r_, z_;
	//double phi_x, phi_y, phi_z;
	double e_dot_mu;
	// add other variable if needed
};

#endif //SIMULSPC_PARTICLE_H
