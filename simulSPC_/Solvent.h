//
// Created by MatthieuL on 07/11/2018.
//

#ifndef SIMULSPC_SOLVENT_H
#define SIMULSPC_SOLVENT_H




class Particle;
class Experiment;

/**
 * \brief   Basic characterisitics of the medium where the particle are diffusing.
 *
 * \details
 * \sa	Experiment.h
 */
class Solvent {
public:
	Experiment *exp_;
	double viscosity_; //!<
	double temperature_; //!< in \b Kelvin
	double box_size_radial_, box_size_axial_; //!< in \b nanometer. Size of the box where to perform the simulation. We use periodic boundary condition.
	double refractive_index_;

	double w_microchanel_;
	double h_microchanel_;
	double meanFlow_;

	double** velocity_map_;

public:
	Solvent(Experiment *exp_);
	~Solvent();
	/**
 * \brief     Get the flow velocity in the microfluidic channel sensed by the particle. The speed vector has only one composant.
 * \details

	\return The x coordinate of the speed vector (y and z composant are zero).
 *
 */
	double get_convection(Particle *p);

	void create_velocity_map();

};


#endif //SIMULSPC_SOLVENT_H
