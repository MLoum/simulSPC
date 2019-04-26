//
// Created by MatthieuL on 07/11/2018.
//

#ifndef SIMULSPC_SOLVENT_H
#define SIMULSPC_SOLVENT_H

/**
 * \brief   Basic characterisitics of the medium where the particle are diffusing.
 *
 * \details		
 * \sa	Experiment.h
 */
class Solvent {
public:
    double viscosity_; //!<
    double temperature_; //!< in \b Kelvin
    double box_size_radial_, box_size_axial_; //!< in \b nanometer. Size of the box where to perform the simulation. We use periodic boundary condition.
    double refractive_index_; 

public:
    Solvent();
    ~Solvent();
};


#endif //SIMULSPC_SOLVENT_H
