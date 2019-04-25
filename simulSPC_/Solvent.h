//
// Created by MatthieuL on 07/11/2018.
//

#ifndef SIMULSPC_SOLVENT_H
#define SIMULSPC_SOLVENT_H


class Solvent {
public:
    double viscosity_; //SI
    double temperature_; //Kelvin
    double box_size_radial_, box_size_axial_; //in nm
    double refractive_index_;

public:
    Solvent();
    ~Solvent();
};


#endif //SIMULSPC_SOLVENT_H
