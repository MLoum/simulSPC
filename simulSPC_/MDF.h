//
// Created by MatthieuL on 07/11/2018.
//

#ifndef SIMULSPC_MDF_H
#define SIMULSPC_MDF_H

#include "Experiment.h"

class Experiment;

class MDF {
public:
    Experiment *exp_;

    double **mdf_;
    int nb_hpoint_, nb_vpoint_;
    double airy_diameter_; // Airy disk diameter -> 1.22 * lamdba/NA
    double w_xy_, w_z_; // radial and axial waist

    virtual ~MDF();

    MDF(Experiment &exp);

    //MDF(Experiment *exp);

    void create_PSF();
    void create_CEF();
};


#endif //SIMULSPC_MDF_H
