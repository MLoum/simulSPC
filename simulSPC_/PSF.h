//
// Created by MatthieuL on 07/11/2018.
//

#ifndef SIMULSPC_PSF_H
#define SIMULSPC_PSF_H


#include "OpticalSetup.h"

class PSF
{
    OpticalSetup OpticalSetup_;

public:
    PSF(OpticalSetup OpticalSetup_);
    ~PSF();
};





#endif //SIMULSPC_PSF_H
