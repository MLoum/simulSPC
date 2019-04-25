//
// Created by MatthieuL on 07/11/2018.
//

#ifndef SIMULSPC_OPTICALSETUP_H
#define SIMULSPC_OPTICALSETUP_H


struct laser_exc
{
    double wl; //in nm
    double intensity; // in W
};

struct objective
{
    double NA;
    double magnification;
    double focal_length; //in mm
};

struct APD
{
    double after_pulsing_probability;
    double deadtime; // in ns
};

class OpticalSetup
{
public:

    laser_exc laser_exc_;

    objective objective_;

    double confocal_pinhole_diameter; //in µm

    APD APD_;

    OpticalSetup();
    ~OpticalSetup();
};


#endif //SIMULSPC_OPTICALSETUP_H
