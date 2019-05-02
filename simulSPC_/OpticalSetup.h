//
// Created by MatthieuL on 07/11/2018.
//

#ifndef SIMULSPC_OPTICALSETUP_H
#define SIMULSPC_OPTICALSETUP_H

/**
 * \brief     Excitation Laser
 *
 * \sa
 */
struct laser_exc
{
    double wl; //!< in nanometer
    double intensity; //!< in mW. Laser intensity, that is to say the laser power integrated on all the the beam surface 
	//TODO polarisation
};


/**
 * \brief     Microscope objective
 *
 * \sa
 */
struct objective
{
    double NA;
    double magnification;
    double focal_length; //!< in millimeter
	double tube_lens_focal_length; //!< in centimeter, typically 18cm for Olympus microscope (20cm -> Nikon, 16cm-> Zeiss)
};

/**
 * \brief      Avalanche Photodiode module
 *
 * \sa
 */
struct APD
{
    double after_pulsing_probability;
    double deadtime; // in ns
};

/**
 * \brief      Collection of small structure containing all the paramters of the experiments
 *
 * \details    This class gather all the parameter of the optical setup for simulating the experience like for instance the laser wavelength or the numerical aperture of the objective
 * \sa
 */
class OpticalSetup
{
public:

    laser_exc laser_exc_;
    objective objective_;
    double confocal_pinhole_diameter; //!< in micrometer
    APD APD_;

    OpticalSetup();
    ~OpticalSetup();
};


#endif //SIMULSPC_OPTICALSETUP_H
