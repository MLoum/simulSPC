#pragma once

#include "Math3D.h"
#include <cmath>
#include <algorithm>

struct Laser {
    double wavelength_nm{488.0};     // Excitation wavelength in nm
    double power_mW{0.01};            // Laser power at focal plane (10 uW = 0.01 mW)
    Vec3 polarization{1.0, 0.0, 0.0}; // Linear polarization unit vector in laboratory frame
};

struct Objective {
    double NA{1.2};                  // Numerical Aperture
    double magnification{60.0};      // Objective magnification
    double immersion_index{1.333};   // Refractive index of immersion medium

    /**
     * \brief Fraction of 4*pi steradians intercepted by the lens: (1 - cos(theta_max)) / 2
     */
    double geometric_collection_efficiency() const {
        const double sin_theta = std::min(1.0, NA / immersion_index);
        const double cos_theta = std::sqrt(1.0 - sin_theta * sin_theta);
        return 0.5 * (1.0 - cos_theta);
    }
};

struct DetectionChannel {
    double pinhole_diameter_um{50.0}; // Confocal pinhole diameter in um
    double optics_transmission{0.70};  // Dichroic and emission bandpass filters transmission
    double pinhole_efficiency{0.70};   // Geometric/diffractive transmission through pinhole
    double detector_efficiency{0.65};  // SPAD / APD photon detection efficiency (PDE)
    double dark_count_rate_cps{100.0}; // Detector dark counts per second
};

/**
 * \brief Consolidated optical parameters for excitation and detection.
 */
class OpticalSetup {
public:
    Laser laser;
    Objective objective;
    DetectionChannel detection;

    OpticalSetup() = default;

    /**
     * \brief Total photon detection efficiency: eta_geom * T_optics * T_pinhole * PDE
     */
    double total_detection_efficiency() const {
        return objective.geometric_collection_efficiency() *
               detection.optics_transmission *
               detection.pinhole_efficiency *
               detection.detector_efficiency;
    }
};