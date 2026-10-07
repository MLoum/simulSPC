#pragma once

#include "Math3D.h"
#include "FastRNG.h"

namespace PhysicsConsts {
    constexpr double KB = 1.380649e-23;
    constexpr double NA = 6.02214076e23;
    constexpr double C_LIGHT = 299792458.0;
    constexpr double H_PLANCK = 6.62607015e-34;
}

class Particle {
public:
    Vec3 pos_{0.0, 0.0, 0.0};
    Quaternion orientation_{1.0, 0.0, 0.0, 0.0};
    Vec3 abs_dipole_local_{1.0, 0.0, 0.0};
    Vec3 box_size_{1000.0, 1000.0, 2000.0};

    double r_hydro_nm_{10.0};
    double D_trans_nm2_ns_{0.0};
    double D_rot_ns_{0.0};
    double sigma_trans_nm_{0.0};
    double sigma_rot_rad_{0.0};

    double lifetime_ns_{3.8};
    double abs_cross_section_cm2_{1.0e-16};
    double quantum_yield_{1.0};

    bool is_point_{true};
    bool is_isotropic_{false};

public:
    Particle() = default;

    void init_physics(double r_hydro_nm,
                      double temperature_K,
                      double viscosity_Pa_s,
                      double dt_ns,
                      double lifetime_ns,
                      double abs_cross_section_cm2,
                      double quantum_yield,
                      bool is_isotropic);

    void set_box_dimensions(double box_radial_nm, double box_axial_nm);
    void set_local_dipole(const Vec3& dipole);

    Vec3 get_lab_dipole() const {
        return orientation_.rotate(abs_dipole_local_);
    }

    void move(FastRNG& rng, const Vec3& drift_nm = Vec3{0.0, 0.0, 0.0});
    bool detect_photon(FastRNG& rng, double detection_probability, double threshold = 1.0e-7);

    /**
     * \brief Tire le canal (0: parallele, 1: perpendiculaire) selon la projection du dipôle.
     */
    uint8_t sample_polarization_channel(FastRNG& rng) const;

    /**
     * \brief Tire le microtemps TCSPC avec déclin exponentiel et convolution IRF.
     */
    uint16_t sample_microtime(FastRNG& rng, double time_step_ns, double irf_sigma_ns) const;

private:
    void apply_periodic_boundaries();
};