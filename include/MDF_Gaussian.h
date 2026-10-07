#pragma once

#include "MDF.h"

/**
 * \brief Gaussian-Gaussian approximation of the confocal Molecular Detection Efficiency Function.
 * \details Calibrated in rigorous SI physical units for single-molecule photon counting.
 */
class MDF_Gaussian : public MDF {
private:
    OpticalSetup optical_setup_;
    Grid2D mdf_grid_;

    double w_xy_nm_{250.0};            // Radial 1/e^2 beam waist in nm
    double w_z_nm_{1000.0};            // Axial 1/e^2 beam waist in nm
    double peak_flux_density_m2_s_{0.0}; // Peak photon flux density S0 [photons / (s * m^2)]

public:
    MDF_Gaussian() = default;

    /**
     * \brief Initializes Gaussian MDF grid and physical photon calibration.
     * \param opt Optical setup specification.
     * \param aspect_ratio Empirical axial-to-radial waist ratio wz / wxy (typically 3 to 5).
     * \param box_radial_nm Radial simulation domain boundary.
     * \param box_axial_nm Axial simulation domain boundary.
     * \param space_step_nm Grid discretization step in nm.
     */
    void init(const OpticalSetup& opt,
              double aspect_ratio = 4.0,
              double box_radial_nm = 1000.0,
              double box_axial_nm = 2000.0,
              double space_step_nm = 25.0);

    double get_detection_probability(const Particle& p, double dt_ns) const override;

    const Grid2D& get_grid() const override { return mdf_grid_; }

    double w_xy_nm() const { return w_xy_nm_; }
    double w_z_nm() const { return w_z_nm_; }
    double peak_flux_density() const { return peak_flux_density_m2_s_; }
};