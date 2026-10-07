#include "MDF_Gaussian.h"
#include <cmath>

void MDF_Gaussian::init(const OpticalSetup& opt,
                        double aspect_ratio,
                        double box_radial_nm,
                        double box_axial_nm,
                        double space_step_nm) {
    optical_setup_ = opt;

    // Standard confocal lateral waist: w_xy ~ 0.66 * (1.22 * lambda / (2 * NA))
    const double airy_radius_nm = (1.22 * opt.laser.wavelength_nm) / (2.0 * opt.objective.NA);
    w_xy_nm_ = airy_radius_nm * 0.66;
    w_z_nm_ = w_xy_nm_ * aspect_ratio;

    // --- Absolute Physical Photon Flux Calibration (SI units) ---
    // Laser power: P [W] = power_mW * 1e-3
    const double power_w = opt.laser.power_mW * 1.0e-3;
    // Photon energy: E_ph [J] = h * c / lambda
    const double lambda_m = opt.laser.wavelength_nm * 1.0e-9;
    const double photon_energy_j = (PhysicsConsts::H_PLANCK * PhysicsConsts::C_LIGHT) / lambda_m;
    // Total photon rate: Phi [photons / s] = P / E_ph
    const double photon_rate_s = power_w / photon_energy_j;

    // Focal peak intensity flux density for TEM00 Gaussian beam: S0 = 2 * Phi / (pi * w_xy^2)
    const double w_xy_m = w_xy_nm_ * 1.0e-9;
    const double beam_area_m2 = 0.5 * MathConsts::PI * w_xy_m * w_xy_m;
    peak_flux_density_m2_s_ = photon_rate_s / beam_area_m2;

    // --- Tabulate Spatial Profile on Cylindrical Grid (r in [0, R_max], z in [-Z_max/2, Z_max/2]) ---
    const double max_r_nm = std::sqrt(2.0) * box_radial_nm * 0.5;
    const double half_z_nm = box_axial_nm * 0.5;

    const size_t nr = static_cast<size_t>(max_r_nm / space_step_nm) + 2;
    const size_t nz = static_cast<size_t>(box_axial_nm / space_step_nm) + 2;

    mdf_grid_.init(nr, nz, 0.0, max_r_nm, -half_z_nm, half_z_nm);

    const double inv_wxy2 = 1.0 / (w_xy_nm_ * w_xy_nm_);
    const double inv_wz2 = 1.0 / (w_z_nm_ * w_z_nm_);

    for (size_t ir = 0; ir < nr; ++ir) {
        const double r = ir * mdf_grid_.dr();
        const double r_term = -2.0 * r * r * inv_wxy2;
        for (size_t iz = 0; iz < nz; ++iz) {
            const double z = -half_z_nm + iz * mdf_grid_.dz();
            const double z_term = -2.0 * z * z * inv_wz2;
            mdf_grid_.at(ir, iz) = std::exp(r_term + z_term);
        }
    }
}

double MDF_Gaussian::get_detection_probability(const Particle& p, double dt_ns) const {
    const double r_nm = std::sqrt(p.pos_.x * p.pos_.x + p.pos_.y * p.pos_.y);
    const double spatial_weight = mdf_grid_.interpolate(r_nm, p.pos_.z);

    if (spatial_weight <= 0.0) {
        return 0.0;
    }

    // Dipole orientation factor theta_orient
    double theta_orient = 1.0;
    if (!p.is_isotropic_) {
        const Vec3 mu = p.get_lab_dipole();
        const double cos_alpha = optical_setup_.laser.polarization.dot(mu);
        // Factor 3 ensures ensemble average <3 * cos^2(alpha)> = 1 over isotropic distribution
        theta_orient = 3.0 * cos_alpha * cos_alpha;
    }

    // Absorption rate: k_abs [s^-1] = S0 [m^-2 s^-1] * sigma_abs [m^2] * spatial_weight * theta_orient
    // Convert cross section: 1 cm^2 = 1e-4 m^2
    const double sigma_abs_m2 = p.abs_cross_section_cm2_ * 1.0e-4;
    const double k_abs_s = peak_flux_density_m2_s_ * sigma_abs_m2 * spatial_weight * theta_orient;

    // Detected photon rate: R_det = k_abs * quantum_yield * total_detection_efficiency
    const double r_det_s = k_abs_s * p.quantum_yield_ * optical_setup_.total_detection_efficiency();

    // Expected photon count over time-step dt: lambda = R_det [s^-1] * dt [ns] * 1e-9
    return r_det_s * (dt_ns * 1.0e-9);
}