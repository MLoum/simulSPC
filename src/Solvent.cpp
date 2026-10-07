#include "Solvent.h"
#include <cmath>

void Solvent::init(double viscosity,
                   double temperature,
                   double refractive_idx,
                   double width_nm,
                   double height_nm,
                   double mean_flow,
                   double grid_step_nm) {
    viscosity_Pa_s = viscosity;
    temperature_K = temperature;
    refractive_index = refractive_idx;
    width_y_nm = width_nm;
    height_z_nm = height_nm;
    mean_flow_uL_min = mean_flow;

    if (mean_flow_uL_min > 0.0) {
        compute_poiseuille_map(grid_step_nm);
    }
}

void Solvent::compute_poiseuille_map(double grid_step_nm) {
    const size_t ny = static_cast<size_t>(width_y_nm / grid_step_nm) + 1;
    const size_t nz = static_cast<size_t>(height_z_nm / grid_step_nm) + 1;

    // Centered coordinates: y in [-W/2, W/2], z in [-H/2, H/2]
    velocity_map_.init(ny, nz, -0.5 * width_y_nm, 0.5 * width_y_nm,
                               -0.5 * height_z_nm, 0.5 * height_z_nm);

    const double h_m = width_y_nm * 1.0e-9;
    const double l_m = height_z_nm * 1.0e-9;

    // Analytical series expansion for rectangular Poiseuille flow
    double sum_flow = (std::pow(h_m, 3) * l_m) / (12.0 * viscosity_Pa_s);
    for (int n = 1; n <= 30; ++n) {
        const double alpha_n = (2 * n - 1) * MathConsts::PI / h_m;
        const double term = (16.0 * std::pow(h_m, 4)) /
                            (std::pow(MathConsts::PI, 5) * viscosity_Pa_s * std::pow(2 * n - 1, 5));
        sum_flow -= term * ((std::cosh(alpha_n * l_m) - 1.0) / std::sinh(alpha_n * l_m));
    }

    // Volumetric flow Q: 1 uL/min = 1e-9 m^3 / 60 s
    const double flow_m3_s = (mean_flow_uL_min * 1.0e-9) / 60.0;
    const double grad_p = (sum_flow > 0.0) ? (flow_m3_s / sum_flow) : 0.0;

    for (size_t iy = 0; iy < ny; ++iy) {
        const double y_pos = -0.5 * width_y_nm + static_cast<double>(iy) * velocity_map_.dr();
        const double ym = (y_pos + 0.5 * width_y_nm) * 1.0e-9;

        for (size_t iz = 0; iz < nz; ++iz) {
            // Enforce exact Dirichlet no-slip boundary condition on channel walls
            if (iy == 0 || iy == ny - 1 || iz == 0 || iz == nz - 1) {
                velocity_map_.at(iy, iz) = 0.0;
                continue;
            }

            const double z_pos = -0.5 * height_z_nm + static_cast<double>(iz) * velocity_map_.dz();
            const double zm = (z_pos + 0.5 * height_z_nm) * 1.0e-9;

            double v = (grad_p * ym * (h_m - ym)) / (2.0 * viscosity_Pa_s);
            for (int n = 1; n <= 30; ++n) {
                const double alpha_n = (2 * n - 1) * MathConsts::PI / h_m;
                const double coeff = (4.0 * grad_p * h_m * h_m) /
                                     (viscosity_Pa_s * std::pow(MathConsts::PI, 3) * std::pow(2 * n - 1, 3));
                v -= coeff * ((std::sinh(alpha_n * zm) + std::sinh(alpha_n * (l_m - zm))) / std::sinh(alpha_n * l_m)) *
                     std::sin(alpha_n * ym);
            }
            velocity_map_.at(iy, iz) = std::max(0.0, v);
        }
    }
}