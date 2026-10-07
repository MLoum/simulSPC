#pragma once

#include "Math3D.h"
#include "Grid2D.h"

/**
 * \brief Fluid medium properties and microfluidic Hagen-Poiseuille velocity mapping.
 */
class Solvent {
public:
    double viscosity_Pa_s{1.0e-3};   // Dynamic viscosity in Pa.s (water ~ 1e-3 Pa.s)
    double temperature_K{293.15};     // Absolute temperature in K
    double refractive_index{1.333};   // Solvent refractive index

    // Microchannel cross-section geometry in nm
    double width_y_nm{1000.0};       // Channel dimension along Y
    double height_z_nm{2000.0};      // Channel dimension along Z
    double mean_flow_uL_min{0.0};    // Volumetric flow rate in uL/min

    Grid2D velocity_map_;             // Precomputed convective velocity map vx(y, z) in m/s

public:
    Solvent() = default;

    void init(double viscosity,
              double temperature,
              double refractive_idx,
              double width_nm,
              double height_nm,
              double mean_flow_uL_min = 0.0,
              double grid_step_nm = 25.0);

    /**
     * \brief Computes convective displacement in nanometers over duration dt_ns.
     */
    inline Vec3 get_drift_nm(const Vec3& pos, double dt_ns) const {
        if (mean_flow_uL_min <= 0.0) {
            return Vec3{0.0, 0.0, 0.0};
        }
        // Convective velocity along X based on local (y, z) position
        const double vx_m_s = velocity_map_.interpolate(pos.y, pos.z);
        // Conversion: vx [m/s] * dt [ns] = vx [nm/ns] * dt [ns] = vx * dt [nm]
        return Vec3{vx_m_s * dt_ns, 0.0, 0.0};
    }

private:
    void compute_poiseuille_map(double grid_step_nm);
};