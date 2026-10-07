#include "Particle.h"
#include <cmath>
#include <algorithm>

void Particle::init_physics(double r_hydro_nm,
                            double temperature_K,
                            double viscosity_Pa_s,
                            double dt_ns,
                            double lifetime_ns,
                            double abs_cross_section_cm2,
                            double quantum_yield,
                            bool is_isotropic) {
    r_hydro_nm_ = r_hydro_nm;
    lifetime_ns_ = lifetime_ns;
    abs_cross_section_cm2_ = abs_cross_section_cm2;
    quantum_yield_ = quantum_yield;
    is_isotropic_ = is_isotropic;

    const double r_m = r_hydro_nm * 1.0e-9;
    const double D_t_si = (PhysicsConsts::KB * temperature_K) / (6.0 * MathConsts::PI * viscosity_Pa_s * r_m);
    D_trans_nm2_ns_ = D_t_si * 1.0e9;
    sigma_trans_nm_ = std::sqrt(2.0 * D_trans_nm2_ns_ * dt_ns);

    const double r_m3 = r_m * r_m * r_m;
    const double D_r_si = (PhysicsConsts::KB * temperature_K) / (8.0 * MathConsts::PI * viscosity_Pa_s * r_m3);
    D_rot_ns_ = D_r_si * 1.0e-9;
    sigma_rot_rad_ = std::sqrt(2.0 * D_rot_ns_ * dt_ns);
}

void Particle::set_box_dimensions(double box_radial_nm, double box_axial_nm) {
    box_size_ = Vec3{box_radial_nm, box_radial_nm, box_axial_nm};
}

void Particle::set_local_dipole(const Vec3& dipole) {
    abs_dipole_local_ = dipole.normalized();
}

void Particle::apply_periodic_boundaries() {
    const double hx = 0.5 * box_size_.x;
    if (pos_.x < -hx) pos_.x += box_size_.x;
    else if (pos_.x >= hx) pos_.x -= box_size_.x;

    const double hy = 0.5 * box_size_.y;
    if (pos_.y < -hy) pos_.y += box_size_.y;
    else if (pos_.y >= hy) pos_.y -= box_size_.y;

    const double hz = 0.5 * box_size_.z;
    if (pos_.z < -hz) pos_.z += box_size_.z;
    else if (pos_.z >= hz) pos_.z -= box_size_.z;
}

void Particle::move(FastRNG& rng, const Vec3& drift_nm) {
    pos_.x += sigma_trans_nm_ * rng.normal() + drift_nm.x;
    pos_.y += sigma_trans_nm_ * rng.normal() + drift_nm.y;
    pos_.z += sigma_trans_nm_ * rng.normal() + drift_nm.z;

    apply_periodic_boundaries();

    if (!is_isotropic_) {
        const Vec3 dtheta{
            sigma_rot_rad_ * rng.normal(),
            sigma_rot_rad_ * rng.normal(),
            sigma_rot_rad_ * rng.normal()
        };
        const double angle_rad = dtheta.norm();
        if (angle_rad > 1.0e-12) {
            const Quaternion dq = Quaternion::from_axis_angle(dtheta / angle_rad, angle_rad);
            orientation_ = (dq * orientation_).normalized();
        }
    }
}

bool Particle::detect_photon(FastRNG& rng, double detection_probability, double threshold) {
    if (detection_probability > threshold) {
        return rng.poisson(detection_probability) > 0;
    }
    return false;
}

uint8_t Particle::sample_polarization_channel(FastRNG& rng) const {
    if (is_isotropic_) {
        return (rng.uniform() < 0.5) ? 0 : 1;
    }
    const Vec3 mu = get_lab_dipole();
    const double i_par = mu.x * mu.x;
    const double i_perp = mu.y * mu.y;
    const double sum = i_par + i_perp;
    if (sum <= 1.0e-12) return 0;
    return (rng.uniform() < (i_par / sum)) ? 0 : 1;
}

uint16_t Particle::sample_microtime(FastRNG& rng, double time_step_ns, double irf_sigma_ns) const {
    // Déclin exponentiel : t = -tau * ln(U)
    double t_ns = -lifetime_ns_ * std::log(rng.uniform() + 1.0e-18);
    // Convolution IRF gaussienne
    if (irf_sigma_ns > 0.0) {
        t_ns += rng.normal(0.0, irf_sigma_ns);
    }
    // Repliement sur la fenêtre de répétition [0, time_step_ns[
    while (t_ns < 0.0) t_ns += time_step_ns;
    while (t_ns >= time_step_ns) t_ns -= time_step_ns;

    // Discrétisation 12 bits (canaux 0 à 4095)
    uint32_t channel_idx = static_cast<uint32_t>((t_ns / time_step_ns) * 4096.0);
    return static_cast<uint16_t>(std::min<uint32_t>(4095, channel_idx));
}