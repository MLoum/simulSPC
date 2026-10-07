#include <iostream>
#include <cassert>
#include <chrono>
#include <cmath>
#include "Grid2D.h"
#include "OpticalSetup.h"
#include "Solvent.h"
#include "MDF_Gaussian.h"

void test_grid_interpolation() {
    std::cout << "[Test 1/4] Validation interpolation bilineaire Grid2D... ";

    Grid2D grid(101, 101, 0.0, 100.0, -50.0, 50.0);
    // Fill with linear function f(r, z) = 2.5 * r - 1.5 * z + 3.0
    for (size_t ir = 0; ir < grid.nr(); ++ir) {
        double r = ir * grid.dr();
        for (size_t iz = 0; iz < grid.nz(); ++iz) {
            double z = -50.0 + iz * grid.dz();
            grid.at(ir, iz) = 2.5 * r - 1.5 * z + 3.0;
        }
    }

    // Test intermediate non-grid points
    double val = grid.interpolate(33.3, 12.7);
    double expected = 2.5 * 33.3 - 1.5 * 12.7 + 3.0;
    assert(std::abs(val - expected) < 1e-10);

    // Test out of bounds returns zero
    assert(grid.interpolate(150.0, 0.0) == 0.0);
    assert(grid.interpolate(50.0, -80.0) == 0.0);

    std::cout << "OK" << std::endl;
}

void test_solvent_poiseuille() {
    std::cout << "[Test 2/4] Validation ecoulement de Poiseuille rectangulaire... ";

    Solvent solvent;
    solvent.init(1.0e-3, 293.15, 1.333, 1000.0, 2000.0, 0.1, 20.0); // 0.1 uL/min

    // Channel walls must have zero velocity (no-slip condition)
    Vec3 drift_wall_y = solvent.get_drift_nm(Vec3{0.0, -500.0, 0.0}, 50.0);
    Vec3 drift_wall_z = solvent.get_drift_nm(Vec3{0.0, 0.0, -1000.0}, 50.0);
    assert(std::abs(drift_wall_y.x) < 1e-6);
    assert(std::abs(drift_wall_z.x) < 1e-6);

    // Maximum velocity at channel center (y = 0, z = 0)
    Vec3 drift_center = solvent.get_drift_nm(Vec3{0.0, 0.0, 0.0}, 50.0);
    assert(drift_center.x > 0.0);

    std::cout << "OK (Deplacement au centre sur 50 ns : " << drift_center.x << " nm)" << std::endl;
}

void test_mdf_physics_and_photon_rate() {
    std::cout << "[Test 3/4] Validation physique MDF et echelle de comptage photon... ";

    OpticalSetup opt;
    opt.laser.wavelength_nm = 488.0;
    opt.laser.power_mW = 0.01; // 10 uW: standard excitation power in single-molecule FCS
    opt.laser.polarization = Vec3{1.0, 0.0, 0.0};
    opt.objective.NA = 1.2;
    opt.objective.immersion_index = 1.333;
    opt.detection.optics_transmission = 0.70;
    opt.detection.pinhole_efficiency = 0.70;
    opt.detection.detector_efficiency = 0.65;

    MDF_Gaussian mdf;
    mdf.init(opt, 4.0, 1000.0, 2000.0, 25.0);

    // Validate lateral waist: w_xy ~ 0.66 * (1.22 * 488 / (2 * 1.2)) ~ 163.7 nm
    assert(mdf.w_xy_nm() > 150.0 && mdf.w_xy_nm() < 180.0);

    // Particle at focal center
    Particle p;
    p.pos_ = Vec3{0.0, 0.0, 0.0};
    p.abs_cross_section_cm2_ = 3.0e-16; // 3e-16 cm^2 (~ 80 000 M^-1 cm^-1)
    p.quantum_yield_ = 0.80;
    p.is_isotropic_ = true;

    const double dt_ns = 50.0;
    double lambda_per_tick = mdf.get_detection_probability(p, dt_ns);

    // At 10 uW with full collection pipeline (eta_tot ~ 9%), count rate at focus is ~ 1 MHz
    double count_rate_cps = lambda_per_tick / (dt_ns * 1.0e-9);
    assert(count_rate_cps > 1.0e5 && count_rate_cps < 3.0e6);

    std::cout << "OK (Taux au centre : " << count_rate_cps * 1e-3
              << " kHz, proba/tick de 50 ns = " << lambda_per_tick << ")" << std::endl;
}

void benchmark_mdf_throughput() {
    std::cout << "\n--- Benchmark Grille MDF (10^7 evaluations) ---" << std::endl;

    OpticalSetup opt;
    opt.laser.power_mW = 0.01;
    MDF_Gaussian mdf;
    mdf.init(opt, 4.0, 1000.0, 2000.0, 10.0);

    Particle p;
    p.is_isotropic_ = true;
    p.abs_cross_section_cm2_ = 3.0e-16;
    p.quantum_yield_ = 0.8;

    const int N_EVALS = 10000000;
    double sum = 0.0;

    auto t0 = std::chrono::steady_clock::now();
    for (int i = 0; i < N_EVALS; ++i) {
        p.pos_.x = (i % 200) - 100.0;
        p.pos_.y = ((i / 200) % 200) - 100.0;
        p.pos_.z = (i % 400) - 200.0;
        sum += mdf.get_detection_probability(p, 50.0);
    }
    auto t1 = std::chrono::steady_clock::now();

    std::chrono::duration<double, std::milli> ms = t1 - t0;
    double throughput = (static_cast<double>(N_EVALS) / (ms.count() * 1.0e-3)) / 1.0e6;

    std::cout << "Duree pour 10^7 interpolations : " << ms.count() << " ms" << std::endl;
    std::cout << "Debit : " << throughput << " millions d'evaluations MDF/sec (mono-coeur, dummy = "
              << sum << ")" << std::endl;
}

int main() {
    test_grid_interpolation();
    test_solvent_poiseuille();
    test_mdf_physics_and_photon_rate();
    benchmark_mdf_throughput();

    std::cout << "\nTous les tests de l'Etape 3 ont reussi." << std::endl;
    return 0;
}