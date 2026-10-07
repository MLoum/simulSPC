#include <iostream>
#include <cassert>
#include <chrono>
#include <vector>
#include <cmath>
#include "Particle.h"

void test_particle_translation() {
    std::cout << "[Test 1/4] Validation hydrodynamique translationnelle (MSD)... ";

    FastRNG rng(1234);
    const double r_h = 5.0;            // 5 nm
    const double T = 293.15;           // 20 deg C
    const double eta = 1.0e-3;         // 1 mPa.s
    const double dt_ns = 50.0;         // 50 ns
    const int N_PARTICLES = 5000;
    const int N_STEPS = 100;
    const double total_time_ns = N_STEPS * dt_ns;

    Particle p_base;
    p_base.init_physics(r_h, T, eta, dt_ns);
    p_base.box_size_ = Vec3{1e9, 1e9, 1e9}; // Large box to avoid boundary wrap during MSD

    double sum_r2 = 0.0;
    for (int i = 0; i < N_PARTICLES; ++i) {
        Particle p = p_base; // Copy state
        Vec3 initial_pos = p.pos_;
        for (int s = 0; s < N_STEPS; ++s) {
            p.move(rng);
        }
        sum_r2 += (p.pos_ - initial_pos).norm_sq();
    }

    // Ensemble-averaged MSD: <|r(t) - r(0)|^2>
    double measured_msd = sum_r2 / N_PARTICLES;
    double theoretical_msd = 6.0 * p_base.D_trans_nm2_ns_ * total_time_ns;
    double relative_error = std::abs(measured_msd - theoretical_msd) / theoretical_msd;

    // Relative standard error of chi2(3) mean over 5000 samples is ~1.15%
    assert(relative_error < 0.04);

    std::cout << "OK (MSD error: " << relative_error * 100.0 << " %)" << std::endl;
}

void test_particle_rotation() {
    std::cout << "[Test 2/4] Validation diffusion rotationnelle SO(3)... ";

    FastRNG rng(5678);
    const double r_h = 2.0;            // 2 nm
    const double T = 293.15;
    const double eta = 1.0e-3;
    const double dt_ns = 1.0;          // 1 ns step
    const int N_TRAJECTORIES = 20000;

    Particle p;
    p.init_physics(r_h, T, eta, dt_ns);
    p.set_local_dipole(Vec3{0.0, 0.0, 1.0});

    double sum_cos_theta = 0.0;

    // Evaluate single-step ensemble average: <mu(0) . mu(dt)>
    for (int i = 0; i < N_TRAJECTORIES; ++i) {
        Particle p_trial = p; // Value copy
        Vec3 mu_0 = p_trial.get_lab_dipole();
        p_trial.move(rng);
        Vec3 mu_t = p_trial.get_lab_dipole();
        sum_cos_theta += mu_0.dot(mu_t);
    }

    double measured_autocorr = sum_cos_theta / N_TRAJECTORIES;
    // Theoretical 1st Legendre correlation: exp(-2 * D_rot * dt)
    double theoretical_autocorr = std::exp(-2.0 * p.D_rot_ns_ * dt_ns);
    double error = std::abs(measured_autocorr - theoretical_autocorr);

    assert(error < 0.01);
    std::cout << "OK (Correlation mesuree = " << measured_autocorr
              << ", theorie = " << theoretical_autocorr << ")" << std::endl;
}

void test_rule_of_zero() {
    std::cout << "[Test 3/4] Semantique de copie et absence de fuite memoire... ";
    std::vector<Particle> particles;
    particles.reserve(1000);

    for (int i = 0; i < 1000; ++i) {
        Particle p;
        p.pos_ = Vec3{static_cast<double>(i), 0.0, 0.0};
        particles.push_back(p); // Test copy constructor
    }

    // Test copy assignment
    Particle p_copy = particles[500];
    assert(p_copy.pos_.x == 500.0);

    // Test vector resize/reallocation without memory leaks
    particles.clear();
    particles.shrink_to_fit();

    std::cout << "OK" << std::endl;
}

void benchmark_particle_throughput() {
    std::cout << "\n--- Benchmark Moteur Particule (10 000 particules x 1 000 ticks) ---" << std::endl;
    FastRNG rng(42);
    const int N_PARTICLES = 10000;
    const int N_TICKS = 1000;

    std::vector<Particle> particles(N_PARTICLES);
    for (auto& p : particles) {
        p.init_physics(10.0, 293.15, 1.0e-3, 50.0);
        p.set_box_dimensions(1000.0, 2000.0);
    }

    auto t0 = std::chrono::steady_clock::now();
    for (int t = 0; t < N_TICKS; ++t) {
        for (int i = 0; i < N_PARTICLES; ++i) {
            particles[i].move(rng);
        }
    }
    auto t1 = std::chrono::steady_clock::now();

    std::chrono::duration<double, std::milli> ms = t1 - t0;
    double total_evals = static_cast<double>(N_PARTICLES) * N_TICKS;
    double throughput = (total_evals / (ms.count() * 1.0e-3)) / 1.0e6;

    std::cout << "Duree totale : " << ms.count() << " ms" << std::endl;
    std::cout << "Debit : " << throughput << " millions de pas particule/sec (mono-coeur)" << std::endl;
}

int main() {
    test_particle_translation();
    test_particle_rotation();
    test_rule_of_zero();
    benchmark_particle_throughput();

    std::cout << "\nTous les tests physiques de l'Etape 2 ont reussi." << std::endl;
    return 0;
}