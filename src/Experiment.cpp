#include "Experiment.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <omp.h>

bool Experiment::init_from_config(const std::string& config_file) {
    Config cfg;
    if (!cfg.load(config_file)) {
        std::cerr << "Erreur: Impossible de lire le fichier de configuration: " << config_file << std::endl;
        return false;
    }

    time_step_ns_ = cfg.get_double("time_step_ns", 50.0);
    simulation_time_s_ = cfg.get_double("simulation_time_s", 0.01);
    total_ticks_ = static_cast<uint64_t>(std::round(simulation_time_s_ / (time_step_ns_ * 1.0e-9)));
    macro_clock_ = 0;

    box_radial_nm_ = cfg.get_double("box_size_radial_nm", 1000.0);
    box_axial_nm_ = cfg.get_double("box_size_axial_nm", 2000.0);
    deadtime_ns_ = cfg.get_double("deadtime_ns", 35.0);

    double irf_fwhm_ps = cfg.get_double("irf_fwhm_ps", 250.0);
    irf_sigma_ns_ = (irf_fwhm_ps * 1.0e-3) / 2.35482; // FWHM to standard deviation

    last_arrival_ch0_ns_ = -1.0e9;
    last_arrival_ch1_ns_ = -1.0e9;
    photon_stream_.clear();

    int max_threads = omp_get_max_threads();
    int req_threads = cfg.get_int("num_threads", 0);
    num_threads_ = (req_threads > 0) ? std::min(req_threads, max_threads) : max_threads;
    thread_rngs_.resize(num_threads_);
    for (int i = 0; i < num_threads_; ++i) {
        thread_rngs_[i].set_seed(10007ULL * (i + 1) + 42ULL);
    }

    // 1. Solvant
    double flow = cfg.get_double("mean_flow_uL_min", 0.0);
    solvent_.init(1.0e-3, 293.15, 1.333, box_radial_nm_, box_axial_nm_, flow, 20.0);

    // 2. Optique
    optical_setup_.laser.wavelength_nm = 488.0;
    optical_setup_.laser.power_mW = cfg.get_double("laser_power_mW", 0.015);
    optical_setup_.laser.polarization = Vec3{1.0, 0.0, 0.0};
    optical_setup_.objective.NA = 1.2;
    optical_setup_.objective.immersion_index = 1.333;
    optical_setup_.detection.optics_transmission = 0.70;
    optical_setup_.detection.pinhole_efficiency = 0.70;
    optical_setup_.detection.detector_efficiency = 0.65;

    // 3. MDF
    mdf_.init(optical_setup_, 4.0, box_radial_nm_, box_axial_nm_, 20.0);

    // 4. Instanciation des espèces
    particles_.clear();
    const double volume_liters = (box_radial_nm_ * box_radial_nm_ * box_axial_nm_) * 1.0e-24;
    FastRNG& rng = thread_rngs_[0];
    const double hr = 0.5 * box_radial_nm_;
    const double hz = 0.5 * box_axial_nm_;

    for (const auto& sp : cfg.species()) {
        const double total_moles = (sp.concentration_nM * 1.0e-9) * volume_liters;
        size_t n_part = static_cast<size_t>(std::round(total_moles * PhysicsConsts::NA));

        for (size_t i = 0; i < n_part; ++i) {
            Particle p;
            double r_val = sp.r_hydro_nm;
            if (sp.r_hydro_dispersion_pct > 0.0) {
                double sigma_r = sp.r_hydro_nm * (sp.r_hydro_dispersion_pct / 100.0);
                r_val = std::max(0.5, rng.normal(sp.r_hydro_nm, sigma_r));
            }

            double tau_val = sp.lifetime_ns;
            if (sp.lifetime_dispersion_pct > 0.0) {
                double sigma_tau = sp.lifetime_ns * (sp.lifetime_dispersion_pct / 100.0);
                tau_val = std::max(0.1, rng.normal(sp.lifetime_ns, sigma_tau));
            }

            p.init_physics(r_val, solvent_.temperature_K, solvent_.viscosity_Pa_s,
                           time_step_ns_, tau_val, sp.cross_section_cm2, sp.quantum_yield, sp.is_isotropic);
            p.set_box_dimensions(box_radial_nm_, box_axial_nm_);
            p.pos_ = Vec3{rng.uniform(-hr, hr), rng.uniform(-hr, hr), rng.uniform(-hz, hz)};
            p.set_local_dipole(Vec3{rng.normal(), rng.normal(), rng.normal()});

            particles_.push_back(p);
        }
    }

    if (particles_.empty()) {
        std::cerr << "Attention: Aucune particule instanciee. Verifiez le volume ou la concentration." << std::endl;
    }
    return true;
}

void Experiment::tick() {
    struct CandidatePhoton {
        double sub_tick_time_ns;
        uint16_t microtime;
        uint8_t channel;
    };

    std::vector<std::vector<CandidatePhoton>> thread_candidates(num_threads_);

    #pragma omp parallel num_threads(num_threads_)
    {
        int tid = omp_get_thread_num();
        FastRNG& rng = thread_rngs_[tid];

        #pragma omp for schedule(static)
        for (int i = 0; i < static_cast<int>(particles_.size()); ++i) {
            Particle& p = particles_[i];

            Vec3 drift = solvent_.get_drift_nm(p.pos_, time_step_ns_);
            p.move(rng, drift);

            double lambda = mdf_.get_detection_probability(p, time_step_ns_);
            if (p.detect_photon(rng, lambda, 1.0e-7)) {
                double sub_tick_ns = rng.uniform(0.0, time_step_ns_);
                uint8_t ch = p.sample_polarization_channel(rng);
                uint16_t micro = p.sample_microtime(rng, time_step_ns_, irf_sigma_ns_);
                thread_candidates[tid].push_back({sub_tick_ns, micro, ch});
            }
        }
    }

    std::vector<CandidatePhoton> all_candidates;
    for (const auto& tc : thread_candidates) {
        all_candidates.insert(all_candidates.end(), tc.begin(), tc.end());
    }

    if (all_candidates.empty()) {
        macro_clock_++;
        return;
    }

    std::sort(all_candidates.begin(), all_candidates.end(),
              [](const CandidatePhoton& a, const CandidatePhoton& b) {
                  return a.sub_tick_time_ns < b.sub_tick_time_ns;
              });

    const double tick_start_ns = static_cast<double>(macro_clock_) * time_step_ns_;

    for (const auto& cand : all_candidates) {
        double arrival_abs_ns = tick_start_ns + cand.sub_tick_time_ns;
        double& last_arrival = (cand.channel == 0) ? last_arrival_ch0_ns_ : last_arrival_ch1_ns_;

        if (arrival_abs_ns - last_arrival >= deadtime_ns_) {
            photon_stream_.emplace_back(macro_clock_, cand.microtime, cand.channel);
            last_arrival = arrival_abs_ns;
        }
    }

    macro_clock_++;
}

void Experiment::run(bool show_progress) {
    const uint64_t step_notify = std::max<uint64_t>(1, total_ticks_ / 10);

    for (uint64_t t = 0; t < total_ticks_; ++t) {
        tick();

        if (show_progress && (t + 1) % step_notify == 0) {
            int pct = static_cast<int>(((t + 1) * 100) / total_ticks_);
            std::cout << "Progression: " << pct << "% ("
                      << photon_stream_.size() << " photons detectes)" << std::endl;
        }
    }
}

void Experiment::export_photon_stream(const std::string& binary_filepath) const {
    std::ofstream out(binary_filepath, std::ios::binary);
    if (!out.is_open()) return;

    std::vector<uint64_t> packed(photon_stream_.size());
    for (size_t i = 0; i < photon_stream_.size(); ++i) {
        packed[i] = photon_stream_[i].pack();
    }

    out.write(reinterpret_cast<const char*>(packed.data()), packed.size() * sizeof(uint64_t));
}