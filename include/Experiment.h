#pragma once

#include <vector>
#include <string>
#include <cstdint>
#include "Particle.h"
#include "Solvent.h"
#include "OpticalSetup.h"
#include "MDF_Gaussian.h"
#include "Photon.h"
#include "FastRNG.h"
#include "Config.h"

class Experiment {
public:
    double time_step_ns_{50.0};
    double simulation_time_s_{0.01};
    uint64_t total_ticks_{0};
    int64_t macro_clock_{0};

    Solvent solvent_;
    OpticalSetup optical_setup_;
    MDF_Gaussian mdf_;

    std::vector<Particle> particles_;
    double box_radial_nm_{1000.0};
    double box_axial_nm_{2000.0};

    double deadtime_ns_{35.0};
    double irf_sigma_ns_{0.106}; // FWHM = 250 ps -> sigma ~ 106 ps
    double last_arrival_ch0_ns_{-1.0e9};
    double last_arrival_ch1_ns_{-1.0e9};
    std::vector<Photon> photon_stream_;

    int num_threads_{1};
    std::vector<FastRNG> thread_rngs_;

public:
    Experiment() = default;

    bool init_from_config(const std::string& config_file);
    void tick();
    void run(bool show_progress = true);
    void export_photon_stream(const std::string& binary_filepath) const;
};