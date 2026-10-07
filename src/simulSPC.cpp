#include <iostream>
#include <chrono>
#include <string>
#include "Experiment.h"

int main(int argc, char* argv[]) {
    std::string config_file = "parameters.cfg";
    std::string output_file = "photon_list.ptn";

    if (argc >= 2) {
        config_file = argv[1];
    }
    if (argc >= 3) {
        output_file = argv[2];
    }

    std::cout << "========================================================\n";
    std::cout << "           SimulSPC - Single Particle Spectroscopy       \n";
    std::cout << "========================================================\n";
    std::cout << "Fichier de configuration : " << config_file << "\n";
    std::cout << "Fichier de sortie        : " << output_file << "\n";

    Experiment exp;
    auto t0 = std::chrono::steady_clock::now();

    if (!exp.init_from_config(config_file)) {
        std::cerr << "Erreur critique : initialisation impossible.\n";
        return 1;
    }

    auto t1 = std::chrono::steady_clock::now();
    std::cout << "Particules instanciees   : " << exp.particles_.size() << "\n";
    std::cout << "Threads OpenMP actifs    : " << exp.num_threads_ << "\n";
    std::cout << "Duree simulee            : " << exp.simulation_time_s_ * 1000.0 << " ms ("
              << exp.total_ticks_ << " ticks de " << exp.time_step_ns_ << " ns)\n";
    std::cout << "Initialisation terminee en "
              << std::chrono::duration<double, std::milli>(t1 - t0).count() << " ms.\n\n";

    std::cout << "Demarrage de la simulation...\n";
    auto t_run0 = std::chrono::steady_clock::now();
    exp.run(true);
    auto t_run1 = std::chrono::steady_clock::now();

    std::chrono::duration<double> duration_sec = t_run1 - t_run0;
    std::cout << "\nSimulation terminee avec succes.\n";
    std::cout << "Duree d'execution        : " << duration_sec.count() << " s\n";
    std::cout << "Photons collectes        : " << exp.photon_stream_.size() << "\n";

    if (duration_sec.count() > 0.0) {
        double rate_khz = (exp.photon_stream_.size() / duration_sec.count()) * 1e-3;
        std::cout << "Debit de generation      : " << rate_khz << " kPhotons/s (temps machine)\n";
    }

    exp.export_photon_stream(output_file);
    std::cout << "Fichier binaire exporte  : " << output_file << "\n";
    std::cout << "========================================================\n";

    return 0;
}