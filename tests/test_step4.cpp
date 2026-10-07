#include <iostream>
#include <cassert>
#include <chrono>
#include "Experiment.h"

void test_photon_packing() {
    std::cout << "[Test 1/3] Encodage/decodage binaire 64-bit Photon... ";

    Photon p1(123456789LL, 42000, 2);
    uint64_t raw = p1.pack();

    Photon p2 = Photon::unpack(raw);
    assert(p2.macrotime == p1.macrotime);
    assert(p2.microtime == p1.microtime);
    assert(p2.channel == p1.channel);

    std::cout << "OK" << std::endl;
}

void test_config_parser() {
    std::cout << "[Test 2/3] Robustesse Config parser... ";

    std::ofstream tmp("test_param.cfg");
    tmp << "# Commentaire test\n"
        << "concentration_nM = 15.5\n"
        << "r_hydro_nm = 7.2\n"
        << "is_isotropic = true\n"
        << "laser_power_mW = 0.02\n";
    tmp.close();

    Config cfg;
    assert(cfg.load("test_param.cfg"));
    assert(std::abs(cfg.get_double("concentration_nM", 0.0) - 15.5) < 1e-6);
    assert(std::abs(cfg.get_double("r_hydro_nm", 0.0) - 7.2) < 1e-6);
    assert(cfg.get_bool("is_isotropic", false) == true);

    std::cout << "OK" << std::endl;
}

void test_full_simulation_run() {
    std::cout << "[Test 3/3] Simulation complete 2 ms (40 000 ticks) avec OpenMP...\n";

    Experiment exp;
    // 10 nM, r = 5 nm, Rhodamine 6G cross-section, 20 uW laser, dt = 50 ns, sim_time = 2 ms
    exp.init(10.0, 5.0, 3.0e-16, 0.8, false,
             0.02, 50.0, 0.002, 1000.0, 2000.0, 0.0, 0);

    std::cout << "Nombre de particules en simulation : " << exp.particles_.size() << std::endl;
    std::cout << "Nombre de threads actifs : " << exp.num_threads_ << std::endl;

    auto t0 = std::chrono::steady_clock::now();
    exp.run(true);
    auto t1 = std::chrono::steady_clock::now();

    std::chrono::duration<double, std::milli> ms = t1 - t0;
    std::cout << "Duree simulation : " << ms.count() << " ms" << std::endl;
    std::cout << "Nombre total de photons enregistres : " << exp.photon_stream_.size() << std::endl;

    assert(!exp.photon_stream_.empty());

    // Export validation
    exp.export_photon_stream("test_stream.ptn");
    std::ifstream in("test_stream.ptn", std::ios::binary | std::ios::ate);
    assert(in.is_open());
    size_t file_size = in.tellg();
    assert(file_size == exp.photon_stream_.size() * sizeof(uint64_t));

    std::cout << "Fichier binaire photon_stream valide (" << file_size << " octets).\n";
}

int main() {
    test_photon_packing();
    test_config_parser();
    test_full_simulation_run();

    std::cout << "\nTous les tests de l'Etape 4 ont reussi avec succes." << std::endl;
    return 0;
}