#include <iostream>
#include <cassert>
#include <chrono>
#include <vector>
#include <numeric>
#include <cmath>
#include "Math3D.h"
#include "FastRNG.h"

namespace {
    constexpr double PI = 3.14159265358979323846;
}

void test_vec3() {
    std::cout << "[Test 1/5] Geometrie vectorielle Vec3... ";
    Vec3 a{1.0, 2.0, 3.0};
    Vec3 b{4.0, -5.0, 6.0};

    assert(std::abs(a.dot(b) - 12.0) < 1e-12);

    Vec3 c = a.cross(b);
    assert(std::abs(c.x - 27.0) < 1e-12);
    assert(std::abs(c.y - 6.0) < 1e-12);
    assert(std::abs(c.z - (-13.0)) < 1e-12);

    Vec3 d{0.0, 3.0, 4.0};
    assert(std::abs(d.norm() - 5.0) < 1e-12);
    assert(std::abs(d.normalized().norm() - 1.0) < 1e-12);

    std::cout << "OK" << std::endl;
}

void test_quaternion() {
    std::cout << "[Test 2/5] Rotations Quaternion... ";
    Vec3 v{1.0, 0.0, 0.0};

    // Rotation de 90 degres autour de Z -> doit donner (0, 1, 0)
    Quaternion qz = Quaternion::from_axis_angle(Vec3{0.0, 0.0, 1.0}, PI * 0.5);
    Vec3 rotated = qz.rotate(v);
    assert(std::abs(rotated.x - 0.0) < 1e-12);
    assert(std::abs(rotated.y - 1.0) < 1e-12);
    assert(std::abs(rotated.z - 0.0) < 1e-12);

    // Stabilite apres 100 000 micro-rotations aleatoires
    Quaternion q_accum;
    FastRNG rng(42);
    for (int i = 0; i < 100000; ++i) {
        Vec3 axis{rng.normal(), rng.normal(), rng.normal()};
        Quaternion dq = Quaternion::from_axis_angle(axis, 1e-4);
        q_accum = dq * q_accum;
    }
    // Verifie la derive de norme avant renormalisation
    double drift = std::abs(q_accum.norm() - 1.0);
    assert(drift < 1e-7);
    q_accum.normalize();
    assert(std::abs(q_accum.norm() - 1.0) < 1e-14);

    std::cout << "OK (derive apres 100k pas : " << drift << ")" << std::endl;
}

void test_rng_uniform() {
    std::cout << "[Test 3/5] Distribution uniforme... ";
    FastRNG rng(12345);
    const size_t N = 5000000;
    double sum = 0.0;
    double sum_sq = 0.0;

    for (size_t i = 0; i < N; ++i) {
        double u = rng.uniform();
        sum += u;
        sum_sq += u * u;
    }

    double mean = sum / N;
    double variance = (sum_sq / N) - (mean * mean);

    // Moyenne attendue: 0.5, Variance attendue: 1/12 ~ 0.083333
    assert(std::abs(mean - 0.5) < 1e-3);
    assert(std::abs(variance - (1.0 / 12.0)) < 1e-3);

    std::cout << "OK (moyenne = " << mean << ", var = " << variance << ")" << std::endl;
}

void test_rng_gaussian() {
    std::cout << "[Test 4/5] Distribution gaussienne Ziggurat... ";
    FastRNG rng(67890);
    const size_t N = 5000000;
    double sum = 0.0;
    double sum_sq = 0.0;

    for (size_t i = 0; i < N; ++i) {
        double g = rng.normal();
        sum += g;
        sum_sq += g * g;
    }

    double mean = sum / N;
    double variance = (sum_sq / N) - (mean * mean);

    // Moyenne attendue: 0.0, Variance attendue: 1.0
    assert(std::abs(mean) < 2e-3);
    assert(std::abs(variance - 1.0) < 2e-3);

    std::cout << "OK (moyenne = " << mean << ", var = " << variance << ")" << std::endl;
}

void test_rng_poisson() {
    std::cout << "[Test 5/5] Distribution de Poisson... ";
    FastRNG rng(111213);
    const size_t N = 2000000;
    const double lambda = 0.045; // Flux photonique typique sur 50 ns

    double sum = 0.0;
    double sum_sq = 0.0;

    for (size_t i = 0; i < N; ++i) {
        double k = static_cast<double>(rng.poisson(lambda));
        sum += k;
        sum_sq += k * k;
    }

    double mean = sum / N;
    double variance = (sum_sq / N) - (mean * mean);

    // Moyenne et variance attendues = lambda
    assert(std::abs(mean - lambda) < 1e-3);
    assert(std::abs(variance - lambda) < 1e-3);

    std::cout << "OK (moyenne = " << mean << ", var = " << variance << ")" << std::endl;
}

void benchmark_speed() {
    std::cout << "\n--- Benchmark de performance (10^7 tirages) ---" << std::endl;
    FastRNG rng(999);
    const size_t N = 10000000;

    auto t0 = std::chrono::steady_clock::now();
    double dummy = 0.0;
    for (size_t i = 0; i < N; ++i) {
        dummy += rng.normal();
    }
    auto t1 = std::chrono::steady_clock::now();

    std::chrono::duration<double, std::milli> ms = t1 - t0;
    double throughput = (static_cast<double>(N) / (ms.count() * 1e-3)) / 1e6;
    std::cout << "Temps Ziggurat N(0,1) : " << ms.count() << " ms ("
              << throughput << " millions d'echantillons/sec, dummy = " << dummy << ")" << std::endl;
}

int main() {
    test_vec3();
    test_quaternion();
    test_rng_uniform();
    test_rng_gaussian();
    test_rng_poisson();
    benchmark_speed();

    std::cout << "\nTous les tests unitaires ont reussi sans aucune erreur." << std::endl;
    return 0;
}