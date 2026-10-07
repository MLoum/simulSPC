#pragma once

#include <cstdint>
#include <cmath>
#include <mutex>
#include <random>

/**
 * \brief High-performance PRNG combining xoshiro256** and the 256-layer Ziggurat algorithm.
 * \details Designed for thread-local instances to avoid lock contention and memory races.
 */
class FastRNG {
private:
    static constexpr double PI = 3.14159265358979323846;
    static constexpr double SQRT1_2 = 0.70710678118654752440;

    uint64_t s_[4];

    static inline uint64_t rotl(const uint64_t x, int k) {
        return (x << k) | (x >> (64 - k));
    }

    static inline uint64_t splitmix64(uint64_t& x) {
        uint64_t z = (x += 0x9E3779B97F4A7C15ULL);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        return z ^ (z >> 31);
    }

    // --- Ziggurat Lookup Tables (256 layers) ---
    static constexpr int ZIG_N = 256;
    static inline double zig_x[ZIG_N];
    static inline double zig_y[ZIG_N];
    static inline double zig_r_y0_ratio;
    static inline double zig_r;
    static inline std::once_flag zig_init_flag;

    /**
     * \brief Precomputes Ziggurat tables using binary search for layer convergence up to y = 1.0.
     */
    static void init_ziggurat_tables() {
        double r_low = 3.0, r_high = 4.0;
        double r = 3.65, v = 0.0;

        for (int iter = 0; iter < 80; ++iter) {
            r = 0.5 * (r_low + r_high);
            double y0 = std::exp(-0.5 * r * r);
            double tail = std::sqrt(PI * 0.5) * std::erfc(r * SQRT1_2);
            v = r * y0 + tail;

            double y = y0;
            double x = r;
            bool overshoot = false;
            for (int i = 1; i < ZIG_N; ++i) {
                y += v / x;
                if (y >= 1.0) {
                    overshoot = true;
                    break;
                }
                x = std::sqrt(-2.0 * std::log(y));
            }
            if (overshoot) {
                r_low = r; // Volume too large -> r too small
            } else {
                r_high = r; // Volume too small -> r too large
            }
        }

        zig_r = r;
        double y0 = std::exp(-0.5 * r * r);
        double tail = std::sqrt(PI * 0.5) * std::erfc(r * SQRT1_2);
        v = r * y0 + tail;

        zig_x[0] = r;
        zig_y[0] = y0;
        zig_r_y0_ratio = (r * y0) / v;

        for (int i = 1; i < ZIG_N - 1; ++i) {
            zig_y[i] = zig_y[i - 1] + v / zig_x[i - 1];
            zig_x[i] = std::sqrt(-2.0 * std::log(zig_y[i]));
        }

        // Top layer boundary
        zig_y[ZIG_N - 1] = 1.0;
        zig_x[ZIG_N - 1] = 0.0;
    }

public:
    explicit FastRNG(uint64_t seed = 0) {
        std::call_once(zig_init_flag, init_ziggurat_tables);
        set_seed(seed);
    }

    /**
     * \brief Seeds the internal state using SplitMix64.
     */
    void set_seed(uint64_t seed) {
        if (seed == 0) {
            std::random_device rd;
            seed = (static_cast<uint64_t>(rd()) << 32) | rd();
        }
        uint64_t sm_state = seed;
        s_[0] = splitmix64(sm_state);
        s_[1] = splitmix64(sm_state);
        s_[2] = splitmix64(sm_state);
        s_[3] = splitmix64(sm_state);
    }

    /**
     * \brief Generates uniform 64-bit unsigned integers via xoshiro256**.
     */
    inline uint64_t next_u64() {
        const uint64_t result = rotl(s_[1] * 5, 7) * 9;
        const uint64_t t = s_[1] << 17;

        s_[2] ^= s_[0];
        s_[3] ^= s_[1];
        s_[1] ^= s_[2];
        s_[0] ^= s_[3];

        s_[2] ^= t;
        s_[3] = rotl(s_[3], 45);

        return result;
    }

    /**
     * \brief Returns a uniform double in [0, 1) with 53 bits of precision.
     */
    inline double uniform() {
        return (next_u64() >> 11) * 0x1.0p-53;
    }

    inline double uniform(double min, double max) {
        return min + (max - min) * uniform();
    }

    /**
     * \brief Samples standard normal distribution N(0, 1) using Ziggurat method.
     */
    double normal() {
        while (true) {
            uint64_t u = next_u64();
            uint8_t i = u & 0xFF;                  // Layer index (0 to 255)
            bool sign = (u & 0x100) != 0;          // Bit 8 encodes sign
            double u_val = ((u >> 9) | 1) * 0x1.0p-55; // 55 random bits in (0, 1)

            // Layer 0: base rectangle or exponential tail
            if (i == 0) {
                if (u_val < zig_r_y0_ratio) {
                    double x = (u_val / zig_r_y0_ratio) * zig_r;
                    return sign ? -x : x;
                }
                // Marsaglia tail sampling
                double x_tail, y_tail;
                do {
                    x_tail = -std::log(uniform() + 1e-18) / zig_r;
                    y_tail = -std::log(uniform() + 1e-18);
                } while (2.0 * y_tail < x_tail * x_tail);
                double x = zig_r + x_tail;
                return sign ? -x : x;
            }

            // Layers 1 to 255: fast inner-box rejection test (>98% accepted here)
            double x = u_val * zig_x[i - 1];
            if (x < zig_x[i]) {
                return sign ? -x : x;
            }

            // Precise wedge boundary evaluation
            double y = zig_y[i - 1] + uniform() * (zig_y[i] - zig_y[i - 1]);
            if (y < std::exp(-0.5 * x * x)) {
                return sign ? -x : x;
            }
        }
    }

    inline double normal(double mean, double stddev) {
        return mean + stddev * normal();
    }

    /**
     * \brief Poisson distribution generator.
     * \details Uses Knuth's exact multiplication method for small lambda (< 30)
     *          and Gaussian approximation for higher photon fluxes.
     */
    uint64_t poisson(double lambda) {
        if (lambda <= 0.0) return 0;
        if (lambda < 30.0) {
            double L = std::exp(-lambda);
            double k = 0;
            double p = 1.0;
            do {
                k += 1.0;
                p *= (uniform() + 1e-18);
            } while (p > L);
            return static_cast<uint64_t>(k - 1.0);
        }
        double val = std::round(lambda + std::sqrt(lambda) * normal());
        return (val > 0.0) ? static_cast<uint64_t>(val) : 0;
    }
};