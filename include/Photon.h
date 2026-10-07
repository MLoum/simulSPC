#pragma once

#include <cstdint>

/**
 * \brief Compact binary photon record compatible with TCSPC hardware file structures.
 */
struct Photon {
    int64_t macrotime{0};    // Global clock tick index (macro-time)
    uint16_t microtime{0};   // Sub-tick arrival time / TCSPC ADC channel (0 to 65535)
    uint8_t channel{0};      // Detector channel ID (0: parallel, 1: perpendicular, etc.)

    constexpr Photon() = default;
    constexpr Photon(int64_t macro, uint16_t micro, uint8_t ch = 0)
        : macrotime(macro), microtime(micro), channel(ch) {}

    /**
     * \brief Packs record into a 64-bit word:
     *        Bits 0-39: macrotime (up to 15 hours at 50 ns resolution)
     *        Bits 40-55: microtime (16 bits)
     *        Bits 56-63: channel ID (8 bits)
     */
    inline uint64_t pack() const {
        return (static_cast<uint64_t>(macrotime) & 0xFFFFFFFFFFULL) |
               (static_cast<uint64_t>(microtime) << 40) |
               (static_cast<uint64_t>(channel) << 56);
    }

    static inline Photon unpack(uint64_t raw) {
        Photon p;
        p.macrotime = static_cast<int64_t>(raw & 0xFFFFFFFFFFULL);
        p.microtime = static_cast<uint16_t>((raw >> 40) & 0xFFFFULL);
        p.channel = static_cast<uint8_t>((raw >> 56) & 0xFFULL);
        return p;
    }
};