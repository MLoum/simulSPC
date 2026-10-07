#pragma once

#include "Particle.h"
#include "OpticalSetup.h"
#include "Grid2D.h"

/**
 * \brief Abstract interface for the Molecular Detection Efficiency Function (MDF).
 */
class MDF {
public:
    virtual ~MDF() = default;

    /**
     * \brief Calculates the expected number of detected photons (Poisson lambda) over dt_ns.
     */
    virtual double get_detection_probability(const Particle& p, double dt_ns) const = 0;

    virtual const Grid2D& get_grid() const = 0;
};