#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <cstddef>

/**
 * \brief Contiguous 1D buffer representing a regular 2D grid with bilinear interpolation.
 */
class Grid2D {
private:
    std::vector<double> data_;
    size_t nr_{0};
    size_t nz_{0};
    double r_min_{0.0};
    double r_max_{0.0};
    double z_min_{0.0};
    double z_max_{0.0};
    double dr_{0.0};
    double dz_{0.0};
    double inv_dr_{0.0};
    double inv_dz_{0.0};

public:
    Grid2D() = default;

    Grid2D(size_t nr, size_t nz, double r_min, double r_max, double z_min, double z_max) {
        init(nr, nz, r_min, r_max, z_min, z_max);
    }

    void init(size_t nr, size_t nz, double r_min, double r_max, double z_min, double z_max) {
        nr_ = nr;
        nz_ = nz;
        r_min_ = r_min;
        r_max_ = r_max;
        z_min_ = z_min;
        z_max_ = z_max;
        dr_ = (nr_ > 1) ? (r_max_ - r_min_) / static_cast<double>(nr_ - 1) : 0.0;
        dz_ = (nz_ > 1) ? (z_max_ - z_min_) / static_cast<double>(nz_ - 1) : 0.0;
        inv_dr_ = (dr_ > 0.0) ? (1.0 / dr_) : 0.0;
        inv_dz_ = (dz_ > 0.0) ? (1.0 / dz_) : 0.0;
        data_.assign(nr_ * nz_, 0.0);
    }

    inline size_t nr() const { return nr_; }
    inline size_t nz() const { return nz_; }
    inline double dr() const { return dr_; }
    inline double dz() const { return dz_; }
    inline double r_min() const { return r_min_; }
    inline double r_max() const { return r_max_; }
    inline double z_min() const { return z_min_; }
    inline double z_max() const { return z_max_; }

    inline double& at(size_t ir, size_t iz) {
        return data_[ir * nz_ + iz];
    }

    inline double at(size_t ir, size_t iz) const {
        return data_[ir * nz_ + iz];
    }

    /**
     * \brief Continuous bilinear interpolation on the 2D grid. Returns 0.0 outside boundaries.
     */
    inline double interpolate(double r, double z) const {
        if (r < r_min_ || r > r_max_ || z < z_min_ || z > z_max_) {
            return 0.0;
        }

        const double u = (r - r_min_) * inv_dr_;
        const double v = (z - z_min_) * inv_dz_;

        const size_t ir = std::min(static_cast<size_t>(u), nr_ - 2);
        const size_t iz = std::min(static_cast<size_t>(v), nz_ - 2);

        const double du = u - static_cast<double>(ir);
        const double dv = v - static_cast<double>(iz);

        const double c00 = at(ir, iz);
        const double c01 = at(ir, iz + 1);
        const double c10 = at(ir + 1, iz);
        const double c11 = at(ir + 1, iz + 1);

        return (1.0 - du) * ((1.0 - dv) * c00 + dv * c01) +
               du * ((1.0 - dv) * c10 + dv * c11);
    }

    void export_txt(const std::string& filename) const {
        std::ofstream out(filename);
        for (size_t ir = 0; ir < nr_; ++ir) {
            for (size_t iz = 0; iz < nz_; ++iz) {
                out << at(ir, iz) << ((iz + 1 == nz_) ? "" : " ");
            }
            out << "\n";
        }
    }
};