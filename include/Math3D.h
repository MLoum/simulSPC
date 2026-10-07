#pragma once

#include <cmath>

namespace MathConsts {
    constexpr double PI = 3.14159265358979323846;
}

/**
 * \brief Lightweight 3D vector allocated on the stack (zero dynamic allocation).
 */
struct Vec3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    constexpr Vec3() = default;
    constexpr Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

    constexpr Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    constexpr Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    constexpr Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
    constexpr Vec3 operator/(double s) const { double inv = 1.0 / s; return {x * inv, y * inv, z * inv}; }

    Vec3& operator+=(const Vec3& o) { x += o.x; y += o.y; z += o.z; return *this; }
    Vec3& operator-=(const Vec3& o) { x -= o.x; y -= o.y; z -= o.z; return *this; }
    Vec3& operator*=(double s) { x *= s; y *= s; z *= s; return *this; }
    Vec3& operator/=(double s) { double inv = 1.0 / s; x *= inv; y *= inv; z *= inv; return *this; }

    constexpr double dot(const Vec3& o) const { return x * o.x + y * o.y + z * o.z; }

    constexpr Vec3 cross(const Vec3& o) const {
        return {
            y * o.z - z * o.y,
            z * o.x - x * o.z,
            x * o.y - y * o.x
        };
    }

    double norm_sq() const { return x * x + y * y + z * z; }
    double norm() const { return std::sqrt(norm_sq()); }

    Vec3 normalized() const {
        double n = norm();
        return (n > 0.0) ? (*this / n) : Vec3{0.0, 0.0, 0.0};
    }
};

inline constexpr Vec3 operator*(double s, const Vec3& v) { return v * s; }

/**
 * \brief Unit quaternion representation for SO(3) 3D rotations without gimbal lock or matrix drift.
 */
struct Quaternion {
    double w = 1.0;
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    constexpr Quaternion() = default;
    constexpr Quaternion(double w_, double x_, double y_, double z_) : w(w_), x(x_), y(y_), z(z_) {}

    static Quaternion from_axis_angle(const Vec3& axis, double angle_rad) {
        Vec3 u = axis.normalized();
        double half = 0.5 * angle_rad;
        double s = std::sin(half);
        return {std::cos(half), u.x * s, u.y * s, u.z * s};
    }

    constexpr Quaternion operator*(const Quaternion& q) const {
        return {
            w * q.w - x * q.x - y * q.y - z * q.z,
            w * q.x + x * q.w + y * q.z - z * q.y,
            w * q.y - x * q.z + y * q.w + z * q.x,
            w * q.z + x * q.y - y * q.x + z * q.w
        };
    }

    double norm_sq() const { return w * w + x * x + y * y + z * z; }
    double norm() const { return std::sqrt(norm_sq()); }

    Quaternion normalized() const {
        double n = norm();
        if (n > 0.0) {
            double inv = 1.0 / n;
            return {w * inv, x * inv, y * inv, z * inv};
        }
        return {1.0, 0.0, 0.0, 0.0};
    }

    void normalize() {
        *this = normalized();
    }

    /**
     * \brief Rotates vector v using Rodrigues' formula: v' = v + 2 * u x (u x v + w * v)
     */
    Vec3 rotate(const Vec3& v) const {
        Vec3 u{x, y, z};
        Vec3 uv = u.cross(v);
        Vec3 uuv = u.cross(uv);
        return v + ((uv * w) + uuv) * 2.0;
    }
};