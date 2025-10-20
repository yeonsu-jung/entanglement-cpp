#include "linking_number.h"
#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include <random>
#include <iomanip>
#include <chrono>

#ifdef _OPENMP
#include <omp.h>
#endif

// ---- Vector3D method implementations ----
Vector3D::Vector3D() : x(0), y(0), z(0) {}
Vector3D::Vector3D(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

Vector3D Vector3D::operator-(const Vector3D& other) const { return Vector3D(x - other.x, y - other.y, z - other.z); }
Vector3D Vector3D::operator+(const Vector3D& other) const { return Vector3D(x + other.x, y + other.y, z + other.z); }
Vector3D Vector3D::operator*(double scalar) const { return Vector3D(x * scalar, y * scalar, z * scalar); }

Vector3D Vector3D::cross(const Vector3D& other) const {
    return Vector3D(
        y * other.z - z * other.y,
        z * other.x - x * other.z,
        x * other.y - y * other.x
    );
}

double Vector3D::dot(const Vector3D& other) const { return x * other.x + y * other.y + z * other.z; }

double Vector3D::norm() const { return std::sqrt(x * x + y * y + z * z); }

Vector3D Vector3D::normalize(double tol) const {
    double m = norm();
    if (m < tol) return Vector3D(0,0,0);
    return Vector3D(x / m, y / m, z / m);
}

void Vector3D::print() const { std::cout << "(" << x << ", " << y << ", " << z << ")"; }

// ---- API implementations ----

double compute_linking_number_cartesian(const Vector3D& p_i, const Vector3D& p_ii,
                                        const Vector3D& p_j, const Vector3D& p_jj) {
    Vector3D r_ij  = p_i  - p_j;
    Vector3D r_ijj = p_i  - p_jj;
    Vector3D r_iij = p_ii - p_j;
    Vector3D r_iijj= p_ii - p_jj;

    const double tol = 1e-6;

    Vector3D n1 = r_ij.cross(r_ijj);   n1 = n1 * (1.0 / (n1.norm() + tol));
    Vector3D n2 = r_ijj.cross(r_iijj); n2 = n2 * (1.0 / (n2.norm() + tol));
    Vector3D n3 = r_iijj.cross(r_iij); n3 = n3 * (1.0 / (n3.norm() + tol));
    Vector3D n4 = r_iij.cross(r_ij);   n4 = n4 * (1.0 / (n4.norm() + tol));

    auto clip = [](double v, double lo, double hi){ return std::max(lo, std::min(hi, v)); };
    const double clip_tol = 0.0;

    double dot1 = clip(n1.dot(n2), -1.0 + clip_tol, 1.0 - clip_tol);
    double dot2 = clip(n2.dot(n3), -1.0 + clip_tol, 1.0 - clip_tol);
    double dot3 = clip(n3.dot(n4), -1.0 + clip_tol, 1.0 - clip_tol);
    double dot4 = clip(n4.dot(n1), -1.0 + clip_tol, 1.0 - clip_tol);

    double sum_angles = std::asin(dot1) + std::asin(dot2) + std::asin(dot3) + std::asin(dot4);
    return 1.0 / (4.0 * M_PI) * sum_angles;
}

Vector3D spherical_to_cartesian(double phi, double theta) {
    return Vector3D(std::sin(phi) * std::cos(theta),
                    std::sin(phi) * std::sin(theta),
                    std::cos(phi));
}

#ifndef ENT_NO_MAIN
int main(int, char**) {
    std::cerr << "This build is intended as a library. Define ENT_NO_MAIN=0 to enable CLI.\n";
    return 0;
}
#endif
