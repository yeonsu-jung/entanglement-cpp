#ifndef LINKING_NUMBER_H
#define LINKING_NUMBER_H

#include <cmath>

/**
 * @brief A 3D vector class with basic operations for geometric computations
 */
class Vector3D {
public:
    double x, y, z;
    
    // Constructors
    Vector3D();
    Vector3D(double x, double y, double z);
    
    // Arithmetic operations
    Vector3D operator-(const Vector3D& other) const;
    Vector3D operator+(const Vector3D& other) const;
    Vector3D operator*(double scalar) const;
    
    // Vector operations
    Vector3D cross(const Vector3D& other) const;
    double dot(const Vector3D& other) const;
    double norm() const;
    Vector3D normalize(double tol = 1e-6) const;
    
    // Utility
    void print() const;
};

/**
 * @brief Compute the linking number between two line segments in 3D space
 * 
 * This function implements the Gauss linking integral formula to compute
 * the linking number between two line segments defined by their endpoints.
 * 
 * @param p_i Start point of the first segment
 * @param p_ii End point of the first segment
 * @param p_j Start point of the second segment
 * @param p_jj End point of the second segment
 * @return The linking number (non-negative value)
 */
double compute_linking_number_cartesian(const Vector3D& p_i, const Vector3D& p_ii, 
                                       const Vector3D& p_j, const Vector3D& p_jj);

/**
 * @brief Convert spherical coordinates to Cartesian unit vector
 * 
 * @param phi Polar angle (0 to π)
 * @param theta Azimuthal angle (0 to 2π)
 * @return Unit vector in Cartesian coordinates
 */
Vector3D spherical_to_cartesian(double phi, double theta);

#endif // LINKING_NUMBER_H
