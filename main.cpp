#include "linking_number.h"
#include <iostream>

int main() {
    // Example usage
    std::cout << "Linking Number Computation Example\n";
    std::cout << "==================================\n\n";
    
    // Example 1: Two parallel line segments
    Vector3D p1_start(0.0, 0.0, 0.0);
    Vector3D p1_end(1.0, 0.0, 0.0);
    Vector3D p2_start(0.0, 1.0, 1.0);
    Vector3D p2_end(1.0, 1.0, 1.0);
    
    double linking_num1 = compute_linking_number_cartesian(p1_start, p1_end, p2_start, p2_end);
    
    std::cout << "Example 1 - Parallel segments:\n";
    std::cout << "Segment 1: ";
    p1_start.print();
    std::cout << " to ";
    p1_end.print();
    std::cout << "\n";
    std::cout << "Segment 2: ";
    p2_start.print();
    std::cout << " to ";
    p2_end.print();
    std::cout << "\n";
    std::cout << "Linking number: " << linking_num1 << "\n\n";
    
    // Example 2: Two perpendicular segments
    Vector3D p3_start(-1.0, 0.0, 0.0);
    Vector3D p3_end(1.0, 0.0, 0.0);
    Vector3D p4_start(0.0, -1.0, 1.0);
    Vector3D p4_end(0.0, 1.0, 1.0);
    
    double linking_num2 = compute_linking_number_cartesian(p3_start, p3_end, p4_start, p4_end);
    
    std::cout << "Example 2 - Perpendicular segments:\n";
    std::cout << "Segment 1: ";
    p3_start.print();
    std::cout << " to ";
    p3_end.print();
    std::cout << "\n";
    std::cout << "Segment 2: ";
    p4_start.print();
    std::cout << " to ";
    p4_end.print();
    std::cout << "\n";
    std::cout << "Linking number: " << linking_num2 << "\n\n";
    
    // Example 3: Using spherical coordinates to create segments
    double phi1 = M_PI / 4;    // 45 degrees
    double theta1 = 0;         // 0 degrees
    double phi2 = M_PI / 2;    // 90 degrees
    double theta2 = M_PI / 2;  // 90 degrees
    double length = 2.0;
    
    Vector3D center1(0.0, 0.0, 0.0);
    Vector3D center2(3.0, 0.0, 0.0);
    
    Vector3D u1 = spherical_to_cartesian(phi1, theta1);
    Vector3D u2 = spherical_to_cartesian(phi2, theta2);
    
    Vector3D p5_start = center1;
    Vector3D p5_end = center1 + u1 * length;
    Vector3D p6_start = center2;
    Vector3D p6_end = center2 + u2 * length;
    
    double linking_num3 = compute_linking_number_cartesian(p5_start, p5_end, p6_start, p6_end);
    
    std::cout << "Example 3 - Segments from spherical coordinates:\n";
    std::cout << "Segment 1: ";
    p5_start.print();
    std::cout << " to ";
    p5_end.print();
    std::cout << "\n";
    std::cout << "Segment 2: ";
    p6_start.print();
    std::cout << " to ";
    p6_end.print();
    std::cout << "\n";
    std::cout << "Linking number: " << linking_num3 << "\n";
    
    return 0;
}
