# Linking Number Computation in C++

This project implements the linking number computation algorithm originally written in Python/JAX, translated to C++. The linking number is a topological invariant that measures how much two closed curves are linked together in 3D space.

## Overview

The program computes the linking number between two line segments in 3D space using the Gauss linking integral formula. This is particularly useful in computational topology, knot theory, and applications involving DNA topology or fluid dynamics.

## Algorithm

The implementation follows these steps:

1. **Input**: Four 3D points defining two line segments:
   - Segment 1: from `p_i` to `p_ii`
   - Segment 2: from `p_j` to `p_jj`

2. **Compute relative position vectors**:
   - `r_ij = p_i - p_j`
   - `r_ijj = p_i - p_jj`
   - `r_iij = p_ii - p_j`
   - `r_iijj = p_ii - p_jj`

3. **Calculate normal vectors** for each face of the quadrilateral formed by the four points:
   - `n1 = r_ij × r_ijj` (normalized)
   - `n2 = r_ijj × r_iijj` (normalized)
   - `n3 = r_iijj × r_iij` (normalized)
   - `n4 = r_iij × r_ij` (normalized)

4. **Compute the linking number**:
   ```
   LN = -1/(4π) * |arcsin(n1·n2) + arcsin(n2·n3) + arcsin(n3·n4) + arcsin(n4·n1)|
   ```

## Features

- **Vector3D class**: Complete 3D vector implementation with:
  - Basic arithmetic operations (+, -, *)
  - Cross product and dot product
  - Vector normalization
  - Pretty printing

- **Robust computation**: Includes tolerance handling and value clipping to ensure numerical stability

- **Helper functions**: Conversion from spherical to Cartesian coordinates

- **Example usage**: Three different test cases demonstrating various geometric configurations

## Building the Project

### Using CMake (Recommended)
```bash
mkdir build
cd build
cmake ..
make
./linking_number
```

### Using Make
```bash
make
./linking_number
```

### Manual compilation
```bash
g++ -std=c++17 -Wall -Wextra -O2 -o linking_number linking_number.cpp -lm
./linking_number
```

## Requirements

- C++17 compatible compiler (GCC, Clang, or MSVC)
- CMake 3.10+ (optional, for CMake build)
- Math library (automatically linked)

## Examples

The program includes three example computations:

1. **Parallel segments**: Two parallel line segments with zero linking number
2. **Perpendicular segments**: Two perpendicular segments with non-zero linking number
3. **Spherical coordinates**: Segments created using spherical coordinate parameterization

## Mathematical Background

The linking number is a topological invariant that quantifies the linking of two closed curves in 3D space. For line segments, this implementation uses the Gauss linking integral formula, which involves:

- Computing the solid angle subtended by the quadrilateral formed by the segment endpoints
- Normalizing by 4π to get the linking number
- Taking the absolute value as specified in the original algorithm

## Applications

This algorithm is useful in:
- Computational topology and knot theory
- DNA topology analysis
- Fluid dynamics (vortex linking)
- Computer graphics and geometric modeling
- Polymer physics

## Notes

- The algorithm includes numerical tolerance handling to avoid division by zero and invalid arcsin arguments
- The implementation closely follows the original Python/JAX version for consistency
- All computations use double precision floating-point arithmetic
