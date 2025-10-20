#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <random>
#include <iomanip>
#include <chrono>
#include <array>

#ifdef _OPENMP
#include <omp.h>
#endif

class Vector3D {
public:
    double x, y, z;
    
    Vector3D() : x(0), y(0), z(0) {}
    Vector3D(double x, double y, double z) : x(x), y(y), z(z) {}
    
    Vector3D operator-(const Vector3D& other) const {
        return Vector3D(x - other.x, y - other.y, z - other.z);
    }
    
    Vector3D operator+(const Vector3D& other) const {
        return Vector3D(x + other.x, y + other.y, z + other.z);
    }
    
    Vector3D operator*(double scalar) const {
        return Vector3D(x * scalar, y * scalar, z * scalar);
    }
    
    // Cross product
    Vector3D cross(const Vector3D& other) const {
        return Vector3D(
            y * other.z - z * other.y,
            z * other.x - x * other.z,
            x * other.y - y * other.x
        );
    }
    
    // Dot product
    double dot(const Vector3D& other) const {
        return x * other.x + y * other.y + z * other.z;
    }
    
    // Magnitude/norm
    double norm() const {
        return std::sqrt(x * x + y * y + z * z);
    }
    
    // Normalize vector
    Vector3D normalize(double tol = 1e-6) const {
        double magnitude = norm();
        if (magnitude < tol) {
            return Vector3D(0, 0, 0);
        }
        return Vector3D(x / magnitude, y / magnitude, z / magnitude);
    }
    
    void print() const {
        std::cout << "(" << x << ", " << y << ", " << z << ")";
    }
};

double compute_linking_number_cartesian(const Vector3D& p_i, const Vector3D& p_ii, 
                                       const Vector3D& p_j, const Vector3D& p_jj) {
    // Compute the relative position vectors
    Vector3D r_ij = p_i - p_j;
    Vector3D r_ijj = p_i - p_jj;
    Vector3D r_iij = p_ii - p_j;
    Vector3D r_iijj = p_ii - p_jj;
    
    const double tol = 1e-6;
    
    // Compute normal vectors for each face of the quadrilateral
    Vector3D n1 = r_ij.cross(r_ijj);
    n1 = n1 * (1.0 / (n1.norm() + tol));
    
    Vector3D n2 = r_ijj.cross(r_iijj);
    n2 = n2 * (1.0 / (n2.norm() + tol));
    
    Vector3D n3 = r_iijj.cross(r_iij);
    n3 = n3 * (1.0 / (n3.norm() + tol));
    
    Vector3D n4 = r_iij.cross(r_ij);
    n4 = n4 * (1.0 / (n4.norm() + tol));
    
    const double clip_tol = 0.0;
    
    // Compute dot products and clip them to valid range for arcsin
    auto clip = [](double value, double min_val, double max_val) {
        return std::max(min_val, std::min(max_val, value));
    };
    
    double dot1 = clip(n1.dot(n2), -1.0 + clip_tol, 1.0 - clip_tol);
    double dot2 = clip(n2.dot(n3), -1.0 + clip_tol, 1.0 - clip_tol);
    double dot3 = clip(n3.dot(n4), -1.0 + clip_tol, 1.0 - clip_tol);
    double dot4 = clip(n4.dot(n1), -1.0 + clip_tol, 1.0 - clip_tol);
    
    // Compute the linking number using the Gauss linking integral formula
    double sum_angles = std::asin(dot1) + std::asin(dot2) + std::asin(dot3) + std::asin(dot4);
    
    return 1.0 / (4.0 * M_PI) * (sum_angles);
}

// Helper function to create unit vectors from spherical coordinates
Vector3D spherical_to_cartesian(double phi, double theta) {
    return Vector3D(
        std::sin(phi) * std::cos(theta),
        std::sin(phi) * std::sin(theta),
        std::cos(phi)
    );
}

// Generate N random rods (each of length `length`) and return as an (N,6) array
// Each entry is {ax, ay, az, bx, by, bz}
std::vector<std::array<double,6>> generate_random_rods_array(int N, double length, double box_size, unsigned long long seed) {
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> pos_dist(-box_size, box_size);
    std::uniform_real_distribution<double> unit_dist(0.0, 1.0);
    std::uniform_real_distribution<double> ang_dist(0.0, 2.0 * M_PI);

    std::vector<std::array<double,6>> arr;
    if (N > 0) arr.reserve(static_cast<size_t>(N));

    for (int i = 0; i < N; ++i) {
        double cx = pos_dist(rng);
        double cy = pos_dist(rng);
        double cz = pos_dist(rng);
        // Uniform random direction on sphere
        double z = 2.0 * unit_dist(rng) - 1.0; // cos(phi)
        double theta = ang_dist(rng);
        double sinphi = std::sqrt(std::max(0.0, 1.0 - z * z));
        Vector3D u(sinphi * std::cos(theta), sinphi * std::sin(theta), z);

        Vector3D center(cx, cy, cz);
        Vector3D a = center - u * (length * 0.5);
        Vector3D b = center + u * (length * 0.5);

        std::array<double,6> entry = {a.x, a.y, a.z, b.x, b.y, b.z};
        arr.push_back(entry);
    }

    return arr;
}

// Compute shortest distance between two 3D segments p1q1 and p2q2
static double segment_segment_distance(const Vector3D& p1, const Vector3D& q1,
                                       const Vector3D& p2, const Vector3D& q2) {
    const double EPS = 1e-12;
    Vector3D u = q1 - p1;
    Vector3D v = q2 - p2;
    Vector3D w = p1 - p2;
    double a = u.dot(u);
    double b = u.dot(v);
    double c = v.dot(v);
    double d = u.dot(w);
    double e = v.dot(w);
    double D = a * c - b * b;

    double sN, sD = D;
    double tN, tD = D;

    if (D < EPS) {
        // Lines almost parallel
        sN = 0.0; sD = 1.0; // force s = 0
        tN = e;  tD = c;    // t = e/c
    } else {
        sN = (b * e - c * d);
        tN = (a * e - b * d);
        if (sN < 0.0) { sN = 0.0; tN = e; tD = c; }
        else if (sN > sD) { sN = sD; tN = e + b; tD = c; }
    }

    if (tN < 0.0) {
        tN = 0.0;
        if (-d < 0.0) { sN = 0.0; sD = 1.0; }
        else if (-d > a) { sN = sD; }
        else { sN = -d; sD = a; }
    } else if (tN > tD) {
        tN = tD;
        if ((-d + b) < 0.0) { sN = 0.0; sD = 1.0; }
        else if ((-d + b) > a) { sN = sD; }
        else { sN = (-d + b); sD = a; }
    }

    double sc = (std::abs(sN) < EPS ? 0.0 : sN / sD);
    double tc = (std::abs(tN) < EPS ? 0.0 : tN / tD);

    Vector3D dP = w + u * sc - v * tc;  // =  p1 + sc*u - (p2 + tc*v)
    return dP.norm();
}

// Compute total sum of |linking| with pruning for rods given as an (N,6) array
// cutoff <= 0 disables pruning. Optionally returns number of evaluated pairs.
double pairwise_abs_linking_sum_with_cutoff(
    const std::vector<std::array<double,6>>& rods_array,
    double cutoff,
    long long* out_pairs,
    int num_threads
) {
    const int N = static_cast<int>(rods_array.size());

    // Fast path: no cutoff => no precompute, no extra checks
    if (!(cutoff > 0.0)) {
        long long pairs = 0;
        double total_abs_linking = 0.0;
        #pragma omp parallel for schedule(dynamic,64) reduction(+: total_abs_linking, pairs) if(num_threads>1)
        for (int i = 0; i < N; ++i) {
            Vector3D ai{rods_array[i][0], rods_array[i][1], rods_array[i][2]};
            Vector3D bi{rods_array[i][3], rods_array[i][4], rods_array[i][5]};
            for (int j = i + 1; j < N; ++j) {
                Vector3D aj{rods_array[j][0], rods_array[j][1], rods_array[j][2]};
                Vector3D bj{rods_array[j][3], rods_array[j][4], rods_array[j][5]};
                double lk = compute_linking_number_cartesian(ai, bi, aj, bj);
                total_abs_linking += std::abs(lk);
                ++pairs;
            }
        }
        if (out_pairs) *out_pairs = pairs;
        return total_abs_linking;
    }

    // Cutoff path: build centers once and prune
    std::vector<Vector3D> centers;
    centers.reserve(N);
    for (int i = 0; i < N; ++i) {
        centers.push_back(Vector3D{
            0.5 * (rods_array[i][0] + rods_array[i][3]),
            0.5 * (rods_array[i][1] + rods_array[i][4]),
            0.5 * (rods_array[i][2] + rods_array[i][5])
        });
    }

    long long pairs = 0;
    double total_abs_linking = 0.0;

    #pragma omp parallel for schedule(dynamic,64) reduction(+: total_abs_linking, pairs) if(num_threads>1)
    for (int i = 0; i < N; ++i) {
        Vector3D ai{rods_array[i][0], rods_array[i][1], rods_array[i][2]};
        Vector3D bi{rods_array[i][3], rods_array[i][4], rods_array[i][5]};
        for (int j = i + 1; j < N; ++j) {
            double centerLower = std::max(0.0, (centers[i] - centers[j]).norm() - 1.0);
            if (centerLower > cutoff) continue;
            Vector3D aj{rods_array[j][0], rods_array[j][1], rods_array[j][2]};
            Vector3D bj{rods_array[j][3], rods_array[j][4], rods_array[j][5]};
            double sdist = segment_segment_distance(ai, bi, aj, bj);
            if (sdist > cutoff) continue;
            double lk = compute_linking_number_cartesian(ai, bi, aj, bj);
            total_abs_linking += std::abs(lk);
            ++pairs;
        }
    }

    if (out_pairs) *out_pairs = pairs;
    return total_abs_linking;
}

int main(int argc, char** argv) {
    // Default parameters
    int N = 50;                     // number of rods
    const double length = 1.0;      // rod length (fixed = 1)
    double box_size = 5.0;          // half-size of the cubic domain for centers
    unsigned long long seed = 123456789ULL; // fixed seed for reproducibility
    double cutoff = -1.0;           // single cutoff (<=0 disables)

    // Sweep parameters (if provided): [cutoff_start, cutoff_end, cutoff_step]
    bool do_sweep = false;
    double cut_start = 0.0, cut_end = 0.0, cut_step = 0.0;

    // Simple CLI parsing: ./linking_number [N] [box_size] [seed] [cutoff]
    // Or: ./linking_number [N] [box_size] [seed] [cut_start] [cut_end] [cut_step]
    if (argc > 1) { try { N = std::max(0, std::stoi(argv[1])); } catch (...) { std::cerr << "Invalid N value: " << argv[1] << "\n"; return 1; } }
    if (argc > 2) { try { box_size = std::stod(argv[2]); } catch (...) { std::cerr << "Invalid box_size value: " << argv[2] << "\n"; return 1; } }
    if (argc > 3) { try { seed = std::stoull(argv[3]); } catch (...) { std::cerr << "Invalid seed value: " << argv[3] << "\n"; return 1; } }

    if (argc == 5) { // single cutoff
        try { cutoff = std::stod(argv[4]); } catch (...) { std::cerr << "Invalid cutoff value: " << argv[4] << "\n"; return 1; }
    } else if (argc >= 7) { // sweep
        try {
            cut_start = std::stod(argv[4]);
            cut_end   = std::stod(argv[5]);
            cut_step  = std::stod(argv[6]);
            do_sweep = (cut_step > 0 && cut_end >= cut_start);
        } catch (...) {
            std::cerr << "Invalid sweep values. Usage: [cut_start] [cut_end] [cut_step]\n"; return 1;
        }
    }

    if (N <= 0) { std::cerr << "N must be positive.\n"; return 1; }

#ifdef _OPENMP
    int num_threads = omp_get_max_threads();
#else
    int num_threads = 1;
#endif

    std::cout << "Linking Number Pairwise Sum\n";
    std::cout << "---------------------------\n";
    std::cout << "Number of rods: " << N << "\n";
    std::cout << "Rod length (fixed): " << length << "\n";
    std::cout << "Box half-size: " << box_size << " (centers uniformly in [-box_size,box_size]^3)\n";
    std::cout << "RNG seed: " << seed << "\n";
    if (do_sweep) {
        std::cout << "Cutoff sweep: start=" << cut_start << ", end=" << cut_end << ", step=" << cut_step << "\n";
    } else {
        std::cout << "Cutoff: " << cutoff << (cutoff > 0.0 ? " (enabled)" : " (disabled)") << "\n";
    }
    std::cout << "Using threads: " << num_threads << "\n\n";

    // Generate rods once (array form) for consistent benchmarking
    auto rods_array = generate_random_rods_array(N, length, box_size, seed);

    if (!do_sweep) {
        long long pairs = 0;
        auto t0 = std::chrono::high_resolution_clock::now();
        double total_abs_linking = pairwise_abs_linking_sum_with_cutoff(rods_array, cutoff, &pairs, num_threads);
        auto t1 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = t1 - t0;

        std::cout << std::fixed << std::setprecision(8);
        std::cout << "Pairs evaluated: " << pairs << "\n";
        std::cout << "Total sum of |linking|: " << total_abs_linking << "\n";
        std::cout << "Average |linking| per pair: " << (pairs ? total_abs_linking / pairs : 0.0) << "\n";
        std::cout << "Elapsed time (s): " << elapsed.count() << "\n";
    } else {
        // Header
        std::cout << std::left
                  << std::setw(12) << "cutoff"
                  << std::setw(18) << "pairs"
                  << std::setw(22) << "sum_|linking|"
                  << std::setw(16) << "time_s"
                  << "\n";
        std::cout << std::string(12+18+22+16, '-') << "\n";

        for (double c = cut_start; c <= cut_end + 1e-12; c += cut_step) {
            long long pairs = 0;
            auto t0 = std::chrono::high_resolution_clock::now();
            double total_abs_linking = pairwise_abs_linking_sum_with_cutoff(rods_array, c, &pairs, num_threads);
            auto t1 = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = t1 - t0;

            std::cout << std::left
                      << std::setw(12) << std::setprecision(6) << std::fixed << c
                      << std::setw(18) << pairs
                      << std::setw(22) << std::setprecision(8) << total_abs_linking
                      << std::setw(16) << std::setprecision(6) << elapsed.count()
                      << "\n";
        }
    }

    return 0;
}
