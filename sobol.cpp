#include "vdc_sobol_matrices52.h"

#include <iostream>

typedef unsigned int uint32;
typedef unsigned long long uint64;

uint64* c0_lower = 0;
uint64* c0_upper = 0;
uint64* c0_lower_inverse = 0;
uint64* c1_lower = 0;
uint64* c1_upper = 0;
uint64* c1_upper_inverse = 0;

// Van der Corput radical inverse in base 2 with 52 bits precision.
inline uint64 van_der_corput(uint64 bits, const uint64 scramble)
{
    bits = (bits << 32) | (bits >> 32);
    bits = ((bits & 0x0000ffff0000ffffULL) << 16) |
           ((bits & 0xffff0000ffff0000ULL) >> 16);
    bits = ((bits & 0x00ff00ff00ff00ffULL) << 8) |
           ((bits & 0xff00ff00ff00ff00ULL) >> 8);
    bits = ((bits & 0x0f0f0f0f0f0f0f0fULL) << 4) |
           ((bits & 0xf0f0f0f0f0f0f0f0ULL) >> 4);
    bits = ((bits & 0x3333333333333333ULL) << 2) |
           ((bits & 0xccccccccccccccccULL) >> 2);
    bits = ((bits & 0x5555555555555555ULL) << 1) |
           ((bits & 0xaaaaaaaaaaaaaaaaULL) >> 1);
    return (scramble ^ bits) >> (64 - 52); // Account for 52 bits precision.
}

// Sobol' radical inverse in base 2 with 52 bits precision.
inline uint64 sobol(uint64 i, const uint64 scramble)
{
    uint64 r = scramble >> (64 - 52);
    for (uint64 v = 1ULL << (52 - 1); i; i >>= 1, v ^= v >> 1)
        if (i & 1)
            r ^= v;
    return r;
}

// Return the index of the frame-th sample falling
// into the square elementary interval (px, py),
// without using look-up tables.
inline uint64 look_up(
    const uint32 m,
    uint32 frame,
    uint32 px,
    uint32 py,
    const uint64 scramble_x,
    const uint64 scramble_y)
{
    const uint32 m2 = m << 1;
    uint64 index = uint64(frame) << m2;

    // Note: the delta value only depends on frame
    // and m, thus it can be cached across multiple
    // function calls, if desired.
    uint64 delta = 0;
    for (uint32 c = 0; frame; frame >>= 1, ++c)
        if (frame & 1) // Add flipped column m + c + 1.
            delta ^= vdc_sobol_matrices[m - 1][c];

    px ^= scramble_x >> (64 - m);
    py ^= scramble_y >> (64 - m);
    uint64 b = ((uint64(px) << m) | py) ^ delta; // flipped b

    for (uint32 c = 0; b; b >>= 1, ++c)
        if (b & 1) // Add column 2 * m - c.
            index ^= vdc_sobol_matrices_inv[m - 1][c];

    return index;
}

// Compute the floating-point raster position for the
// frame-th sample falling into the pixel (px, py),
// and return the corresponding sample index.
// 2^m must be >= max{width, height}.
inline uint64 get_raster_pos(
    const uint32 m,
    const uint32 frame,
    const uint32 px,
    const uint32 py,
    const uint64 scramble_x,
    const uint64 scramble_y,
    double& rx,
    double& ry)
{
    const uint64 index = look_up(m, frame, px, py, scramble_x, scramble_y);
    rx = double(van_der_corput(index, scramble_x)) / (1ULL << (52 - m));
    ry = double(sobol(index, scramble_y)) / (1ULL << (52 - m));
    return index;
}

// Precompute the look-up table for the given
// frame number, with a resolution of 2^m * 2^m.
// c0_lower, c0_upper, c0_lower_inverse,
// c1_lower, c1_upper, and c1_upper_inverse
// denote uint64 arrays of size 2^m.
void precompute_lut(
    const uint32 m,
    const uint32 frame)
{
    const uint64 shift = uint64(frame) << (m << 1);
    for (uint64 i = 0; i < (1ULL << m); ++i) {
        const uint64 j = van_der_corput(i, 0);
        c0_lower[i] = j;
        c0_lower_inverse[j >> (52 - m)] = i;
        c0_upper[i] = van_der_corput((i << m) | shift, 0);

        c1_lower[i] = sobol(i, 0);
        const uint64 k = sobol((i << m) | shift, 0);
        c1_upper[i] = k;
        c1_upper_inverse[k >> (52 - m)] = i;
    }
}

// Look up the sample for the square elementary
// interval (px, py), fills the 32-bit integer
// coordinates (x, y) of the sample, and return
// the index of the sample.
inline uint64 look_up_lut(
    const uint32 m,
    const uint32 frame,
    uint32 px,
    uint32 py,
    const uint64 scramble_x,
    const uint64 scramble_y,
    uint64& x,
    uint64& y)
{
    px ^= scramble_x >> (64 - m);
    py ^= scramble_y >> (64 - m);
    const uint64 lower = c0_lower_inverse[px];
    const uint64 delta = uint64(py) ^ (c1_lower[lower] >> (52 - m));
    const uint64 upper = c1_upper_inverse[delta];
    const uint64 index = ((upper << m) | lower) | (uint64(frame) << (m << 1));
    x = c0_lower[lower] ^ c0_upper[upper] ^ (scramble_x >> (64 - 52));
    y = c1_lower[lower] ^ c1_upper[upper] ^ (scramble_y >> (64 - 52));
    return index;
}

// Compute the floating-point raster position for the
// frame-th sample falling into the pixel (px, py),
// and return the corresponding sample index.
// 2^m must be >= max{width, height} and
// precompute_lut must have been called
// for this frame before.
inline uint64 get_raster_pos_lut(
    const uint32 m,
    const uint32 frame,
    const uint32 px,
    const uint32 py,
    const uint64 scramble_x,
    const uint64 scramble_y,
    double& rx,
    double& ry)
{
    uint64 x, y;
    const uint64 index = look_up_lut(m, frame, px, py, scramble_x, scramble_y, x, y);
    rx = double(x) / (1ULL << (52 - m));
    ry = double(y) / (1ULL << (52 - m));
    return index;
}


int main(int argc, char** argv)
{
    using namespace std;

    // Run some tests and compare the results of
    // both algorithms to verify them.
    bool success = true;
    const uint32 max_m = 18;
    const uint32 max_samples = 1U << 8;
    for (uint32 m = 1; m < max_m; ++m)
    {
        cout << m << "..." << endl;

        const uint32 p = 1U << m; // 2^m

        c0_lower = new uint64[p];
        c0_upper = new uint64[p];
        c0_lower_inverse = new uint64[p];
        c1_lower = new uint64[p];
        c1_upper = new uint64[p];
        c1_upper_inverse = new uint64[p];

        // Scrambling parameters, these could be chosen randomly.
        // Only the 52 most significant bits of these are used.
        const uint64 scramble_x = 0xAABBCCDDEEFFEEFFull;
        const uint64 scramble_y = 0xFEFEFFEEABCDEFBCull;

        for (uint32 frame = 0; frame < max_samples; ++frame)
        {
            precompute_lut(m, frame);

            for (uint32 py = 0; py < p; py += 127) // arbitrary increase
            {
                for (uint32 px = 0; px < p; px += 111) // arbitrary increase
                {
                    uint64 x, y;
                    const uint64 index1 = look_up_lut(m, frame, px, py, scramble_x, scramble_y, x, y);
                    const uint64 index2 = look_up(m, frame, px, py, scramble_x, scramble_y);
                    double rx1, ry1;
                    const uint64 index3 = get_raster_pos_lut(m, frame, px, py, scramble_x, scramble_y, rx1, ry1);
                    double rx2, ry2;
                    const uint64 index4 = get_raster_pos(m, frame, px, py, scramble_x, scramble_y, rx2, ry2);

                    if (index1 != index2 ||
                        x >> (52 - m) != px ||
                        y >> (52 - m) != py ||
                        van_der_corput(index1, scramble_x) != x ||
                        sobol(index1, scramble_y) != y ||
                        index3 != index1 ||
                        index4 != index2 ||
                        rx1 < double(px) ||
                        rx1 >= double(px + 1) ||
                        ry1 < double(py) ||
                        ry1 >= double(py + 1) ||
                        rx1 != rx2 ||
                        ry1 != ry2)
                    {
                        cout << "failure detected: " << endl;
                        cout << "m = " << m << ", frame = " << frame << endl;
                        cout << "px = " << px << ", py = " << py << endl;
                        cout << "index: " << index1 << " vs. " << index2 << endl;
                        success = false;
                    }
                }
            }
        }

        delete [] c0_lower;
        delete [] c0_upper;
        delete [] c0_lower_inverse;
        delete [] c1_lower;
        delete [] c1_upper;
        delete [] c1_upper_inverse;
    }

    if (success)
    {
        cout << "All tests successful!" << endl;
        return 0;
    }
    else
    {
        cerr << "Failure detected!" << endl;
        return 1;
    }
}

