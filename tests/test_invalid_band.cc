// Verifies that an invalid band index terminates the process with a
// failure status rather than a success status.
#include <iostream>
#include "rad_planck_mean.h"

int main() {
    rad_planck_mean pm;
    double k = 0.0;
    double a = 0.0;
    // Band 1 does not exist for the Planck-mean model (single band 0).
    pm.get_k_a_oneband(k, a, 1, 1200.0, 101325.0, 0.0, 0.1, 0.1, 0.0, 0.0);
    std::cerr << "FAILED: get_k_a_oneband returned after an invalid band\n";
    return 0;
}
