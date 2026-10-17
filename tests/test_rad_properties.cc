// Numerical property and regression tests for the radlib radiation models.
//
// These tests check structural invariants of the Planck-mean, WSGG, and
// RCSLW models and exact analytic behavior of the interpolation kernels.
// They do not re-derive the fitted correlations; instead they pin down
// reference values as regressions and verify physical constraints
// (partition of unity of the gas weights, non-negative absorption,
// uniform soot addition across bands, temperature clamping).
#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "multilinear_interpolation.h"
#include "rad.h"
#include "rad_planck_mean.h"
#include "rad_rcslw.h"
#include "rad_wsgg.h"

namespace {

void require(bool ok, const char* msg) {
    if (!ok) {
        std::fprintf(stderr, "FAILED: %s\n", msg);
        std::exit(EXIT_FAILURE);
    }
}

bool near(double a, double b, double rtol, double atol) {
    const double scale = rtol * (std::fabs(a) > std::fabs(b) ? std::fabs(a) : std::fabs(b));
    const double tol = (scale > atol ? scale : atol);
    return std::fabs(a - b) <= tol;
}

}  // namespace

int main() {
    // ------------------------------------------------------------------
    // Multilinear interpolation kernels: exactness on linear fields.
    // ------------------------------------------------------------------
    {
        const double xs[3] = {0.0, 1.0, 2.5};
        const double fs[3] = {1.0, 3.0, 6.0};          // f(x) = 1 + 2x on the grid
        require(near(LI_1D(3, xs, fs, 0.625), 1.0 + 2.0 * 0.625, 1e-14, 1e-14),
                "LI_1D reproduces a linear function");
        // Outside the grid the kernels continue the end interval linearly
        // (rad_rcslw clamps its own table inputs before calling them).
        require(near(LI_1D(3, xs, fs, 3.5), 1.0 + 2.0 * 3.5, 1e-14, 1e-14),
                "LI_1D continues the end interval above the grid");
        require(near(LI_1D(3, xs, fs, -5.0), 1.0 + 2.0 * (-5.0), 1e-14, 1e-14),
                "LI_1D continues the first interval below the grid");

        const double gx[3] = {0.0, 1.0, 2.0};
        const double gy[3] = {0.0, 1.0, 3.0};
        std::vector<double> f2(9);
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                f2[i * 3 + j] = 1.0 + 2.0 * gx[i] - 0.5 * gy[j];
        require(near(LI_2D(3, 3, gx, gy, f2.data(), 0.625, 2.3), 1.0 + 2.0 * 0.625 - 0.5 * 2.3, 1e-13, 1e-13),
                "LI_2D reproduces a bilinear (linear) function");

        const double gz[3] = {0.0, 0.5, 2.0};
        std::vector<double> f3(27);
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                for (int k = 0; k < 3; ++k)
                    f3[(i * 3 + j) * 3 + k] = 2.0 + 1.5 * gx[i] - gy[j] + 0.25 * gz[k];
        require(near(LI_3D(3, 3, 3, gx, gy, gz, f3.data(), 1.3, 0.2, 0.9),
                     2.0 + 1.5 * 1.3 - 0.2 + 0.25 * 0.9, 1e-13, 1e-13),
                "LI_3D reproduces a trilinear (linear) function");

        const double gw[3] = {0.0, 1.0, 4.0};
        std::vector<double> f4(81);
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                for (int k = 0; k < 3; ++k)
                    for (int l = 0; l < 3; ++l)
                        f4[((i * 3 + j) * 3 + k) * 3 + l] = -1.0 + gx[i] + 2.0 * gy[j] - gz[k] + 0.5 * gw[l];
        require(near(LI_4D(3, 3, 3, 3, gx, gy, gz, gw, f4.data(), 1.7, 2.9, 1.1, 2.5),
                     -1.0 + 1.7 + 2.0 * 2.9 - 1.1 + 0.5 * 2.5, 1e-13, 1e-13),
                "LI_4D reproduces a multilinear function");
    }

    // ------------------------------------------------------------------
    // Planck-mean model.
    // ------------------------------------------------------------------
    {
        rad_planck_mean pm;
        std::vector<double> k, a;

        // Soot-only state: k = 1817 * fvsoot * T exactly for in-range T.
        pm.get_k_a(k, a, 1200.0, 101325.0, 1.0e-6, 0.0, 0.0, 0.0, 0.0);
        require(near(k[0], 1817.0 * 1.0e-6 * 1200.0, 1e-14, 1e-14),
                "Planck soot-only absorption matches F_s*fv*T");
        require(a[0] == 1.0, "Planck weight is unity");

        // Below 300 K the curve fits are clamped, so 200 K == 299 K.
        std::vector<double> k1, a1, k2, a2;
        pm.get_k_a(k1, a1, 200.0, 101325.0, 1.0e-7, 0.1, 0.1, 0.05, 0.01);
        pm.get_k_a(k2, a2, 299.0, 101325.0, 1.0e-7, 0.1, 0.1, 0.05, 0.01);
        require(k1[0] == k2[0], "Planck values clamp below 300 K");
        pm.get_k_a(k1, a1, 2500.0, 101325.0, 0.0, 0.1, 0.1, 0.05, 0.01);
        pm.get_k_a(k2, a2, 4000.0, 101325.0, 0.0, 0.1, 0.1, 0.05, 0.01);
        require(k1[0] == k2[0], "Planck values clamp above 2500 K");

        // Linearity in mole fraction: k(xCO2) is proportional to xCO2.
        pm.get_k_a(k1, a1, 1500.0, 101325.0, 0.0, 0.0, 0.1, 0.0, 0.0);
        pm.get_k_a(k2, a2, 1500.0, 101325.0, 0.0, 0.0, 0.2, 0.0, 0.0);
        require(near(k2[0], 2.0 * k1[0], 1e-14, 1e-14),
                "Planck absorption is linear in xCO2");

        // Positivity across a representative combustion sweep.
        for (double T = 300.0; T <= 2500.0; T += 275.0) {
            pm.get_k_a(k, a, T, 101325.0, 1.0e-7, 0.15, 0.10, 0.02, 0.005);
            require(k[0] > 0.0, "Planck absorption stays positive");
        }

        // Regression anchors (values from the reference implementation).
        pm.get_k_a(k, a, 1500.0, 101325.0, 1.0e-6, 0.10, 0.10, 0.0, 0.0);
        require(near(k[0], 4.1956217992230442, 1e-13, 1e-13),
                "Planck-mean anchor at 1500 K");
        require(a[0] == 1.0, "Planck-mean weight anchor at 1500 K");
    }

    // ------------------------------------------------------------------
    // WSGG model (Bordbar 2020): four gray gases + one clear gas.
    // ------------------------------------------------------------------
    {
        rad_wsgg wsgg;
        std::vector<double> k, a;

        const double T_values[3] = {600.0, 1200.0, 2000.0};
        // States chosen to cover the Mr interpolation branches:
        // mid-range, pure CO2, dilute H2O (Mr < 0.01), H2O-rich (Mr > 4),
        // and effectively pure H2O (zero CO2).
        const double xH2O_values[5] = {0.05,  0.0,   0.001, 0.30,  0.20};
        const double xCO2_values[5] = {0.10,  0.10,  0.20,  0.03,  0.0 };

        for (int s = 0; s < 5; ++s) {
            for (int t = 0; t < 3; ++t) {
                wsgg.get_k_a(k, a, T_values[t], 101325.0, 0.0,
                             xH2O_values[s], xCO2_values[s], 0.0, 0.0);
                require((int)k.size() == 5 && (int)a.size() == 5, "WSGG returns nGGa entries");
                double sum = 0.0;
                for (int i = 0; i < 5; ++i) {
                    sum += a[i];
                    require(a[i] > -1e-12, "WSGG weights are non-negative");
                    require(k[i] >= 0.0, "WSGG absorption is non-negative");
                }
                require(near(sum, 1.0, 1e-12, 1e-12), "WSGG weights sum to one");
                require(k[0] == 0.0, "WSGG clear gas has zero absorption without soot");

                // Soot adds the same gray increment to every band.
                std::vector<double> ks, as;
                const double fv = 1.0e-6;
                wsgg.get_k_a(ks, as, T_values[t], 101325.0, fv,
                             xH2O_values[s], xCO2_values[s], 0.0, 0.0);
                const double expected = 1817.0 * fv * T_values[t];
                for (int i = 0; i < 5; ++i)
                    require(near(ks[i] - k[i], expected, 1e-12, 1e-12),
                            "WSGG soot increments every band by F_s*fv*T");
            }
        }

        // Temperature clamping at the correlation bound of 2400 K.
        wsgg.get_k_a(k, a, 2400.0, 101325.0, 0.0, 0.1, 0.1, 0.0, 0.0);
        std::vector<double> k2, a2;
        wsgg.get_k_a(k2, a2, 3000.0, 101325.0, 0.0, 0.1, 0.1, 0.0, 0.0);
        for (int i = 0; i < 5; ++i)
            require(k[i] == k2[i] && a[i] == a2[i], "WSGG values clamp above 2400 K");

        // Regression anchors (values from the reference implementation).
        wsgg.get_k_a(k, a, 1200.0, 101325.0, 0.0, 0.10, 0.10, 0.0, 0.0);
        require(near(k[2], 0.14219203798006155, 1e-13, 1e-13),
                "WSGG absorption anchor at 1200 K");
        require(near(a[1], 0.30968436702830071, 1e-13, 1e-13),
                "WSGG weight anchor at 1200 K");
    }

    // ------------------------------------------------------------------
    // RCSLW model (requires the generated binary ALBDF tables; run
    // "cmake --install" before ctest so the installed data dir is filled).
    // ------------------------------------------------------------------
    {
        rad_rcslw rc(4, 1200.0, 101325.0, 1.0e-6, 0.10, 0.10, 0.0);
        std::vector<double> k, a;
        rc.get_k_a(k, a, 1200.0, 101325.0, 1.0e-6, 0.10, 0.10, 0.0, 0.0);
        require((int)k.size() == 5 && (int)a.size() == 5, "RCSLW returns nGGa entries");
        double sum = 0.0;
        for (int i = 0; i < 5; ++i) {
            sum += a[i];
            require(k[i] >= 0.0, "RCSLW absorption is non-negative");
        }
        require(near(sum, 1.0, 1e-3, 1e-3), "RCSLW weights approximately sum to one");
        require(k[0] == 0.0, "RCSLW clear gas has zero absorption");

        // Gray-gas absorption grows with soot loading. Soot is fixed at
        // construction in the RCSLW model, so compare two instances.
        rad_rcslw rc_soot(4, 1200.0, 101325.0, 1.0e-5, 0.10, 0.10, 0.0);
        std::vector<double> ks, as;
        rc_soot.get_k_a(ks, as, 1200.0, 101325.0, 1.0e-5, 0.10, 0.10, 0.0, 0.0);
        double gray_plain = 0.0;
        double gray_sooty = 0.0;
        for (int i = 1; i < 5; ++i) {
            require(ks[i] >= k[i], "RCSLW gray absorption does not drop with soot");
            gray_plain += k[i];
            gray_sooty += ks[i];
        }
        require(gray_sooty > gray_plain, "RCSLW total gray absorption increases with soot");

        // Regression anchors (values from the reference implementation and
        // its generated ALBDF tables).
        require(near(k[1], 1.2567417656281081, 1e-12, 1e-12),
                "RCSLW gray-gas absorption anchor");
        require(near(k[4], 11.692078359380767, 1e-12, 1e-12),
                "RCSLW strongest-band absorption anchor");
        require(near(a[0], 0.00050505488390405185, 1e-12, 1e-15),
                "RCSLW clear-gas weight anchor");
    }

    std::printf("radlib property tests passed\n");
    return EXIT_SUCCESS;
}
