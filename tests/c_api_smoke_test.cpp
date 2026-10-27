// Smoke test for the 0.1.0 C API. Exercises the SHARED library — the same
// binary the Python bindings load.
//
// Covers the cache lifecycle with the status of every create rejection, the
// cached computes with their explicit sizes and the usage-code order
// (2 > 4 > 3 > the value codes), NULL and negative-size arguments, a valid and
// a separated shape through every output, the one-shot calls with their
// nullable trailing status, and the static status-message table.
#include "tqs_parameterization.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <numbers>
#include <vector>

namespace {
int failures = 0;

void check(const bool ok, const char* label) {
    if (!ok) {
        std::printf("FAIL: %s\n", label);
        ++failures;
    }
}

constexpr int n_theta = 64;
constexpr int n_z = 101;

const std::vector<double> sym_necked{1.6, -0.6, 1.0, 0.0, 0.0};
const std::vector<double> separated{3.0, -1.0, 1.0, 0.0, 0.0};

// Uniform theta nodes over [0, pi], endpoints pinned: under fast-math the
// product (n-1)*pi/(n-1) can land one ulp above pi, which the domain check
// would then reject.
std::vector<double> theta_grid(const int n) {
    std::vector<double> thetas(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        thetas[static_cast<std::size_t>(i)] =
                static_cast<double>(i) * std::numbers::pi / static_cast<double>(n - 1);
    }
    thetas[0] = 0.0;
    thetas[static_cast<std::size_t>(n - 1)] = std::numbers::pi;
    return thetas;
}

bool all_finite_positive(const std::vector<double>& v) {
    for (const double x : v) {
        if (!std::isfinite(x) || x <= 0.0) return false;
    }
    return true;
}

bool all_zero(const std::vector<double>& v) {
    for (const double x : v) {
        if (x != 0.0) return false;
    }
    return true;
}

void test_create(const std::vector<double>& thetas) {
    int32_t status = -1;
    tqs_param_cache_t* c = tqs_param_cache_create(4, thetas.data(), n_theta, &status);
    check(c == nullptr && status == TQS_ERROR_INVALID_INIT, "create max_params 4 -> NULL, 5");
    c = tqs_param_cache_create(6, thetas.data(), n_theta, &status);
    check(c == nullptr && status == TQS_ERROR_TOO_MANY_PARAMS, "create max_params 6 -> NULL, 1");
    c = tqs_param_cache_create(5, nullptr, 0, &status);
    check(c == nullptr && status == TQS_ERROR_INVALID_GRID, "create with no thetas -> NULL, 3");
    c = tqs_param_cache_create(5, nullptr, -3, &status);
    check(c == nullptr && status == TQS_ERROR_INVALID_GRID, "create with negative n_thetas -> NULL, 3");
    c = tqs_param_cache_create(5, thetas.data(), n_theta, nullptr);
    check(c != nullptr, "create with a NULL status pointer");
    tqs_param_cache_destroy(c);
    tqs_param_cache_destroy(nullptr);
}

void test_usage_codes(const tqs_param_cache_t* cache) {
    std::vector<double> radii(n_theta), dr(n_theta), z(n_z), rho(n_z), drho(n_z);
    double a = 0.0, b = 0.0, c = 0.0, d = 0.0, e = 0.0;
    int32_t flag = -1;
    tqs_param_geometry_t geom{};
    const double* p = sym_necked.data();

    check(tqs_param_cache_radius_grid(nullptr, p, 5, radii.data(), n_theta) == 2, "NULL handle radius_grid -> 2");
    check(tqs_param_cache_radius_and_derivative(nullptr, p, 5, radii.data(), dr.data(), n_theta) == 2,
          "NULL handle radius_and_derivative -> 2");
    check(tqs_param_cache_radius_and_derivative_at_thetas(nullptr, p, 5, radii.data(), n_theta, radii.data(),
                                                          dr.data()) == 2, "NULL handle at_thetas -> 2");
    check(tqs_param_cache_shape(nullptr, p, 5, &a, &b) == 2, "NULL handle shape -> 2");
    check(tqs_param_cache_rho_z_grid(nullptr, p, 5, n_z, z.data(), rho.data(), drho.data()) == 2,
          "NULL handle rho_z_grid -> 2");
    check(tqs_param_cache_rho_z_grid_unchecked(nullptr, p, 5, n_z, z.data(), rho.data(), drho.data()) == 2,
          "NULL handle rho_z_grid_unchecked -> 2");
    check(tqs_param_cache_neck(nullptr, p, 5, &a, &b, &flag) == 2, "NULL handle neck -> 2");
    check(tqs_param_cache_geometry(nullptr, p, 5, &geom) == 2, "NULL handle geometry -> 2");
    check(tqs_param_cache_moller_coordinates(nullptr, p, 5, &a, &b, &c, &d, &e) == 2, "NULL handle moller -> 2");
    check(tqs_param_cache_scission_sigma1(nullptr, p, 5, &a, &b, &flag) == 2, "NULL handle scission -> 2");

    const std::vector<double> six{1.6, -0.6, 1.0, 0.0, 0.0, 0.0};
    check(tqs_param_cache_radius_grid(cache, p, 4, radii.data(), n_theta) == 4, "n_params 4 -> 4");
    check(tqs_param_cache_radius_grid(cache, six.data(), 6, radii.data(), n_theta) == 4, "n_params 6 -> 4");
    check(tqs_param_cache_radius_grid(cache, p, -1, radii.data(), n_theta) == 4, "negative n_params -> 4");
    check(tqs_param_cache_radius_grid(cache, nullptr, 0, radii.data(), n_theta) == 4,
          "NULL params with n_params 0 -> 4");
    check(tqs_param_cache_radius_grid(cache, p, 5, radii.data(), n_theta - 1) == 3, "n_radii 63 -> 3");
    check(tqs_param_cache_radius_grid(cache, p, 5, radii.data(), -1) == 3, "negative n_radii -> 3");
    check(tqs_param_cache_rho_z_grid(cache, p, 5, -5, z.data(), rho.data(), drho.data()) == 3,
          "negative n_z -> 3");
    check(tqs_param_cache_rho_z_grid(cache, p, 5, 2, z.data(), rho.data(), drho.data()) == 3,
          "n_z below the floor -> 3");
    check(tqs_param_cache_radius_and_derivative_at_thetas(cache, p, 5, nullptr, -2, nullptr, nullptr) == 3,
          "at_thetas negative n_thetas -> 3");
    check(tqs_param_cache_radius_grid(cache, p, 4, radii.data(), -1) == 4, "4 outranks 3");
}

void test_valid_shape(const tqs_param_cache_t* cache, const std::vector<double>& thetas) {
    std::vector<double> radii(n_theta), radii_d(n_theta), dr(n_theta), at_r(n_theta), at_dr(n_theta);
    std::vector<double> z(n_z), rho(n_z), drho(n_z);
    double rn = 0.0, rs = 0.0, zn = 0.0, rhon = 0.0, q2 = 0.0, d = 0.0, ag = 0.0, e1 = 0.0, e2 = 0.0;
    double s1t = 0.0, s1s = 0.0;
    int32_t found = -1, defined = -1;
    tqs_param_geometry_t geom{};
    const double* p = sym_necked.data();

    check(tqs_param_cache_radius_grid(cache, p, 5, radii.data(), n_theta) == 0, "SYM_NECKED radius_grid 0");
    check(all_finite_positive(radii), "SYM_NECKED radii finite positive");
    check(tqs_param_cache_radius_and_derivative(cache, p, 5, radii_d.data(), dr.data(), n_theta) == 0,
          "SYM_NECKED radius_and_derivative 0");
    check(tqs_param_cache_radius_and_derivative_at_thetas(cache, p, 5, thetas.data(), n_theta, at_r.data(),
                                                          at_dr.data()) == 0, "SYM_NECKED at_thetas 0");
    check(std::memcmp(at_r.data(), radii_d.data(), sizeof(double) * n_theta) == 0,
          "SYM_NECKED at_thetas == radius_and_derivative bytes");
    check(tqs_param_cache_shape(cache, p, 5, &rn, &rs) == 0 && rn > 0.0 && rs > 0.0, "SYM_NECKED shape 0");
    check(tqs_param_cache_rho_z_grid(cache, p, 5, n_z, z.data(), rho.data(), drho.data()) == 0,
          "SYM_NECKED rho_z_grid 0");
    check(tqs_param_cache_rho_z_grid_unchecked(cache, p, 5, n_z, z.data(), rho.data(), drho.data()) == 0,
          "SYM_NECKED rho_z_grid_unchecked 0");
    check(tqs_param_cache_neck(cache, p, 5, &zn, &rhon, &found) == 0 && found == 1 && rhon > 0.0,
          "SYM_NECKED neck found with rho_neck > 0");
    check(tqs_param_cache_geometry(cache, p, 5, &geom) == 0 && geom.mid_kind == TQS_MID_ONE_SHEET
          && geom.waist_inside == 1 && geom.is_separated == 0 && geom.is_pinched == 0,
          "SYM_NECKED geometry one_sheet, waist inside");
    check(geom.r_north == rn && geom.r_south == rs, "SYM_NECKED geometry poles == shape poles");
    check(tqs_param_cache_moller_coordinates(cache, p, 5, &q2, &d, &ag, &e1, &e2) == 0 && q2 > 0.0 && d > 0.0,
          "SYM_NECKED moller 0");
    check(tqs_param_cache_scission_sigma1(cache, p, 5, &s1t, &s1s, &defined) == 0 && defined == 1 && s1s > s1t,
          "SYM_NECKED scission defined");
}

void test_separated(const tqs_param_cache_t* cache) {
    std::vector<double> radii(n_theta, 7.0), z(n_z, 7.0), rho(n_z, 7.0), drho(n_z, 7.0);
    double rn = 7.0, rs = 7.0;
    tqs_param_geometry_t geom{};
    const double* p = separated.data();

    check(tqs_param_cache_radius_grid(cache, p, 5, radii.data(), n_theta) == TQS_ERROR_SEPARATED,
          "SEPARATED radius_grid -> 103");
    check(all_zero(radii), "SEPARATED radii zero-filled");
    check(tqs_param_cache_rho_z_grid(cache, p, 5, n_z, z.data(), rho.data(), drho.data()) == TQS_ERROR_SEPARATED,
          "SEPARATED rho_z_grid -> 103");
    check(all_zero(z) && all_zero(rho) && all_zero(drho), "SEPARATED rho_z_grid zero-filled");
    check(tqs_param_cache_rho_z_grid_unchecked(cache, p, 5, n_z, z.data(), rho.data(), drho.data()) == 0,
          "SEPARATED unchecked -> 0");
    check(rho[n_z / 2] == 0.0 && rho[1] > 0.0 && rho[n_z - 2] > 0.0, "SEPARATED unchecked draws two fragments");
    check(tqs_param_cache_geometry(cache, p, 5, &geom) == 0 && geom.is_separated == 1 && geom.r_north > 0.0,
          "SEPARATED geometry is_separated");
    check(tqs_param_cache_shape(cache, p, 5, &rn, &rs) == TQS_ERROR_SEPARATED && rn == 0.0 && rs == 0.0,
          "SEPARATED shape -> 103, poles zero-filled");
}

void test_one_shot(const std::vector<double>& thetas) {
    std::vector<double> radii(n_theta), dr(n_theta), z(n_z), rho(n_z), drho(n_z);
    double a = 0.0, b = 0.0, c = 0.0, d = 0.0, e = 0.0;
    int32_t flag = -1, status = -1;
    tqs_param_geometry_t geom{};
    const double* p = sym_necked.data();
    const std::vector<double> six{1.6, -0.6, 1.0, 0.0, 0.0, 0.0};

    // A NULL status pointer is allowed on every entry.
    tqs_param_radius_grid(p, 5, thetas.data(), n_theta, radii.data(), nullptr);
    check(all_finite_positive(radii), "one-shot radius_grid with NULL status computes");
    tqs_param_radius_and_derivative(p, 5, thetas.data(), n_theta, radii.data(), dr.data(), nullptr);
    tqs_param_shape(p, 5, &a, &b, nullptr);
    check(a > 0.0 && b > 0.0, "one-shot shape with NULL status computes");
    tqs_param_rho_z_grid(p, 5, n_z, z.data(), rho.data(), drho.data(), nullptr);
    tqs_param_rho_z_grid_unchecked(p, 5, n_z, z.data(), rho.data(), drho.data(), nullptr);
    tqs_param_neck(p, 5, &a, &b, &flag, nullptr);
    check(flag == 1, "one-shot neck with NULL status computes");
    tqs_param_geometry(p, 5, &geom, nullptr);
    tqs_param_moller_coordinates(p, 5, &a, &b, &c, &d, &e, nullptr);
    tqs_param_scission_sigma1(p, 5, &a, &b, &flag, nullptr);

    tqs_param_radius_grid(p, 5, thetas.data(), n_theta, radii.data(), &status);
    check(status == 0, "one-shot radius_grid 0");
    tqs_param_radius_grid(separated.data(), 5, thetas.data(), n_theta, radii.data(), &status);
    check(status == TQS_ERROR_SEPARATED && all_zero(radii), "one-shot SEPARATED radius_grid -> 103");
    tqs_param_rho_z_grid(separated.data(), 5, n_z, z.data(), rho.data(), drho.data(), &status);
    check(status == TQS_ERROR_SEPARATED, "one-shot SEPARATED rho_z_grid -> 103");
    tqs_param_rho_z_grid_unchecked(separated.data(), 5, n_z, z.data(), rho.data(), drho.data(), &status);
    check(status == 0, "one-shot SEPARATED unchecked -> 0");
    tqs_param_geometry(separated.data(), 5, &geom, &status);
    check(status == 0 && geom.is_separated == 1, "one-shot SEPARATED geometry");

    tqs_param_radius_grid(p, -1, thetas.data(), n_theta, radii.data(), &status);
    check(status == TQS_ERROR_INVALID_INIT, "one-shot negative n_params -> 5");
    tqs_param_radius_grid(p, 5, thetas.data(), -1, nullptr, &status);
    check(status == TQS_ERROR_INVALID_INIT, "one-shot negative n_thetas -> 5");
    tqs_param_rho_z_grid(p, 5, -1, nullptr, nullptr, nullptr, &status);
    check(status == TQS_ERROR_INVALID_INIT, "one-shot negative n_z -> 5");
    tqs_param_shape(six.data(), 6, &a, &b, &status);
    check(status == TQS_ERROR_TOO_MANY_PARAMS, "one-shot n_params 6 -> 1");
    tqs_param_shape(nullptr, 0, &a, &b, &status);
    check(status == TQS_ERROR_WRONG_PARAM_COUNT, "one-shot n_params 0 -> 4");
    tqs_param_scission_sigma1(p, 4, &a, &b, &flag, &status);
    check(status == TQS_ERROR_WRONG_PARAM_COUNT, "one-shot scission n_params 4 -> 4");
}

void test_messages() {
    const int32_t codes[] = {0, 3, 100, 101, 102, 103, 999};
    for (const int32_t code : codes) {
        const char* m = tqs_param_status_message(code);
        check(m != nullptr && std::strlen(m) > 0, "status message non-empty");
    }
    check(tqs_param_status_message(101) == tqs_param_status_message(101), "message pointer stable");
    check(std::strcmp(tqs_param_status_message(103), tqs_param_status_message(102)) != 0,
          "distinct codes, distinct messages");
}
}  // namespace

int main() {
    const std::vector<double> thetas = theta_grid(n_theta);
    test_create(thetas);

    int32_t status = -1;
    tqs_param_cache_t* cache = tqs_param_cache_create(5, thetas.data(), n_theta, &status);
    check(cache != nullptr && status == 0, "create max_params 5 with 64 thetas -> handle, 0");
    if (cache == nullptr) {
        std::printf("c_api_smoke_test: no cache, %d failure(s)\n", failures);
        return 1;
    }

    test_usage_codes(cache);
    test_valid_shape(cache, thetas);
    test_separated(cache);
    test_one_shot(thetas);
    test_messages();

    tqs_param_cache_destroy(cache);
    std::printf("c_api_smoke_test: %d failure(s)\n", failures);
    return failures == 0 ? 0 : 1;
}
