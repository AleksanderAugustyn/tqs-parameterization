// Contract family 3: concurrency. ONE read-only cache, shared by every thread.
//
// A cache is immutable after creation, so any number of threads may compute on
// one handle at once. The main thread first walks every shape of every thread
// serially through the shared handle and records all ten cached outputs. Then
// 8 threads repeat exactly their own walk concurrently on the SAME handle.
// Concurrency is the only difference, and it may not move a bit: every status
// must match and every output must be identical byte for byte.
#include "tqs_parameterization.h"

#include <cstdio>
#include <cstring>
#include <numbers>
#include <thread>
#include <vector>

namespace {
constexpr int n_theta = 64;
constexpr int n_other = 9;
constexpr int n_z = 129;
constexpr int n_threads = 8;
constexpr int n_steps = 100;
constexpr int n_out = 10;

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

// Deterministic per-(thread, step) shape inside the valid box: necked
// one-sheet middles with mildly asymmetric ends.
std::vector<double> shape_for(const int t, const int s) {
    return {1.2 + 0.012 * s + 0.01 * t, -0.6 - 0.002 * s, 1.0, 0.05 * t - 0.2, 0.1 - 0.004 * s};
}

struct StepResult {
    int32_t status[n_out] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
    std::vector<double> radii, radii_d, dr, at_radii, at_dr;
    std::vector<double> z, rho, drho, raw_z, raw_rho, raw_drho;
    double scalars[11] = {};
    int32_t found = -1, defined = -1;
    tqs_param_geometry_t geom{};
};

using Walk = std::vector<StepResult>;

// The one code path both the serial reference and the threads use.
void run_walk(const tqs_param_cache_t* cache, const std::vector<double>& other_thetas,
              const int thread_id, Walk& out) {
    out.assign(static_cast<std::size_t>(n_steps), StepResult{});
    for (int step = 0; step < n_steps; ++step) {
        StepResult& r = out[static_cast<std::size_t>(step)];
        const std::vector<double> p = shape_for(thread_id, step);
        r.radii.assign(n_theta, 0.0);
        r.radii_d.assign(n_theta, 0.0);
        r.dr.assign(n_theta, 0.0);
        r.at_radii.assign(n_other, 0.0);
        r.at_dr.assign(n_other, 0.0);
        r.z.assign(n_z, 0.0);
        r.rho.assign(n_z, 0.0);
        r.drho.assign(n_z, 0.0);
        r.raw_z.assign(n_z, 0.0);
        r.raw_rho.assign(n_z, 0.0);
        r.raw_drho.assign(n_z, 0.0);

        r.status[0] = tqs_param_cache_radius_grid(cache, p.data(), 5, r.radii.data(), n_theta);
        r.status[1] = tqs_param_cache_radius_and_derivative(cache, p.data(), 5, r.radii_d.data(),
                                                            r.dr.data(), n_theta);
        r.status[2] = tqs_param_cache_radius_and_derivative_at_thetas(
                cache, p.data(), 5, other_thetas.data(), n_other, r.at_radii.data(), r.at_dr.data());
        r.status[3] = tqs_param_cache_shape(cache, p.data(), 5, &r.scalars[0], &r.scalars[1]);
        r.status[4] = tqs_param_cache_rho_z_grid(cache, p.data(), 5, n_z, r.z.data(), r.rho.data(),
                                                 r.drho.data());
        r.status[5] = tqs_param_cache_rho_z_grid_unchecked(cache, p.data(), 5, n_z, r.raw_z.data(),
                                                           r.raw_rho.data(), r.raw_drho.data());
        r.status[6] = tqs_param_cache_neck(cache, p.data(), 5, &r.scalars[2], &r.scalars[3], &r.found);
        r.status[7] = tqs_param_cache_geometry(cache, p.data(), 5, &r.geom);
        r.status[8] = tqs_param_cache_moller_coordinates(cache, p.data(), 5, &r.scalars[4], &r.scalars[5],
                                                         &r.scalars[6], &r.scalars[7], &r.scalars[8]);
        r.status[9] = tqs_param_cache_scission_sigma1(cache, p.data(), 5, &r.scalars[9], &r.scalars[10],
                                                      &r.defined);
    }
}

bool same_bytes(const std::vector<double>& a, const std::vector<double>& b) {
    return a.size() == b.size() && std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0;
}

bool same_step(const StepResult& a, const StepResult& b) {
    return std::memcmp(a.status, b.status, sizeof(a.status)) == 0
           && std::memcmp(a.scalars, b.scalars, sizeof(a.scalars)) == 0
           && std::memcmp(&a.geom, &b.geom, sizeof(a.geom)) == 0
           && a.found == b.found && a.defined == b.defined
           && same_bytes(a.radii, b.radii) && same_bytes(a.radii_d, b.radii_d) && same_bytes(a.dr, b.dr)
           && same_bytes(a.at_radii, b.at_radii) && same_bytes(a.at_dr, b.at_dr)
           && same_bytes(a.z, b.z) && same_bytes(a.rho, b.rho) && same_bytes(a.drho, b.drho)
           && same_bytes(a.raw_z, b.raw_z) && same_bytes(a.raw_rho, b.raw_rho)
           && same_bytes(a.raw_drho, b.raw_drho);
}

bool all_valid(const StepResult& r) {
    for (const int32_t s : r.status) {
        if (s != TQS_VALID) return false;
    }
    return true;
}
}  // namespace

int main() {
    const std::vector<double> thetas = theta_grid(n_theta);
    const std::vector<double> other_thetas = theta_grid(n_other);

    int32_t status = -1;
    tqs_param_cache_t* cache = tqs_param_cache_create(5, thetas.data(), n_theta, &status);
    if (cache == nullptr) {
        std::printf("FAIL: cache_create returned NULL (status %d)\n", status);
        return 1;
    }

    std::vector<Walk> serial(static_cast<std::size_t>(n_threads));
    for (int t = 0; t < n_threads; ++t) {
        run_walk(cache, other_thetas, t, serial[static_cast<std::size_t>(t)]);
    }

    std::vector<Walk> threaded(static_cast<std::size_t>(n_threads));
    {
        std::vector<std::jthread> pool;
        pool.reserve(static_cast<std::size_t>(n_threads));
        for (int t = 0; t < n_threads; ++t) {
            pool.emplace_back([cache, &other_thetas, &threaded, t] {
                run_walk(cache, other_thetas, t, threaded[static_cast<std::size_t>(t)]);
            });
        }
    }  // all joined here

    int failures = 0;
    int valid_steps = 0;
    for (int t = 0; t < n_threads; ++t) {
        for (int step = 0; step < n_steps; ++step) {
            const StepResult& a = threaded[static_cast<std::size_t>(t)][static_cast<std::size_t>(step)];
            const StepResult& b = serial[static_cast<std::size_t>(t)][static_cast<std::size_t>(step)];
            if (!same_step(a, b)) {
                std::printf("FAIL: thread %d step %d differs from the serial reference\n", t, step);
                ++failures;
            }
            if (all_valid(b)) ++valid_steps;
        }
    }

    // A rejected shape zero-fills every buffer and two zero buffers always
    // match: the comparison must not pass vacuously.
    if (valid_steps == 0) {
        std::printf("FAIL: no step was fully VALID\n");
        ++failures;
    }

    tqs_param_cache_destroy(cache);
    std::printf("thread_stress_test: %d thread(s) x %d steps on one shared cache, %d of %d steps fully "
                "valid, %d failure(s)\n", n_threads, n_steps, valid_steps, n_threads * n_steps, failures);
    return failures == 0 ? 0 : 1;
}
