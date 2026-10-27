/**
 * @file tqs_parameterization.h
 * @brief C API for the Fortran three-quadratic-surface (3QS) shape library (v0.1.0).
 *
 * The shape is Nix's (Nucl. Phys. A 130 (1969) 241): two end spheroids joined
 * smoothly by a middle body that is a spheroid, cylinder, hyperboloid of one
 * sheet, cone or hyperboloid of two sheets. params = [sigma1, sigma2, sigma3,
 * alpha2, alpha3], Nix's coordinates, exactly TQS_PARAM_N_PARAMS (5) of them.
 * The sphere is [0, 1, 1, 0, 0]. All lengths are in R0 = 1 units with the
 * centre of mass at the origin.
 *
 * Two tiers of computation:
 *   - One-shot (`tqs_param_<compute>`): one call computes a shape and keeps
 *     nothing.
 *   - Read-only cache (`tqs_param_cache_*`): build a cache once with
 *     tqs_param_cache_create(), then request any number of shapes against it.
 *     The cache holds only the primary theta set; nothing derived from a
 *     shape parameter is stored, so no call influences a later one.
 *
 *   The two tiers return the same status and bitwise-identical outputs.
 *
 * Parameters:
 *   tqs is fixed-arity. max_params at creation must be 5 (larger -> 1,
 *   smaller -> 5). A cached call needs n_params == 5, anything else is 4. A
 *   one-shot call with n_params > 5 is 1 and with fewer (0 included) is 4.
 *   Nothing is padded or trimmed.
 *
 * Thread safety:
 *   A cache is IMMUTABLE after creation. Every compute takes it `const`, and
 *   any number of threads may compute on one handle concurrently.
 *   tqs_param_cache_create() and tqs_param_cache_destroy() on a handle must
 *   not race with any other call on it. The one-shot calls hold no state.
 *
 * Diagnostics:
 *   tqs_param_cache_create() returns NULL on failure and reports the
 *   rejecting code through its nullable `status`. The cached computes return
 *   the status directly; the one-shot calls report through a trailing
 *   nullable `int32_t* status`. TQS_VALID (0) means success.
 *   tqs_param_status_message() maps a code to a static, null-terminated
 *   string owned by the library and safe to read from any thread.
 *
 * Order of the checks (the first failing one sets the status):
 *   create        max_params > 5 -> 1; < 5 -> 5; then the theta set -> 3
 *                 (none, or one outside [0, pi])
 *   cached call   NULL handle -> 2; n_params != 5 -> 4; output size -> 3;
 *                 at-thetas only: the theta set -> 3; then the value codes
 *   one-shot      a negative size argument -> 5; n_params > 5 -> 1; < 5 -> 4;
 *                 output size and theta set -> 3; then the value codes
 *   Value codes, in order: 100 (coordinates), 101 (no geometry), then 103
 *   (separated) and 102 (pinched or not star-convex) where an output gates
 *   them: the R(theta) outputs and shape apply both, the checked rho(z) grid
 *   applies 103, the rest neither; scission_sigma1 fails only with 100.
 *
 * Failure behavior:
 *   On any nonzero status every output is zero-filled over its stated
 *   extent, `found` and `defined` are 0, and a geometry record is all zero
 *   (mid_kind TQS_MID_NONE).
 *
 * Buffer sizes:
 *   Every compute takes the size of its output buffers explicitly: n_radii
 *   must equal the handle's n_thetas, the at-thetas outputs are n_thetas long,
 *   and a rho(z) grid is n_z >= TQS_PARAM_N_POINTS_FLOOR nodes. A wrong size
 *   is 3 with nothing written past the stated size. A stated size MUST be
 *   the actual extent of the buffer passed: the library zero-fills exactly
 *   that many elements. A negative output size in a cached call is 3, a
 *   negative n_params there is an empty vector (4); in a one-shot call a
 *   negative n_params, n_thetas or n_z is a bad argument (5). A NULL pointer
 *   is safe wherever its stated size is 0 or negative.
 *
 * Precondition — finite input:
 *   `params` and `thetas` must be finite. Release builds assume no NaN, so a
 *   non-finite input is undefined behavior.
 */

#ifndef TQS_PARAMETERIZATION_H
#define TQS_PARAMETERIZATION_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* --- Limits --- */
/** The fixed parameter count: [sigma1, sigma2, sigma3, alpha2, alpha3]. */
#define TQS_PARAM_N_PARAMS 5
/** Fewest nodes of a rho(z) grid: the two poles and one interior node. */
#define TQS_PARAM_N_POINTS_FLOOR 3

/* --- Shared contract status codes (0-99, identical in every
 *     shape-parameterization library) --- */
#define TQS_VALID                          0
#define TQS_ERROR_TOO_MANY_PARAMS          1
#define TQS_ERROR_CACHE_NOT_INITIALIZED    2
#define TQS_ERROR_INVALID_GRID             3
#define TQS_ERROR_WRONG_PARAM_COUNT        4
#define TQS_ERROR_INVALID_INIT             5

/* --- tqs-specific status codes (>= 100) --- */
#define TQS_ERROR_INVALID_COORDINATES      100
#define TQS_ERROR_NO_GEOMETRY              101
#define TQS_ERROR_NOT_STAR_CONVEX          102
#define TQS_ERROR_SEPARATED                103

/* --- Middle-body kinds (tqs_param_geometry_t.mid_kind) --- */
#define TQS_MID_NONE                       0
#define TQS_MID_SPHEROID                   1
#define TQS_MID_CYLINDER                   2
#define TQS_MID_ONE_SHEET                  3
#define TQS_MID_CONE                       4
#define TQS_MID_TWO_SHEET                  5

/** Opaque read-only cache. Created once, shared across threads. */
typedef struct tqs_param_cache tqs_param_cache_t;

/**
 * The resolved geometry. Segment 0 (C index) is the south end spheroid, 1 the
 * north one, 2 the middle body: rho^2 = a_sq[i] - k[i] (z - l[i])^2. c[0..1]
 * are the end semi-axes along z, z_join[0..1] the joins, r_north and r_south
 * the pole distances from the centre of mass, alpha1 Nix's centre-of-mass
 * coordinate. The flags are 0 or 1.
 */
typedef struct {
    double u;
    double a_sq[3];
    double k[3];
    double c[2];
    double l[3];
    double z_join[2];
    double r_north;
    double r_south;
    double alpha1;
    int32_t mid_kind;
    int32_t waist_inside;
    int32_t is_separated;
    int32_t is_pinched;
} tqs_param_geometry_t;

/* ===================================================================== */
/* Diagnostics                                                            */
/* ===================================================================== */

/** Static, null-terminated description of a status code. Never NULL. */
const char *tqs_param_status_message(int32_t status);

/* ===================================================================== */
/* Cache lifecycle                                                        */
/* ===================================================================== */

/**
 * Build a read-only cache on a primary theta set.
 *
 * @param max_params  TQS_PARAM_N_PARAMS (5)
 * @param thetas      n_thetas polar angles in [0, pi], any order
 * @param n_thetas    At least 1
 * @param status      Receives 0 or the rejecting code (1, 5 or 3). May be NULL.
 * @return            Handle, or NULL on failure
 */
tqs_param_cache_t *tqs_param_cache_create(int32_t max_params, const double *thetas,
                                          int32_t n_thetas, int32_t *status);

/** Release a cache. NULL-safe. Must not race with any other call on it. */
void tqs_param_cache_destroy(tqs_param_cache_t *cache);

/* ===================================================================== */
/* Cached computes — const handle, safe to call concurrently              */
/* ===================================================================== */

/** R(theta) at the handle's own thetas. n_radii must equal n_thetas. */
int32_t tqs_param_cache_radius_grid(const tqs_param_cache_t *cache, const double *params,
                                    int32_t n_params, double *radii, int32_t n_radii);

/** R and dR/dtheta at the handle's own thetas. n_radii must equal n_thetas. */
int32_t tqs_param_cache_radius_and_derivative(const tqs_param_cache_t *cache,
                                              const double *params, int32_t n_params,
                                              double *radii, double *dr_dtheta,
                                              int32_t n_radii);

/** R and dR/dtheta at caller thetas (each in [0, pi]); outputs n_thetas long. */
int32_t tqs_param_cache_radius_and_derivative_at_thetas(const tqs_param_cache_t *cache,
                                                        const double *params,
                                                        int32_t n_params,
                                                        const double *thetas,
                                                        int32_t n_thetas, double *radii,
                                                        double *dr_dtheta);

/** Pole radii with the R(theta) verdict (102 / 103 zero-fill them). */
int32_t tqs_param_cache_shape(const tqs_param_cache_t *cache, const double *params,
                              int32_t n_params, double *r_north, double *r_south);

/** rho(z), drho/dz on n_z uniform nodes over [z_south, z_north]; 103 if separated. */
int32_t tqs_param_cache_rho_z_grid(const tqs_param_cache_t *cache, const double *params,
                                   int32_t n_params, int32_t n_z, double *z, double *rho,
                                   double *drho_dz);

/** As tqs_param_cache_rho_z_grid, never 103: rho = 0 inside the void. */
int32_t tqs_param_cache_rho_z_grid_unchecked(const tqs_param_cache_t *cache,
                                             const double *params, int32_t n_params,
                                             int32_t n_z, double *z, double *rho,
                                             double *drho_dz);

/** Neck: found = 1 for a waist inside a sigma2 < 0 middle; rho_neck = 0 when split. */
int32_t tqs_param_cache_neck(const tqs_param_cache_t *cache, const double *params,
                             int32_t n_params, double *z_neck, double *rho_neck,
                             int32_t *found);

/** The resolved geometry; valid for separated and pinched shapes too. */
int32_t tqs_param_cache_geometry(const tqs_param_cache_t *cache, const double *params,
                                 int32_t n_params, tqs_param_geometry_t *geometry);

/** Moller's grid coordinates as diagnostics: q2 = <r^2 P2> (R0 = 1), d, alpha_g, eps. */
int32_t tqs_param_cache_moller_coordinates(const tqs_param_cache_t *cache,
                                           const double *params, int32_t n_params,
                                           double *q2, double *d, double *alpha_g,
                                           double *eps_f1, double *eps_f2);

/** Scission map: sigma1_t always; sigma1_star with defined = 1 when sigma2 < 0. */
int32_t tqs_param_cache_scission_sigma1(const tqs_param_cache_t *cache,
                                        const double *params, int32_t n_params,
                                        double *sigma1_t, double *sigma1_star,
                                        int32_t *defined);

/* ===================================================================== */
/* One-shot computes — status is a trailing nullable out-parameter        */
/* ===================================================================== */

void tqs_param_radius_grid(const double *params, int32_t n_params, const double *thetas,
                           int32_t n_thetas, double *radii, int32_t *status);

void tqs_param_radius_and_derivative(const double *params, int32_t n_params,
                                     const double *thetas, int32_t n_thetas,
                                     double *radii, double *dr_dtheta, int32_t *status);

void tqs_param_shape(const double *params, int32_t n_params, double *r_north,
                     double *r_south, int32_t *status);

void tqs_param_rho_z_grid(const double *params, int32_t n_params, int32_t n_z, double *z,
                          double *rho, double *drho_dz, int32_t *status);

void tqs_param_rho_z_grid_unchecked(const double *params, int32_t n_params, int32_t n_z,
                                    double *z, double *rho, double *drho_dz,
                                    int32_t *status);

void tqs_param_neck(const double *params, int32_t n_params, double *z_neck,
                    double *rho_neck, int32_t *found, int32_t *status);

void tqs_param_geometry(const double *params, int32_t n_params,
                        tqs_param_geometry_t *geometry, int32_t *status);

void tqs_param_moller_coordinates(const double *params, int32_t n_params, double *q2,
                                  double *d, double *alpha_g, double *eps_f1,
                                  double *eps_f2, int32_t *status);

void tqs_param_scission_sigma1(const double *params, int32_t n_params, double *sigma1_t,
                               double *sigma1_star, int32_t *defined, int32_t *status);

#ifdef __cplusplus
}
#endif

#endif /* TQS_PARAMETERIZATION_H */
