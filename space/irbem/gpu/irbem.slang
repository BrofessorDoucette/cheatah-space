// Copyright (c) 2026 BigBrain LLC. MIT-licensed (see LICENSE).
//
// irbem.slang — space.irbem's OWN device kernels, written once in Slang and compiled to SPIR-V by
// slangc at build time, exactly the way cheatah-gpu-linalg compiles kernels/linalg.slang. This
// file is deliberately separate from linalg.slang: linear algebra is a library, magnetospheric
// physics is not, and a field model does not belong in a BLAS.
//
// ABI — the SAME contract cheatah-gpu-linalg's contexts implement, because the dispatch seam in
// dispatch.hpp drives those contexts and nothing else:
//   - buffers are bound at set 0, bindings 0..n-1, in the order dispatch.hpp passes them;
//   - the LAST buffer of every entry is a `uint` dims buffer, spelled out per kernel below;
//   - 1-D entries use [numthreads(256,1,1)] — kernels.hpp's kLocal1d, the workgroup size
//     Context::dispatch_1d rounds the grid up to — so every entry MUST bounds-check against dims;
//   - vectors are stored xyz-interleaved and contiguous, matching the host's std::span<float>.
//
// fp32 is not a shortcut here, it is the lane: ERROR_BUDGET.md §3 measures the discretization floor
// at 1.2e-3 … 4.0e-3 relative while one fp32 operation costs 6e-8, and SPIR-V's GLSL.std.450
// transcendentals exist only for 16- and 32-bit floats anyway. The ACCUMULATION is what fp32 would
// break, and no kernel in this file accumulates: one thread, one point, no reduction.
//
// Every entry point is registered in dispatch.hpp's `registered_kernels`, and a test parses this
// file to assert the two lists agree — a kernel added here without a registry entry, or a registry
// entry with no kernel, fails the suite instead of failing at a customer's first dispatch.

// ---------------------------------------------------------------------------------------------
// irbem_dipole_f32 — the centred-dipole magnetic field, one thread per point.
//
// PHYSICS. The geomagnetic scalar potential is the IAGA spherical-harmonic series that defines
// IGRF (Alken et al., "International Geomagnetic Reference Field: the thirteenth generation",
// Earth Planets Space 73:49 (2021), eq. 1):
//
//     V(r, theta, phi) = a * sum_n (a/r)^(n+1) * sum_m [ g_n^m cos(m phi) + h_n^m sin(m phi) ]
//                                                      * P_n^m(cos theta),        B = -grad V
//
// with `a` the IGRF reference radius, theta the geocentric COLATITUDE and phi east longitude.
// Keeping only the degree-1 order-0 term — the centred dipole, which is what the MAG frame is
// built around — leaves V = a^3 g10 cos(theta) / r^2, and hence
//
//     B_r     = -dV/dr        =  2 a^3 g10 cos(theta) / r^3
//     B_theta = -(1/r) dV/dth =    a^3 g10 sin(theta) / r^3
//     B_phi   = 0
//
// In the MAG frame the dipole axis IS z-hat, so with r-hat = r/|r| and cos(theta) = z/r those two
// components are the single vector identity
//
//     B = g10 (a/r)^3 * ( 3 (z-hat . r-hat) r-hat - z-hat )
//
// (radial part 3cos(theta) - cos(theta) = 2cos(theta); theta part -(z-hat . theta-hat) =
// +sin(theta) — the two lines above, recovered). Taking a = 1 Re and r in Re:
//
//     Bx = 3 g10 x z / r^5,   By = 3 g10 y z / r^5,   Bz = g10 (3 z^2 - r^2) / r^5.
//
// g10 is NEGATIVE in the real field (-29404.8 nT for IGRF-13 epoch 2020.0), which is what makes B
// point INTO the ground at the northern dipole pole; the sign is carried by the coefficient, not by
// the formula, so this kernel never needs to know it. Units: positions in Earth radii, g10 and the
// output in nanotesla — the potential is linear in the coefficient, so B comes out in whatever unit
// g10 is given in.
//
// The origin is a genuine pole of the field, not a numerical accident. r^2 == 0 writes an exact
// zero rather than an inf/NaN, and dispatch.hpp's host lane makes the SAME choice, so the branch is
// a defined, testable value on both sides instead of a difference that only shows up in a NaN.
//
// bindings: 0 pos (3N floats, Re) | 1 out (3N floats, nT) | 2 params {g10} | 3 dims {N}
// ---------------------------------------------------------------------------------------------
// Guarded so each entry point compiles alone: every kernel here binds set 0 from slot 0, which is
// the ABI the contexts implement, so two kernels visible at once would overlap. linalg.slang uses
// the same -DKERNEL_<name> convention.
#ifdef KERNEL_DIPOLE
[[vk::binding(0, 0)]] StructuredBuffer<float>   dipole_pos;
[[vk::binding(1, 0)]] RWStructuredBuffer<float> dipole_out;
[[vk::binding(2, 0)]] StructuredBuffer<float>   dipole_params;
[[vk::binding(3, 0)]] StructuredBuffer<uint>    dipole_dims;

[shader("compute")]
[numthreads(256, 1, 1)]
void irbem_dipole_f32(uint3 tid : SV_DispatchThreadID) {
    const uint n = dipole_dims[0];
    const uint i = tid.x;
    if (i >= n) return;   // the grid is rounded up to a multiple of 256; the tail does nothing

    const float x = dipole_pos[3 * i + 0];
    const float y = dipole_pos[3 * i + 1];
    const float z = dipole_pos[3 * i + 2];
    const float g10 = dipole_params[0];

    const float r2 = x * x + y * y + z * z;
    float bx = 0.0f;
    float by = 0.0f;
    float bz = 0.0f;
    if (r2 > 0.0f) {
        // r^5 built from r^2 so there is exactly one sqrt and no pow: r^5 = (r^2)^2 * r. The host
        // reference in dispatch.hpp evaluates the SAME expressions in the SAME order, which is why
        // the two agree exactly on integer-valued geometry and to the fp32 budget elsewhere.
        const float r5 = r2 * r2 * sqrt(r2);
        bx = 3.0f * g10 * x * z / r5;
        by = 3.0f * g10 * y * z / r5;
        bz = g10 * (3.0f * z * z - r2) / r5;
    }

    dipole_out[3 * i + 0] = bx;
    dipole_out[3 * i + 1] = by;
    dipole_out[3 * i + 2] = bz;
}

#endif  // KERNEL_DIPOLE

// ---------------------------------------------------------------------------------------------
// irbem_igrf_f32 — the full IGRF-14 internal field, one thread per point.
//
// This is the kernel every other kernel in this file will call: a field-line trace evaluates it
// five times per RK4 step and ~150 steps per line, so ~10^5 times per L* point. It is also close
// to ideal GPU code — ONE square root, ZERO transcendentals, no data-dependent branch, and a few
// hundred bytes of coefficients that every lane reads identically.
//
// PHYSICS. The IAGA spherical-harmonic series (Alken et al., Earth Planets Space, 2025), evaluated
// in the singularity-free Cartesian form the host lane uses (space/irbem/igrf.hpp): the textbook
// formulation carries sin^m(theta) in P_n^m and then divides by sin(theta) for the phi component,
// which is 0/0 on the axis. Folding the sin^m into the longitude recursion makes every one of
// those cancel algebraically BEFORE any arithmetic happens, so the poles are ordinary points and
// the kernel needs no trig and no branch. That property is what makes it worth porting verbatim
// rather than re-deriving.
//
// TWO DEVIATIONS FROM THE HOST LANE, both deliberate and both about registers:
//
//  1. The host stores the whole triangular A_n^m and dA_n^m — 105 slots each at degree 13, 210
//     floats. That would sit at the edge of the 255-register-per-thread limit and spill. Degree n
//     of the recursion only reads degrees n-1 and n-2, so this keeps a ROLLING THREE-ROW WINDOW
//     (3 x 14 for A, 3 x 14 for dA = 84 floats) and folds the harmonic sum into the same pass.
//     Identical arithmetic, identical order, a third of the storage.
//  2. The coefficients arrive ALREADY INTERPOLATED to the epoch. The 26-epoch table is 5460
//     doubles and interpolating it per thread would be 10^5 redundant copies of a calculation the
//     host does once per batch. The device sees 105 g and 105 h floats and nothing else.
//
// The accumulation is fp32 here where the host uses fp64. That is sound and it is bounded: the sum
// is over 13 degrees, not 10^3 terms, so the error budget's accumulation concern (ERROR_BUDGET.md
// §3) does not bite -- the test measures the actual deviation against the fp64 host lane rather
// than assuming it.
//
// bindings: 0 pos (3N floats, Re, GEO) | 1 out (3N floats, nT, GEO)
//           2 coef (2*SLOTS floats: g[SLOTS] then h[SLOTS])
//           3 norm (2*SLOTS + NMAX+1 floats: e[SLOTS], f[SLOTS], diagonal[NMAX+1])
//           4 dims {N, NMAX}
// ---------------------------------------------------------------------------------------------

// Degree 13 is IGRF-14's full published degree; the arrays below are sized from it. A lower
// truncation is passed at runtime in dims[1] and simply ends the loops early.
// =============================================================================================
// The IGRF-14 internal field, as a SHARED function.
//
// Declared outside every kernel guard and taking its buffers as PARAMETERS, so both the batched
// field kernel and the field-line tracer call one implementation. The tracer evaluates it four
// times per RK4 step and ~120 steps per line — a second copy of this recursion would be a second
// place for the Legendre seeding to be got wrong, which has already happened once (seeding degree 1
// instead of degree 0 into the rolling window produced a 21% error).
//
// See the notes on irbem_igrf_f32 below for the physics and for why the recursion carries a rolling
// three-row window rather than the host's full 105-slot triangle.
// =============================================================================================

static const int kIgrfMaxDegree = 13;
static const int kIgrfMaxSlots = (kIgrfMaxDegree + 1) * (kIgrfMaxDegree + 2) / 2;  // 105
static const int kIgrfRow = kIgrfMaxDegree + 1;                                     // 14

/// Triangular packing, identical to the host's detail::slot_index: n(n+1)/2 + m.
int igrf_slot(int n, int m) { return ((n * (n + 1)) / 2) + m; }

/// The internal field at a GEO point, in nT.
/// @param coef  g[slots] then h[slots], already interpolated to the epoch on the host.
/// @param norm  e[slots], f[slots], diagonal[nmax+1] — the Legendre normalisation.
/// @param nmax  truncation degree. @param slots (nmax+1)(nmax+2)/2. @param p position in Re.
float3 igrf_eval(StructuredBuffer<float> coef, StructuredBuffer<float> norm,
                 int nmax, int slots, float3 p) {
    const float r2 = dot(p, p);
    // The origin is a pole of the field, not a numerical accident: an exact zero, matching the host.
    if (r2 == 0.0f) return float3(0.0f, 0.0f, 0.0f);

    const float r = sqrt(r2);
    const float inv_r = 1.0f / r;
    const float ux = p.x * inv_r, uy = p.y * inv_r, uz = p.z * inv_r;

    float cs[kIgrfRow + 1];
    float ss[kIgrfRow + 1];
    cs[0] = 0.0f; ss[0] = 0.0f; cs[1] = 1.0f; ss[1] = 0.0f;
    for (int m = 1; m <= nmax; ++m) {
        cs[m + 1] = (ux * cs[m]) - (uy * ss[m]);
        ss[m + 1] = (ux * ss[m]) + (uy * cs[m]);
    }

    // Rolling window; the INVARIANT at the top of iteration n is that a1 holds degree n-1 and a2
    // holds degree n-2, so the seed is degree 0 and degree 1 is produced by the loop like any other.
    float a2[kIgrfRow], a1[kIgrfRow], a0[kIgrfRow];
    float d2[kIgrfRow], d1[kIgrfRow], d0[kIgrfRow];
    for (int q = 0; q < kIgrfRow; ++q) {
        a2[q] = 0.0f; a1[q] = 0.0f; a0[q] = 0.0f;
        d2[q] = 0.0f; d1[q] = 0.0f; d0[q] = 0.0f;
    }
    a1[0] = 1.0f;

    float t0 = 0.0f, t1 = 0.0f, t2 = 0.0f, su = 0.0f, sv = 0.0f;
    float rn = inv_r * inv_r * inv_r;

    for (int n = 1; n <= nmax; ++n) {
        for (int q = 0; q < kIgrfRow; ++q) { a0[q] = 0.0f; d0[q] = 0.0f; }
        if (n == 1) {
            a0[0] = uz; d0[0] = 1.0f; a0[1] = norm[(2 * slots) + 1];
        } else {
            const int rbase = igrf_slot(n, 0);
            for (int m = 0; m < n; ++m) {
                const int k = rbase + m;
                const float e = norm[k];
                const float f = norm[slots + k];
                a0[m] = (e * uz * a1[m]) - (f * a2[m]);
                d0[m] = (e * (a1[m] + (uz * d1[m]))) - (f * d2[m]);
            }
            a0[n] = norm[(2 * slots) + n];
        }

        const int base = igrf_slot(n, 0);
        float p0 = 0.0f, p1 = 0.0f, p2 = 0.0f, pu = 0.0f, pv = 0.0f;
        for (int m = 0; m <= n; ++m) {
            const int k = base + m;
            const float gc = coef[k];
            const float hc = coef[slots + k];
            const float bigg = (gc * cs[m + 1]) + (hc * ss[m + 1]);
            const float ga = bigg * a0[m];
            const float fm = float(m);
            p0 += ga;
            p2 += fm * ga;
            p1 += bigg * d0[m];
            const float ma = fm * a0[m];
            pu += ma * ((gc * cs[m]) + (hc * ss[m]));
            pv += ma * ((hc * cs[m]) - (gc * ss[m]));
        }
        t0 += rn * float(n + 1) * p0;
        t1 += rn * p1;
        t2 += rn * p2;
        su += rn * pu;
        sv += rn * pv;
        rn *= inv_r;

        for (int q = 0; q < kIgrfRow; ++q) {
            a2[q] = a1[q]; a1[q] = a0[q];
            d2[q] = d1[q]; d1[q] = d0[q];
        }
    }

    const float w = uz;
    const float radial = t0 + (w * t1) + t2;
    return float3((ux * radial) - su,
                  (uy * radial) - sv,
                  (w * t0) - ((1.0f - (w * w)) * t1) + (w * t2));
}

#ifdef KERNEL_IGRF
[[vk::binding(0, 0)]] StructuredBuffer<float>   igrf_pos;
[[vk::binding(1, 0)]] RWStructuredBuffer<float> igrf_out;
[[vk::binding(2, 0)]] StructuredBuffer<float>   igrf_coef;
[[vk::binding(3, 0)]] StructuredBuffer<float>   igrf_norm;
[[vk::binding(4, 0)]] StructuredBuffer<uint>    igrf_dims;

[shader("compute")]
[numthreads(256, 1, 1)]
void irbem_igrf_f32(uint3 tid : SV_DispatchThreadID) {
    const uint n_points = igrf_dims[0];
    const int nmax = int(igrf_dims[1]);
    const uint i = tid.x;
    if (i >= n_points) return;

    const int slots = ((nmax + 1) * (nmax + 2)) / 2;
    const float3 p = float3(igrf_pos[(3 * i) + 0], igrf_pos[(3 * i) + 1], igrf_pos[(3 * i) + 2]);
    const float3 b = igrf_eval(igrf_coef, igrf_norm, nmax, slots, p);
    igrf_out[(3 * i) + 0] = b.x;
    igrf_out[(3 * i) + 1] = b.y;
    igrf_out[(3 * i) + 2] = b.z;
}
#endif  // KERNEL_IGRF

// =============================================================================================
// irbem_trace_i_f32 — the field-line trace and the second adiabatic invariant, ONE THREAD PER
// FIELD LINE. This is the kernel the whole module exists for.
//
// WHY THE WHOLE LOOP IS IN HERE. The RK4 chain is serial — every step depends on the last — so
// dispatching per field evaluation would cost the ~30 us submit floor for ~37 ns of work, roughly
// 800x WORSE than simply staying on the host. The GPU is not being used to make one trace fast;
// it cannot. It is being used to run 10^5-10^7 INDEPENDENT traces at once, which is the shape the
// problem actually has: ntime x pitch angles x drift-shell azimuths x root-find trials.
//
// WHY IT PAYS, when the dipole kernel did not. Arithmetic intensity, measured on this same seam:
//
//     dipole streaming     0.5 flops/byte     device 243 vs host 352 Mpts/s   ->  0.69x  LOSES
//     IGRF batched         ~20 flops/byte     device 27.3 vs host 3.3 Mpts/s  ->  8.96x  wins
//     this kernel        ~9400 flops/byte     ~120 RK4 steps x 4 IGRF evals x ~500 flops
//                                             for ~28 bytes in and ~16 out
//
// The point is uploaded once and never leaves; only I, Bmin, Bmirr and a status word come back.
// No path is stored — the integral accumulates in registers — which is exactly what keeps this
// compute-bound rather than bandwidth-bound. The APIs that DO return the path (trace_field_line,
// drift_shell) are a different, bandwidth-bound kernel and deliberately not this one.
//
// STEP SIZE IS LOAD-BEARING. ds = L/steps_per_l, with L the start point's dipole L. Arc length AND
// step both scale with L, so the step COUNT stays roughly uniform across lanes carrying an L=2 line
// and an L=8 line — which is what keeps a warp converged. A fixed absolute step is the
// obvious-looking simplification and would make high-L lanes do ~4x the work of low-L ones.
//
// THE CAP IS A MASKED EARLY-OUT, never a divergent branch out of the warp: a line that fails to
// close returns kStatusOpenFieldLine and its lane simply stops contributing, so one pathological
// line cannot stall the other 31.
//
// Mirrors space/irbem/lstar.hpp's host lane step for step, including the two-sided integral: I is
// defined between BOTH mirror points, and integrating one half and doubling it is right for a
// centred dipole and wrong for every real field, which is asymmetric about the equator.
//
// bindings: 0 pos (3N floats, GEO, Re) | 1 pitch (N floats, degrees) | 2 coef | 3 norm
//           4 out (4N floats: I, Bmin, Bmirr, Blocal) | 5 status (N uints) | 6 dims
// dims: {N, nmax, max_steps, steps_per_l_x1000}
// =============================================================================================
#ifdef KERNEL_TRACE

static const uint kStatusOk = 0u;
static const uint kStatusOpenFieldLine = 2u;   // matches Status::OpenFieldLine
static const uint kStatusDomainError = 5u;     // matches Status::DomainError

[[vk::binding(0, 0)]] StructuredBuffer<float>   tr_pos;
[[vk::binding(1, 0)]] StructuredBuffer<float>   tr_pitch;
[[vk::binding(2, 0)]] StructuredBuffer<float>   tr_coef;
[[vk::binding(3, 0)]] StructuredBuffer<float>   tr_norm;
[[vk::binding(4, 0)]] RWStructuredBuffer<float> tr_out;
[[vk::binding(5, 0)]] RWStructuredBuffer<uint>  tr_status;
[[vk::binding(6, 0)]] StructuredBuffer<uint>    tr_dims;

/// The unit tangent B/|B|, or zero where the field vanishes.
float3 tr_unit(float3 b) {
    const float m = length(b);
    return m > 0.0f ? b / m : float3(0.0f, 0.0f, 0.0f);
}

/// One RK4 step along dx/ds = B/|B|. Returns the arrival point and, through `b_next`, the field
/// there — because the next step needs it as its own k1. The reference implementation evaluates a
/// fifth time at the endpoint purely to report |B| and then re-evaluates the same point next call;
/// carrying it forward removes one of every five field evaluations, a 20% cut in the dominant cost.
float3 tr_rk4(int nmax, int slots, float3 p, float3 b_here, float ds, out float3 b_next) {
    const float3 k1 = tr_unit(b_here);
    const float3 k2 = tr_unit(igrf_eval(tr_coef, tr_norm, nmax, slots, p + k1 * (ds * 0.5f)));
    const float3 k3 = tr_unit(igrf_eval(tr_coef, tr_norm, nmax, slots, p + k2 * (ds * 0.5f)));
    const float3 k4 = tr_unit(igrf_eval(tr_coef, tr_norm, nmax, slots, p + k3 * ds));
    const float3 out_p = p + (k1 + k2 * 2.0f + k3 * 2.0f + k4) * (ds / 6.0f);
    b_next = igrf_eval(tr_coef, tr_norm, nmax, slots, out_p);
    return out_p;
}

[shader("compute")]
[numthreads(256, 1, 1)]
void irbem_trace_i_f32(uint3 tid : SV_DispatchThreadID) {
    const uint n_lines = tr_dims[0];
    const int nmax = int(tr_dims[1]);
    const int max_steps = int(tr_dims[2]);
    const float steps_per_l = float(tr_dims[3]) * 0.001f;
    const uint i = tid.x;
    if (i >= n_lines) return;

    const int slots = ((nmax + 1) * (nmax + 2)) / 2;
    const float3 start = float3(tr_pos[3 * i + 0], tr_pos[3 * i + 1], tr_pos[3 * i + 2]);

    tr_out[4 * i + 0] = 0.0f; tr_out[4 * i + 1] = 0.0f;
    tr_out[4 * i + 2] = 0.0f; tr_out[4 * i + 3] = 0.0f;

    const float r0 = length(start);
    if (!(r0 > 1.0f)) { tr_status[i] = kStatusDomainError; return; }

    float3 b = igrf_eval(tr_coef, tr_norm, nmax, slots, start);
    const float b_local = length(b);
    if (!(b_local > 0.0f)) { tr_status[i] = kStatusDomainError; return; }

    const float sin_a = sin(tr_pitch[i] * 0.01745329251994329577f);
    if (!(sin_a > 0.0f)) { tr_status[i] = kStatusDomainError; return; }
    const float b_mirror = b_local / (sin_a * sin_a);

    // Dipole L of the start point, purely to size the step. Floors at 1 so the step cannot collapse.
    const float sin_lat = start.z / r0;
    const float cos2 = 1.0f - sin_lat * sin_lat;
    const float dipole_l = cos2 > 1e-12f ? r0 / cos2 : 1.0e3f;
    const float ds_mag = dipole_l / steps_per_l;

    // --- 1. downhill to the magnetic equator ---------------------------------------------------
    float3 b_probe;
    tr_rk4(nmax, slots, start, b, ds_mag, b_probe);
    float ds = (length(b_probe) < b_local) ? ds_mag : -ds_mag;

    float3 eq = start;
    float3 b_eq = b;
    float b_min = b_local;
    bool open_line = false;
    for (int k = 0; k < max_steps; ++k) {
        float3 b_new;
        const float3 q = tr_rk4(nmax, slots, eq, b_eq, ds, b_new);
        const float bmag = length(b_new);
        if (length(q) < 1.0f) { open_line = true; break; }
        if (bmag >= b_min) break;
        b_min = bmag;
        eq = q;
        b_eq = b_new;
    }

    // --- 2 and 3. both halves of the bounce path, from the equator outward ---------------------
    float total = 0.0f;
    int steps_taken = 0;
    bool closed_both = !open_line;
    for (int half = 0; half < 2; ++half) {
        const float dsh = (half == 0) ? ds : -ds;
        float3 p = eq;
        float3 bb = b_eq;
        bool closed = false;
        for (int k = 0; k < max_steps; ++k) {
            const float b_prev = length(bb);
            float3 b_new;
            const float3 q = tr_rk4(nmax, slots, p, bb, dsh, b_new);
            const float bmag = length(b_new);
            ++steps_taken;
            if (length(q) < 1.0f) break;                 // atmosphere: the line is open
            if (bmag >= b_mirror) {
                if (bmag > b_prev) {
                    const float frac = (b_mirror - b_prev) / (bmag - b_prev);
                    total += 0.5f * sqrt(max(0.0f, 1.0f - (b_prev / b_mirror))) * frac;
                }
                closed = true;
                break;
            }
            total += sqrt(1.0f - (bmag / b_mirror));
            p = q;
            bb = b_new;
        }
        closed_both = closed_both && closed;
    }

    tr_out[4 * i + 0] = total * abs(ds);   // I
    tr_out[4 * i + 1] = b_min;
    tr_out[4 * i + 2] = b_mirror;
    tr_out[4 * i + 3] = b_local;
    tr_status[i] = closed_both ? kStatusOk : kStatusOpenFieldLine;
}

#endif  // KERNEL_TRACE
