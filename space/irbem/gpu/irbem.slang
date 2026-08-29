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

// =============================================================================================
// irbem_shell_foot_f32 — the DRIFT-SHELL FOOTPOINT, one thread per shell azimuth.
//
// WHY THIS KERNEL EXISTS. L* is 2*pi*k0/(Phi*Re) with Phi the flux through the polar cap bounded
// by the drift shell's IONOSPHERIC FOOTPOINTS. Everything upstream of that — the second invariant
// I, the mirror field Bm, the per-azimuth root-find that finds which field line belongs to the
// shell — is served by irbem_trace_i_f32 above. What no kernel yet produced is the one geometric
// fact the flux integral is built on: WHERE each of those field lines meets r = 1.
//
// It cannot be folded into irbem_trace_i_f32. That kernel stops at the MIRROR points, which sit
// well above the atmosphere for every trapped particle, and it stores no path — the whole reason
// it is compute-bound rather than bandwidth-bound. Reaching the surface is a strictly longer walk
// with a different termination test, so it is a different kernel rather than a flag on that one.
//
// WHY IT IS ON THE DEVICE. Per L* point this is Nder = 25 traces of ~250 RK4 steps, four IGRF
// evaluations each: ~25 000 field evaluations, ~7.7 ms on the fp64 host lane and comparable to the
// entire rest of the calculation. It is also perfectly parallel — every azimuth of every point is
// independent — so it is the same bargain the tracer makes: ~16 bytes in, 12 out, ~10^5 flops
// between them.
//
// THE LAST STEP IS REFINED, NOT INTERPOLATED. Marching at a fixed ds and then interpolating the
// r = 1 crossing linearly across the final chord costs O(ds^2) in the footpoint's colatitude,
// which at ds = L/50 is a tenth of a degree — an order of magnitude worse than the 0.25 degree
// quadrature cell it feeds, and a systematic bias rather than noise. So a step that would cross
// the surface is REJECTED and ds is halved instead, repeatedly, until it is below kFootTol. The
// walk therefore lands ON the sphere to ~1e-5 Re along the line, the RK4 truncation stays the only
// error, and the cost is ~2 extra steps per halving — some thirty steps out of ~250.
//
// DIRECTION IS AN INPUT, NOT A GUESS. Which way along B leads to the NORTHERN cap depends on the
// hemisphere the seed sits in and on the sign of g10, and the kernel has no idea where the dipole
// axis points — it works in GEO. The host, which does, passes +1 or -1 per line. That keeps the
// hemisphere choice in one place instead of duplicating the dipole-axis convention on both sides
// of the bus.
//
// bindings: 0 pos (3N floats, GEO seed, Re) | 1 dir (N floats, +1 or -1) | 2 coef | 3 norm
//           4 out (3N floats: the footpoint, GEO, on the UNIT sphere) | 5 status (N uints)
//           6 dims {N, nmax, max_steps, steps_per_l_x1000}
// =============================================================================================
#ifdef KERNEL_SHELL

static const uint kShellOk = 0u;
static const uint kShellOpen = 2u;          // matches Status::OpenFieldLine
static const uint kShellDomain = 5u;        // matches Status::DomainError

/// The step length below which the walk is considered to have landed on the surface, in Re.
/// ~64 m: two orders below the 0.25-degree (28 km) quadrature cell the footpoint feeds, and still
/// far above the fp32 resolution of a coordinate near 1.0 (6e-8 Re).
static const float kFootTol = 1.0e-5f;

[[vk::binding(0, 0)]] StructuredBuffer<float>   sh_pos;
[[vk::binding(1, 0)]] StructuredBuffer<float>   sh_dir;
[[vk::binding(2, 0)]] StructuredBuffer<float>   sh_coef;
[[vk::binding(3, 0)]] StructuredBuffer<float>   sh_norm;
[[vk::binding(4, 0)]] RWStructuredBuffer<float> sh_out;
[[vk::binding(5, 0)]] RWStructuredBuffer<uint>  sh_status;
[[vk::binding(6, 0)]] StructuredBuffer<uint>    sh_dims;

/// The unit tangent B/|B|, or zero where the field vanishes.
float3 sh_unit(float3 b) {
    const float m = length(b);
    return m > 0.0f ? b / m : float3(0.0f, 0.0f, 0.0f);
}

/// One RK4 step along dx/ds = B/|B|, carrying the arrival field forward as the next k1 — the same
/// five-evaluations-into-four saving irbem_trace_i_f32 makes, for the same reason.
float3 sh_rk4(int nmax, int slots, float3 p, float3 b_here, float ds, out float3 b_next) {
    const float3 k1 = sh_unit(b_here);
    const float3 k2 = sh_unit(igrf_eval(sh_coef, sh_norm, nmax, slots, p + k1 * (ds * 0.5f)));
    const float3 k3 = sh_unit(igrf_eval(sh_coef, sh_norm, nmax, slots, p + k2 * (ds * 0.5f)));
    const float3 k4 = sh_unit(igrf_eval(sh_coef, sh_norm, nmax, slots, p + k3 * ds));
    const float3 out_p = p + (k1 + k2 * 2.0f + k3 * 2.0f + k4) * (ds / 6.0f);
    b_next = igrf_eval(sh_coef, sh_norm, nmax, slots, out_p);
    return out_p;
}

[shader("compute")]
[numthreads(256, 1, 1)]
void irbem_shell_foot_f32(uint3 tid : SV_DispatchThreadID) {
    const uint n_lines = sh_dims[0];
    const int nmax = int(sh_dims[1]);
    const int max_steps = int(sh_dims[2]);
    const float steps_per_l = float(sh_dims[3]) * 0.001f;
    const uint i = tid.x;
    if (i >= n_lines) return;

    const int slots = ((nmax + 1) * (nmax + 2)) / 2;
    float3 p = float3(sh_pos[3 * i + 0], sh_pos[3 * i + 1], sh_pos[3 * i + 2]);

    sh_out[3 * i + 0] = 0.0f; sh_out[3 * i + 1] = 0.0f; sh_out[3 * i + 2] = 0.0f;

    const float r0 = length(p);
    if (!(r0 > 1.0f)) { sh_status[i] = kShellDomain; return; }

    float3 b = igrf_eval(sh_coef, sh_norm, nmax, slots, p);
    if (!(length(b) > 0.0f)) { sh_status[i] = kShellDomain; return; }

    // The step is L/steps_per_l for the same reason the tracer's is: arc length and step both
    // scale with L, so lanes carrying an L=2 line and an L=8 line retire together.
    const float sin_lat = p.z / r0;
    const float cos2 = 1.0f - sin_lat * sin_lat;
    const float dipole_l = cos2 > 1e-12f ? r0 / cos2 : 1.0e3f;
    float ds = sh_dir[i] * dipole_l / steps_per_l;

    bool landed = false;
    for (int k = 0; k < max_steps; ++k) {
        float3 b_new;
        const float3 q = sh_rk4(nmax, slots, p, b, ds, b_new);
        if (length(q) < 1.0f) {
            // The surface lies inside this step. Refuse it and halve, until the residual step is
            // smaller than the quadrature can see; only then accept the point below the sphere.
            if (abs(ds) < kFootTol) { p = q; landed = true; break; }
            ds *= 0.5f;
            continue;
        }
        p = q;
        b = b_new;
    }

    // Report the footpoint ON the unit sphere: the caller reads a colatitude and a longitude off
    // it, and normalising here means those are exact rather than carrying the last step's residual.
    const float rp = length(p);
    const float3 foot = rp > 0.0f ? p / rp : float3(0.0f, 0.0f, 1.0f);
    sh_out[3 * i + 0] = foot.x;
    sh_out[3 * i + 1] = foot.y;
    sh_out[3 * i + 2] = foot.z;
    sh_status[i] = landed ? kShellOk : kShellOpen;
}

#endif  // KERNEL_SHELL

// =============================================================================================
// irbem_t89_f32 — Tsyganenko (1989), the external magnetospheric field, one thread per point.
//
// The first EXTERNAL model behind this seam, and the best-shaped one of its family: pure
// straight-line closed-form arithmetic, no loop, no data-dependent branch, ~400 flops with exactly
// one exp, thirteen sqrt and NO trigonometry (the tilt arrives as its sine and cosine, resolved
// once per epoch on the host). Contrast the two kernels above — the streaming dipole at 0.5
// flops/byte LOSES 0.69x, IGRF at ~20 wins 8.96x. T89 is ~11 flops/byte for 24 bytes in and 12
// out, which puts it squarely in the regime where the device pays.
//
// PHYSICS. Tsyganenko, "A magnetospheric magnetic field model with a warped tail current sheet",
// Planet. Space Sci. 37(1):5-20 (1989). Four contributions, each of them the curl of a vector
// potential A = {0, A(rho,z), 0}:
//
//   tail          eqs. (12)-(15)   two disc potentials, thickened, warped, truncated by W(x,y)
//   ring current  eqs. (16)-(17)   the third disc potential, the only one with a finite moment
//   closure       eqs. (18)-(19)   two planar sheets at z = -+R_T, the return current
//   Chapman-      eq.  (20)        a low-order polynomial in y,z under an exp(x/dx) envelope
//   Ferraro
//
// With P = A_phi/rho the curl of A = P*(-y, x, 0) is Bx = -x dP/dz, By = -y dP/dz,
// Bz = 2P + x dP/dx + y dP/dy — which is where the paper's field expressions come from and the
// form evaluated here. Every d/dx and d/dy below is the exact derivative of the line above it;
// dropping any of them breaks div B = 0, which is precisely what the host suite's
// finite-difference divergence test looks for.
//
// TWO COORDINATE SYSTEMS, and the paper says so in as many words: eqs. (12)-(17) are evaluated in
// SM (the dipole axis IS z, because the current sheet is defined by its departure from the dipole
// equator) and eqs. (18)-(20) in GSM (the closure sheets are parallel to the GSM equator). The
// kernel takes GSM in, rotates about y by the tilt, and rotates back. Evaluating the tail term in
// GSM is the obvious-looking simplification and misplaces the sheet by x*tan(psi).
//
// This mirrors space/irbem/ext_t89.hpp's t89_components<float> expression for expression, so a
// disagreement between the two is attributable to the DEVICE (a contraction, a driver's exp) and
// not to the arithmetic having been written differently on the two sides.
//
// bindings: 0 pos (3N floats, GSM, Re) | 1 out (3N floats, nT, GSM)
//           2 params (30 floats: sin psi, cos psi, C1..C19, dx, a_rc, D0, gamma_rc, Rc, G,
//                     a_T, Dy, x0)
//           3 dims {N}
// =============================================================================================
// =============================================================================================
// Tsyganenko 1989, as a SHARED function.
//
// Declared outside every kernel guard and taking its parameter buffer as an argument, so both the
// batched field kernel and the total-field tracer call one implementation of the physics —
// the same arrangement igrf_eval has, for the same reason: the tracer evaluates the field four
// times per RK4 step, and a second copy of eqs. (11)-(20) would be a second place for a sign or a
// slope to be got wrong.
//
// The parameter layout is the one t89_field_host stages: [0]=sin(tilt) [1]=cos(tilt), [2..6] the
// Kp bin's c1..c5, [7..20] the Chapman-Ferraro/Birkeland table, [21..29] delta_x, a_RC, D_0,
// gamma_RC, R_C, G, a_T, D_y, x_0. One bin per buffer — the bin is selected on the HOST, because a
// batch shares one epoch and one Kp, and a per-thread bin branch would diverge the warp for
// nothing.
// =============================================================================================

// The parameters Tsyganenko (1989) sections 3 and 4 hold FIXED across every Kp bin. Spelled out
// here rather than passed, because they are the same in all seven sets and a uniform read of a
// literal costs nothing.
static const float kT89Ly    = 10.0f;   // L_y  : dawn-dusk scale of the sheet's bending, eq. (11)
static const float kT89Dx    = 13.0f;   // D_x  : sunward scale of the truncation factor W
static const float kT89Lrc   = 5.0f;    // L_RC : scale of the ring-current thickness profile
static const float kT89Lt    = 6.3f;    // L_T  : scale of the tail thickness profile
static const float kT89GamT  = 4.0f;    // gamma_T : dayside thickness increment (the pinned one)
static const float kT89Delta = 0.01f;   // delta   : flank thickening rate
static const float kT89Gam1  = 1.0f;    // gamma_1 : far-tail thickening amplitude
static const float kT89Rt    = 30.0f;   // R_T  : the closure sheets sit at z = -+R_T
static const float kT89X0c   = 4.0f;    // x_0c : where W_c makes its transition
static const float kT89Lxc2  = 50.0f;   // L_xc^2, as the paper states it
static const float kT89Dyc2  = 20.0f;   // D_yc^2, likewise
static const float kT89Hinge = 16.0f;   // the (4 Re)^2 rounding the hinge in eq. (11)
static const float kT89H1off = 16.0f;   // the 16 Re offset in h_1
static const float kT89H1sc2 = 36.0f;   // the (6 Re)^2 scale in h_1

/// The T89 external field at a GSM point, in nT. See the layout note above; physics per
/// Tsyganenko (1989), eqs. (11)-(20), with the sheet warp, thickness profiles, truncation factor,
/// tail + ring currents, closure sheets and the Chapman-Ferraro/Birkeland envelope.
float3 t89_eval(StructuredBuffer<float> par, float3 p_gsm) {
    const float sps = par[0];
    const float cps = par[1];
    const float c1 = par[2], c2 = par[3], c3 = par[4];
    const float c4 = par[5], c5 = par[6];
    const float delta_x = par[21], a_rc = par[22], d0 = par[23];
    const float gam_rc  = par[24], rc   = par[25], gg = par[26];
    const float a_t     = par[27], d_y  = par[28], x0 = par[29];

    const float x = p_gsm.x;
    const float y = p_gsm.y;
    const float z = p_gsm.z;
    const float tps = sps / cps;

    // ---- GSM -> SM: one rotation about y by the tilt -----------------------------------------
    const float xs = x * cps - z * sps;
    const float ys = y;
    const float zs = x * sps + z * cps;

    // ---- eq. (11): the warped sheet surface and its slopes -----------------------------------
    const float hinge = xs + rc;
    const float hroot = sqrt(hinge * hinge + kT89Hinge);
    const float y2 = ys * ys;
    const float y4 = y2 * y2;
    const float ly2 = kT89Ly * kT89Ly;
    const float ly4 = ly2 * ly2;
    const float bden = y4 + ly4;
    const float sheet = 0.5f * tps * (hinge - hroot) - gg * sps * y4 / bden;
    const float dsdx = 0.5f * tps * (1.0f - hinge / hroot);
    const float dsdy = -gg * sps * 4.0f * ys * y2 * ly4 / (bden * bden);

    const float zr = zs - sheet;
    const float rho2 = xs * xs + ys * ys;

    // ---- eq. (13): the thickness profiles and their slopes -----------------------------------
    const float lt2 = kT89Lt * kT89Lt;
    const float lrc2 = kT89Lrc * kT89Lrc;
    const float rtr = sqrt(xs * xs + lt2);
    const float rrr = sqrt(xs * xs + lrc2);
    const float xh1 = xs + kT89H1off;
    const float h1r = sqrt(xh1 * xh1 + kT89H1sc2);
    const float h_t  = 0.5f * (1.0f + xs / rtr);
    const float h_rc = 0.5f * (1.0f + xs / rrr);
    const float h_1  = 0.5f * (1.0f - xh1 / h1r);
    const float dh_t  = 0.5f * lt2 / (rtr * rtr * rtr);
    const float dh_rc = 0.5f * lrc2 / (rrr * rrr * rrr);
    const float dh_1  = -0.5f * kT89H1sc2 / (h1r * h1r * h1r);

    const float dT = d0 + kT89Delta * y2 + kT89GamT * h_t + kT89Gam1 * h_1;
    const float dR = d0 + gam_rc * h_rc + kT89Gam1 * h_1;
    const float ddTdx = kT89GamT * dh_t + kT89Gam1 * dh_1;
    const float ddTdy = 2.0f * kT89Delta * ys;
    const float ddRdx = gam_rc * dh_rc + kT89Gam1 * dh_1;

    const float xiT = sqrt(zr * zr + dT * dT);
    const float xiR = sqrt(zr * zr + dR * dR);
    const float uT = a_t + xiT;
    const float uR = a_rc + xiR;
    const float sT = sqrt(rho2 + uT * uT);
    const float sR = sqrt(rho2 + uR * uR);

    // ---- eq. (13): the truncation factor W and its slopes ------------------------------------
    const float xw = xs - x0;
    const float dx2 = kT89Dx * kT89Dx;
    const float wroot = sqrt(xw * xw + dx2);
    const float wx = 0.5f * (1.0f - xw / wroot);
    const float wy = 1.0f / (1.0f + y2 / (d_y * d_y));
    const float W = wx * wy;
    const float dWdx = wy * (-0.5f * dx2 / (wroot * wroot * wroot));
    const float dWdy = wx * (-2.0f * ys / (d_y * d_y)) * wy * wy;
    const float gradW = xs * dWdx + ys * dWdy;

    // ---- eqs. (14)-(15): the tail current, in SM ----------------------------------------------
    const float qT = (W / (xiT * sT)) * (c1 / (sT + uT) + c2 / (sT * sT));
    const float bxT = qT * xs * zr;
    const float byT = qT * ys * zr;
    const float bzT = (W / sT) * (c1 + c2 * uT / (sT * sT))
                    + gradW * (c1 + c2 / sT) / (sT + uT)
                    + bxT * dsdx + byT * dsdy
                    - qT * dT * (xs * ddTdx + ys * ddTdy);

    // ---- eqs. (16)-(17): the ring current, in SM ----------------------------------------------
    const float sR2 = sR * sR;
    const float sR5 = sR2 * sR2 * sR;
    const float qR = 3.0f * c5 * uR / (xiR * sR5);
    const float bxR = qR * xs * zr;
    const float byR = qR * ys * zr;
    const float bzR = c5 * (2.0f * uR * uR - rho2) / sR5
                    + bxR * dsdx + byR * dsdy - qR * dR * xs * ddRdx;

    // ---- SM -> GSM ----------------------------------------------------------------------------
    const float bxsm = bxT + bxR, bysm = byT + byR, bzsm = bzT + bzR;
    float bx = bxsm * cps + bzsm * sps;
    float by = bysm;
    float bz = -bxsm * sps + bzsm * cps;

    // ---- eqs. (18)-(19): the closure currents, in GSM ------------------------------------------
    // The paper's -+ resolves |z -+ R_T| for a point BETWEEN the sheets; sgn(zeta) is identical
    // there and keeps the expression finite outside, where the paper's spelling divides by zero.
    float fx[2], fy[2], fz[2];
    for (int k = 0; k < 2; ++k) {
        const float side = (k == 0) ? 1.0f : -1.0f;
        const float zeta = z + side * kT89Rt;
        const float S = sqrt(zeta * zeta + x * x + y * y);
        const float sgn = (zeta >= 0.0f) ? 1.0f : -1.0f;
        const float den = S + sgn * zeta;
        const float xc = x - kT89X0c;
        const float sc = sqrt(xc * xc + kT89Lxc2);
        const float wcx = 0.5f * (1.0f - xc / sc);
        const float wcy = 1.0f / (1.0f + y * y / kT89Dyc2);
        const float wc = wcx * wcy;
        const float dwcdx = wcy * (-0.5f * kT89Lxc2 / (sc * sc * sc));
        const float dwcdy = wcx * (-2.0f * y / kT89Dyc2) * wcy * wcy;
        const bool ok = den > 0.0f;
        fx[k] = ok ? sgn * wc * x / (S * den) : 0.0f;
        fy[k] = ok ? sgn * wc * y / (S * den) : 0.0f;
        fz[k] = ok ? (wc / S + (x * dwcdx + y * dwcdy) / den) : 0.0f;
    }
    bx += c3 * (fx[0] + fx[1]) + c4 * (fx[0] - fx[1]) * sps;
    by += c3 * (fy[0] + fy[1]) + c4 * (fy[0] - fy[1]) * sps;
    bz += c3 * (fz[0] + fz[1]) + c4 * (fz[0] - fz[1]) * sps;

    // ---- eq. (20): the Chapman-Ferraro and Birkeland systems, in GSM ---------------------------
    const float env = exp(x / delta_x);
    const float z2 = z * z;
    bx += env * (par[7] * z * cps
               + (par[8] + par[9] * y2 + par[10] * z2) * sps);
    by += env * (par[11] * y * z * cps
               + (par[12] * y + par[13] * y * y2 + par[14] * y * z2) * sps);
    bz += env * ((par[15] + par[16] * y2 + par[17] * z2) * cps
               + (par[18] * z + par[19] * z * y2 + par[20] * z * z2) * sps);

    return float3(bx, by, bz);
}

#ifdef KERNEL_T89

[[vk::binding(0, 0)]] StructuredBuffer<float>   t89_pos;
[[vk::binding(1, 0)]] RWStructuredBuffer<float> t89_out;
[[vk::binding(2, 0)]] StructuredBuffer<float>   t89_params;
[[vk::binding(3, 0)]] StructuredBuffer<uint>    t89_dims;

[shader("compute")]
[numthreads(256, 1, 1)]
void irbem_t89_f32(uint3 tid : SV_DispatchThreadID) {
    const uint n = t89_dims[0];
    const uint i = tid.x;
    if (i >= n) return;   // the grid is rounded up to a multiple of 256; the tail does nothing
    const float3 b = t89_eval(t89_params,
                              float3(t89_pos[3 * i + 0], t89_pos[3 * i + 1], t89_pos[3 * i + 2]));
    t89_out[3 * i + 0] = b.x;
    t89_out[3 * i + 1] = b.y;
    t89_out[3 * i + 2] = b.z;
}

#endif  // KERNEL_T89

// =============================================================================================
// irbem_trace_path_f32 — the field line WITH its path, one thread per line.
//
// THE OTHER SIDE OF THE ROOFLINE. irbem_trace_i_f32 above returns four floats per line for ~10^5
// field evaluations — ~9 400 flops per byte moved, which is why it beats the host by 48.9x. This
// kernel runs the SAME RK4 chain and writes every sample: 16 bytes per step instead of 16 bytes
// per line. At the ~4 IGRF evaluations and ~500 flops each that one step costs, that is ~125
// flops/byte — a factor of 75 down, and the crossover moves with it (measured: 4096 lines, versus
// the tracer's 512). Two kernels, not one, because merging them would put a `posit` store in the
// L* hot path and cost the module its whole performance argument.
//
// Implements space/irbem/trace_api.hpp's trace_field_line_toward_earth: a FIXED step (not the
// L-proportional one — a caller asking for a path is sampling for a plot and wants uniform
// spacing), in the direction of increasing |B|, from the input point to the first sample inside
// R0. The direction is chosen by probing BOTH neighbours and taking the larger |B|; comparing one
// neighbour against the start picks the wrong hemisphere within a step of the magnetic equator,
// which the host lane's own history records.
//
// fp32 is a real choice here and not a free one. The tracer above accumulates one scalar; this
// one accumulates a POSITION over hundreds of steps, and each step's rounding is carried forward.
// MEASURED against the fp64 host lane, 512 lines x ~250 steps at ds = 0.02 Re: 5.0e-7 R_E at the
// first sample and 6.2e-6 R_E (~40 metres) at the last, i.e. the error GROWS along the path, and
// |B| agrees to 2.0e-5 relative. Sample counts were identical on every line. That last-sample
// figure is a factor of 28 inside the path-geometry agreement the fp64 host lane itself has with
// the reference implementation (1.8e-4 R_E), so this lane is fit for plotting and for a
// drift-shell sweep -- and it is NOT fit for anything that differentiates the path, where 4e-5
// relative on a 0.02 Re step is first-order noise. Said here so it is not discovered later.
//
// bindings: 0 pos (3N floats, GEO, Re) | 1 coef | 2 norm
//           3 path (3*M*N floats, GEO, Re) | 4 bmag (M*N floats, nT)
//           5 report (2N uints: sample count, status) | 6 dims
// dims: {N, nmax, M, ds_x1e6, r0_x1e6}
// =============================================================================================
#ifdef KERNEL_TRACE_PATH

static const uint kPathStatusOk = 0u;
static const uint kPathStatusNotConverged = 3u;   // matches Status::NotConverged
static const uint kPathStatusOpenFieldLine = 2u;  // matches Status::OpenFieldLine
static const uint kPathStatusDomainError = 5u;    // matches Status::DomainError

[[vk::binding(0, 0)]] StructuredBuffer<float>   tp_pos;
[[vk::binding(1, 0)]] StructuredBuffer<float>   tp_coef;
[[vk::binding(2, 0)]] StructuredBuffer<float>   tp_norm;
[[vk::binding(3, 0)]] RWStructuredBuffer<float> tp_path;
[[vk::binding(4, 0)]] RWStructuredBuffer<float> tp_bmag;
[[vk::binding(5, 0)]] RWStructuredBuffer<uint>  tp_report;
[[vk::binding(6, 0)]] StructuredBuffer<uint>    tp_dims;

/// The unit tangent B/|B|, or zero where the field vanishes.
float3 tp_unit(float3 b) {
    const float m = length(b);
    return m > 0.0f ? b / m : float3(0.0f, 0.0f, 0.0f);
}

/// One RK4 step along dx/ds = B/|B|, carrying the arrival field forward as the next step's k1 —
/// the same 20% saving the invariant tracer takes.
float3 tp_rk4(int nmax, int slots, float3 p, float3 b_here, float ds, out float3 b_next) {
    const float3 k1 = tp_unit(b_here);
    const float3 k2 = tp_unit(igrf_eval(tp_coef, tp_norm, nmax, slots, p + k1 * (ds * 0.5f)));
    const float3 k3 = tp_unit(igrf_eval(tp_coef, tp_norm, nmax, slots, p + k2 * (ds * 0.5f)));
    const float3 k4 = tp_unit(igrf_eval(tp_coef, tp_norm, nmax, slots, p + k3 * ds));
    const float3 out_p = p + (k1 + k2 * 2.0f + k3 * 2.0f + k4) * (ds / 6.0f);
    b_next = igrf_eval(tp_coef, tp_norm, nmax, slots, out_p);
    return out_p;
}

[shader("compute")]
[numthreads(256, 1, 1)]
void irbem_trace_path_f32(uint3 tid : SV_DispatchThreadID) {
    const uint n_lines = tp_dims[0];
    const int nmax = int(tp_dims[1]);
    const uint max_points = tp_dims[2];
    const float ds_mag = float(tp_dims[3]) * 1.0e-6f;
    const float r0 = float(tp_dims[4]) * 1.0e-6f;
    // PathTraceOptions::max_steps. It rides in the dims buffer because the host loop is bounded by
    // it and not by the buffer: without it a caller who caps the step count gets one answer with a
    // device present and another without, which is exactly the silent lane-dependence the batch
    // entry point exists to rule out. Measured before it was carried: at max_steps = 10 all 1024
    // lines of a batch disagreed (device 80 samples / Ok, host 11 / OpenFieldLine).
    const uint max_steps = tp_dims[5];
    const uint i = tid.x;
    if (i >= n_lines) return;

    const int slots = ((nmax + 1) * (nmax + 2)) / 2;
    const uint base = i * max_points;
    const float3 start = float3(tp_pos[3 * i + 0], tp_pos[3 * i + 1], tp_pos[3 * i + 2]);

    tp_report[2 * i + 0] = 0u;
    tp_report[2 * i + 1] = kPathStatusDomainError;
    if (!(length(start) > r0) || max_points == 0u) return;

    float3 b = igrf_eval(tp_coef, tp_norm, nmax, slots, start);
    const float b_start = length(b);
    if (!(b_start > 0.0f)) return;

    // Direction: probe both neighbours, take the larger |B|. See the block comment above.
    float3 up;
    float3 down;
    tp_rk4(nmax, slots, start, b, ds_mag, up);
    tp_rk4(nmax, slots, start, b, -ds_mag, down);
    const float ds = length(up) >= length(down) ? ds_mag : -ds_mag;

    tp_path[3 * base + 0] = start.x;
    tp_path[3 * base + 1] = start.y;
    tp_path[3 * base + 2] = start.z;
    tp_bmag[base] = b_start;
    uint n = 1u;

    float3 p = start;
    uint status = kPathStatusOpenFieldLine;
    for (uint k = 0u; k < max_steps; ++k) {
        if (n == max_points) { status = kPathStatusNotConverged; break; }
        float3 b_new;
        const float3 q = tp_rk4(nmax, slots, p, b, ds, b_new);
        tp_path[3 * (base + n) + 0] = q.x;
        tp_path[3 * (base + n) + 1] = q.y;
        tp_path[3 * (base + n) + 2] = q.z;
        tp_bmag[base + n] = length(b_new);
        ++n;
        if (length(q) < r0) { status = kPathStatusOk; break; }
        p = q;
        b = b_new;
    }

    tp_report[2 * i + 0] = n;
    tp_report[2 * i + 1] = status;
}

#endif  // KERNEL_TRACE_PATH

// =============================================================================================
// irbem_trace_total_f32 — the field-line trace through the TOTAL field: IGRF plus T89, one
// thread per field line.
//
// This is the kernel that closes the gap the host-side TotalFieldT89 documented: until it
// existed, a trace through internal-plus-external ran CPU-only, so every GPU speedup this module
// quoted was for the internal field alone and the storm-time claims ran at host speed.
//
// STRUCTURE. The PHYSICS is shared — igrf_eval and t89_eval are the same functions the batched
// field kernels call, so there is exactly one copy of the Legendre recursion and one copy of
// Tsyganenko's eqs. (11)-(20) in this file. The trace CONTROL FLOW (downhill walk to the equator,
// then both halves of the bounce path) is a copy of irbem_trace_i_f32's, because Slang gives no
// cheap way to parameterize one loop over two field evaluators; that duplication is accepted
// deliberately, it is the smaller risk (simple control flow, no physics), and the differential
// test holds both kernels to the same fp64 host lane.
//
// FRAMES. The trace walks in GEO (the frame the invariants and footpoints are defined in); T89 is
// defined in GSM. The GEO->GSM rotation depends only on the epoch, so the HOST builds it once and
// uploads 9 floats; being orthogonal, its transpose — applied to the returned field vector — is
// the inverse, so one matrix serves both directions. A FieldVector rotates with no origin shift,
// which is why the same 9 floats can carry both the position in and the field back out.
//
// THE EXT BUFFER: [0..29] the T89 parameter block exactly as t89_eval reads it (one Kp bin,
// selected on the HOST — a batch shares one epoch and one Kp, and a per-thread bin branch would
// diverge the warp for nothing); [30..38] gsm_from_geo, column-major.
//
// EIGHT BINDINGS — the shared descriptor layout's cap, and the reason Leases::capacity is 8.
//
// bindings: 0 pos (3N, GEO, Re) | 1 pitch (N, deg) | 2 igrf coef | 3 igrf norm
//           4 ext (39 floats, above) | 5 out (4N: I, Bmin, Bmirr, Blocal) | 6 status (N) | 7 dims
// dims: {N, nmax, max_steps, steps_per_l_x1000}
// =============================================================================================
#ifdef KERNEL_TRACE_TOTAL

static const uint kTotStatusOk = 0u;
static const uint kTotStatusOpenFieldLine = 2u;   // matches Status::OpenFieldLine
static const uint kTotStatusDomainError = 5u;     // matches Status::DomainError

[[vk::binding(0, 0)]] StructuredBuffer<float>   tt_pos;
[[vk::binding(1, 0)]] StructuredBuffer<float>   tt_pitch;
[[vk::binding(2, 0)]] StructuredBuffer<float>   tt_coef;
[[vk::binding(3, 0)]] StructuredBuffer<float>   tt_norm;
[[vk::binding(4, 0)]] StructuredBuffer<float>   tt_ext;
[[vk::binding(5, 0)]] RWStructuredBuffer<float> tt_out;
[[vk::binding(6, 0)]] RWStructuredBuffer<uint>  tt_status;
[[vk::binding(7, 0)]] StructuredBuffer<uint>    tt_dims;

/// The total field at a GEO point: internal plus external, in nT.
float3 tt_field(int nmax, int slots, float3 p_geo) {
    const float3 b_int = igrf_eval(tt_coef, tt_norm, nmax, slots, p_geo);
    // gsm_from_geo, column-major at ext[30..38]: column j holds R[:,j].
    const float3 c0 = float3(tt_ext[30], tt_ext[31], tt_ext[32]);
    const float3 c1 = float3(tt_ext[33], tt_ext[34], tt_ext[35]);
    const float3 c2 = float3(tt_ext[36], tt_ext[37], tt_ext[38]);
    const float3 p_gsm = c0 * p_geo.x + c1 * p_geo.y + c2 * p_geo.z;
    const float3 b_gsm = t89_eval(tt_ext, p_gsm);
    // The transpose IS the inverse (orthogonal): row i of R = (c0[i], c1[i], c2[i]).
    const float3 b_ext = float3(dot(c0, b_gsm), dot(c1, b_gsm), dot(c2, b_gsm));
    return b_int + b_ext;
}

float3 tt_unit(float3 b) {
    const float m = length(b);
    return m > 0.0f ? b / m : float3(0.0f, 0.0f, 0.0f);
}

/// One RK4 step along dx/ds = B/|B| of the TOTAL field, carrying the arrival field forward as the
/// next step's k1 — one of every five evaluations removed, as in the internal-field tracer.
float3 tt_rk4(int nmax, int slots, float3 p, float3 b_here, float ds, out float3 b_next) {
    const float3 k1 = tt_unit(b_here);
    const float3 k2 = tt_unit(tt_field(nmax, slots, p + k1 * (ds * 0.5f)));
    const float3 k3 = tt_unit(tt_field(nmax, slots, p + k2 * (ds * 0.5f)));
    const float3 k4 = tt_unit(tt_field(nmax, slots, p + k3 * ds));
    const float3 out_p = p + (k1 + k2 * 2.0f + k3 * 2.0f + k4) * (ds / 6.0f);
    b_next = tt_field(nmax, slots, out_p);
    return out_p;
}

[shader("compute")]
[numthreads(256, 1, 1)]
void irbem_trace_total_f32(uint3 tid : SV_DispatchThreadID) {
    const uint n_lines = tt_dims[0];
    const int nmax = int(tt_dims[1]);
    const int max_steps = int(tt_dims[2]);
    const float steps_per_l = float(tt_dims[3]) * 0.001f;
    const uint i = tid.x;
    if (i >= n_lines) return;

    const int slots = ((nmax + 1) * (nmax + 2)) / 2;
    const float3 start = float3(tt_pos[3 * i + 0], tt_pos[3 * i + 1], tt_pos[3 * i + 2]);

    tt_out[4 * i + 0] = 0.0f; tt_out[4 * i + 1] = 0.0f;
    tt_out[4 * i + 2] = 0.0f; tt_out[4 * i + 3] = 0.0f;

    const float r0 = length(start);
    if (!(r0 > 1.0f)) { tt_status[i] = kTotStatusDomainError; return; }

    float3 b = tt_field(nmax, slots, start);
    const float b_local = length(b);
    if (!(b_local > 0.0f)) { tt_status[i] = kTotStatusDomainError; return; }

    const float sin_a = sin(tt_pitch[i] * 0.01745329251994329577f);
    if (!(sin_a > 0.0f)) { tt_status[i] = kTotStatusDomainError; return; }
    const float b_mirror = b_local / (sin_a * sin_a);

    // Dipole L of the start, purely to size the step: ds = L/steps_per_l keeps step COUNTS uniform
    // across lanes carrying different shells, which is what keeps a warp converged.
    const float sin_lat = start.z / r0;
    const float cos2 = 1.0f - sin_lat * sin_lat;
    const float dipole_l = cos2 > 1e-12f ? r0 / cos2 : 1.0e3f;
    const float ds_mag = dipole_l / steps_per_l;

    // --- 1. downhill to the magnetic equator ---------------------------------------------------
    float3 b_probe;
    tt_rk4(nmax, slots, start, b, ds_mag, b_probe);
    float ds = (length(b_probe) < b_local) ? ds_mag : -ds_mag;

    float3 eq = start;
    float3 b_eq = b;
    float b_min = b_local;
    bool open_line = false;
    for (int k = 0; k < max_steps; ++k) {
        float3 b_new;
        const float3 q = tt_rk4(nmax, slots, eq, b_eq, ds, b_new);
        const float bmag = length(b_new);
        if (length(q) < 1.0f) { open_line = true; break; }
        if (bmag >= b_min) break;
        b_min = bmag;
        eq = q;
        b_eq = b_new;
    }

    // --- 2 and 3. both halves of the bounce path, from the equator outward ---------------------
    float total = 0.0f;
    bool closed_both = !open_line;
    for (int half = 0; half < 2; ++half) {
        const float dsh = (half == 0) ? ds : -ds;
        float3 p = eq;
        float3 bb = b_eq;
        bool closed = false;
        for (int k = 0; k < max_steps; ++k) {
            const float b_prev = length(bb);
            float3 b_new;
            const float3 q = tt_rk4(nmax, slots, p, bb, dsh, b_new);
            const float bmag = length(b_new);
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

    tt_out[4 * i + 0] = total * abs(ds);   // I
    tt_out[4 * i + 1] = b_min;
    tt_out[4 * i + 2] = b_mirror;
    tt_out[4 * i + 3] = b_local;
    tt_status[i] = closed_both ? kTotStatusOk : kTotStatusOpenFieldLine;
}

#endif  // KERNEL_TRACE_TOTAL

// =============================================================================================
// irbem_opd_f32 — the Olson-Pfitzer DYNAMIC external field (kext = 6), one thread per point.
//
// The documented structure of the 1988 model, as space/irbem/ext_opd.hpp states and cites it:
//
//     B(r) = s^3 * B_q(s r) + C * R(r)
//
// a quiet field compressed self-similarly by s = (P / P_0)^(1/6), plus a flat symmetric ring
// current in the dipole equator whose amplitude C is the pressure-corrected Dst's departure from
// quiet. The quiet field is the SHARED t89_eval on the scaled position — the parameter buffer's
// first thirty floats are exactly the T89 layout for the quietest published set — and the ring is
// the A^(3) disc of Tsyganenko (1989) eq. (9), unwarped, evaluated in SM and rotated back. This
// mirrors ext_opd.hpp's opd_components<float> expression for expression, so a disagreement between
// the two lanes is attributable to the DEVICE.
//
// bindings: 0 pos (3N floats, GSM, Re) | 1 out (3N floats, nT, GSM)
//           2 params (32 floats: [0..29] the T89 quiet block as t89_eval reads it —
//                     sin psi, cos psi, C1..C19, dx, a_rc, D0, gamma_rc, Rc, G, a_T, Dy, x0 —
//                     then [30] s, [31] C)
//           3 dims {N}
// =============================================================================================
#ifdef KERNEL_OPD

[[vk::binding(0, 0)]] StructuredBuffer<float>   opd_pos;
[[vk::binding(1, 0)]] RWStructuredBuffer<float> opd_out;
[[vk::binding(2, 0)]] StructuredBuffer<float>   opd_params;
[[vk::binding(3, 0)]] StructuredBuffer<uint>    opd_dims;

/// The Dst-driven ring per unit amplitude: Tsyganenko (1989) eqs. (16)-(17) with z_s = 0 and a
/// constant thickness, so every slope term of the published form is identically zero. Radius and
/// half-thickness are the quiet set's a_RC and D_0, read from the T89 block.
float3 opd_ring_unit(StructuredBuffer<float> par, float3 p) {
    const float sps = par[0];
    const float cps = par[1];
    const float a = par[22];
    const float d = par[23];
    const float xs = p.x * cps - p.z * sps;
    const float ys = p.y;
    const float zs = p.x * sps + p.z * cps;
    const float rho2 = xs * xs + ys * ys;
    const float xi = sqrt(zs * zs + d * d);
    const float u = a + xi;
    const float S = sqrt(rho2 + u * u);
    const float S2 = S * S;
    const float S5 = S2 * S2 * S;
    const float q = 3.0f * u / (xi * S5);
    const float bxsm = q * xs * zs;
    const float bysm = q * ys * zs;
    const float bzsm = (2.0f * u * u - rho2) / S5;
    return float3(bxsm * cps + bzsm * sps, bysm, -bxsm * sps + bzsm * cps);
}

[shader("compute")]
[numthreads(256, 1, 1)]
void irbem_opd_f32(uint3 tid : SV_DispatchThreadID) {
    const uint n = opd_dims[0];
    const uint i = tid.x;
    if (i >= n) return;   // the grid is rounded up to a multiple of 256; the tail does nothing
    const float s = opd_params[30];
    const float c = opd_params[31];
    const float s3 = s * s * s;
    const float3 p = float3(opd_pos[3 * i + 0], opd_pos[3 * i + 1], opd_pos[3 * i + 2]);
    const float3 quiet = t89_eval(opd_params, float3(s * p.x, s * p.y, s * p.z));
    const float3 ring = opd_ring_unit(opd_params, p);
    opd_out[3 * i + 0] = s3 * quiet.x + c * ring.x;
    opd_out[3 * i + 1] = s3 * quiet.y + c * ring.y;
    opd_out[3 * i + 2] = s3 * quiet.z + c * ring.z;
}

#endif  // KERNEL_OPD

// =============================================================================================
// irbem_mead_f32 — Mead & Fairfield (1975) external field, one thread per point.
//
// PHYSICS. Mead & Fairfield, J. Geophys. Res. 80(4):523 (1975): the external field as a quadratic
// in the ABERRATED solar-magnetic position and a linear function of the dipole tilt in DEGREES,
// four Kp-binned coefficient sets (the paper's Tables 2-3). With (xm, ym, z) the SM position
// rotated about the dipole axis by the 4-degree mean solar-wind aberration and t the tilt in
// degrees:
//
//     Bx = a1 z + a2 xm z + t (a3 + a4 xm + a5 xm^2 + a6 ym^2 + a7 z^2)
//     By = b1 ym z + t (b2 ym + b3 xm ym)
//     Bz = c1 + c2 xm + c3 xm^2 + c4 ym^2 + c5 z^2 + t (c6 z + c7 xm z)
//
// The aberration rotates the POSITION only; the components come out in SM and are rotated back
// to GSM by the tilt alone — the file brief of space/irbem/ext_mead.hpp carries the oracle
// experiment that establishes this. About 50 flops: no loop, no branch, no transcendental.
//
// This mirrors ext_mead.hpp's mead_components<float> expression for expression, so a
// disagreement between the two is attributable to the DEVICE and not to the arithmetic having
// been written differently on the two sides.
//
// bindings: 0 pos (3N floats, GSM, Re) | 1 out (3N floats, nT, GSM)
//           2 params (20 floats: sin psi, cos psi, psi in degrees, a1..a7, b1..b3, c1..c7)
//           3 dims {N}
// =============================================================================================

// sin(4 deg) and cos(4 deg), the one parameter the paper fixed — the same literals ext_mead.hpp
// carries, rounded to float here exactly as its float lane rounds them.
static const float kMeadAbSin = 0.0697564737441253f;
static const float kMeadAbCos = 0.9975640502598242f;

/// The Mead-Fairfield external field at a GSM point, in nT. Shared with any future total-field
/// tracer the way t89_eval is, so the physics exists in this file exactly once. Parameter layout
/// per the block comment above.
float3 mead_eval(StructuredBuffer<float> par, float3 p_gsm) {
    const float sps = par[0];
    const float cps = par[1];
    const float t   = par[2];

    // ---- GSM -> SM: one rotation about y by the tilt -----------------------------------------
    const float xs = p_gsm.x * cps - p_gsm.z * sps;
    const float ys = p_gsm.y;
    const float zs = p_gsm.x * sps + p_gsm.z * cps;

    // ---- the aberrated position: SM rotated about the dipole axis by 4 degrees ---------------
    const float xm = xs * kMeadAbCos - ys * kMeadAbSin;
    const float ym = xs * kMeadAbSin + ys * kMeadAbCos;

    // ---- the three polynomials, in SM components ---------------------------------------------
    const float xm2 = xm * xm;
    const float ym2 = ym * ym;
    const float zs2 = zs * zs;
    const float bx_sm = par[3] * zs + par[4] * xm * zs
                      + t * (par[5] + par[6] * xm + par[7] * xm2 + par[8] * ym2 + par[9] * zs2);
    const float by_sm = par[10] * ym * zs + t * (par[11] * ym + par[12] * xm * ym);
    const float bz_sm = par[13] + par[14] * xm + par[15] * xm2 + par[16] * ym2 + par[17] * zs2
                      + t * (par[18] * zs + par[19] * xm * zs);

    // ---- SM -> GSM ----------------------------------------------------------------------------
    return float3(bx_sm * cps + bz_sm * sps, by_sm, -bx_sm * sps + bz_sm * cps);
}

#ifdef KERNEL_MEAD

[[vk::binding(0, 0)]] StructuredBuffer<float>   mead_pos;
[[vk::binding(1, 0)]] RWStructuredBuffer<float> mead_out;
[[vk::binding(2, 0)]] StructuredBuffer<float>   mead_params;
[[vk::binding(3, 0)]] StructuredBuffer<uint>    mead_dims;

[shader("compute")]
[numthreads(256, 1, 1)]
void irbem_mead_f32(uint3 tid : SV_DispatchThreadID) {
    const uint n = mead_dims[0];
    const uint i = tid.x;
    if (i >= n) return;   // the grid is rounded up to a multiple of 256; the tail does nothing
    const float3 b = mead_eval(mead_params,
                               float3(mead_pos[3 * i + 0], mead_pos[3 * i + 1], mead_pos[3 * i + 2]));
    mead_out[3 * i + 0] = b.x;
    mead_out[3 * i + 1] = b.y;
    mead_out[3 * i + 2] = b.z;
}

#endif  // KERNEL_MEAD

// =============================================================================================
// Olson & Pfitzer (1977) quiet-time external field, as a SHARED function.
//
// Declared outside every kernel guard and taking its parameter buffer as an argument, for the
// same reason igrf_eval and t89_eval are: one implementation of the physics on the device, so a
// tracer that later composes it evaluates the same series the batched kernel does.
//
// The model is the report's subroutine BXYZMU (Olson & Pfitzer, AFOSR F44620-75-C-0033, 1977,
// pp. 64-67): a power series in SOLAR-MAGNETIC position through fifth order with an
// exp(-0.06 r^2) envelope, its coefficients already folded for the epoch's tilt ON THE HOST
// (the fold is a cubic in the tilt in degrees; a batch shares one epoch, so folding per thread
// would be 172 redundant pairs per point). Three radial rules from the report: zero inside 2 Re,
// a linear-in-r^2 taper to 2.5 Re, and the zero template beyond 15 Re. The enumeration of the
// 32 + 22 + 32 terms is the report's own loop with its two truncation counters, mirrored
// statement for statement from space/irbem/ext_opq.hpp's opq_components<float>, so that a
// device-vs-host disagreement is attributable to the DEVICE and not to the arithmetic having
// been written differently on the two sides.
//
// The parameter layout is the one opq_param_block stages: [0]=sin(tilt) [1]=cos(tilt),
// [2..33]=A(32), [34..65]=B(32), [66..87]=C(22), [88..109]=D(22), [110..141]=E(32),
// [142..173]=F(32) — 174 floats. The kernel takes GSM in, rotates about y by the tilt into SM,
// evaluates, and rotates the field back.
//
// bindings: 0 pos (3N floats, GSM, Re) | 1 out (3N floats, nT, GSM)
//           2 params (174 floats, as above) | 3 dims {N}
// =============================================================================================

static const float kOpqExpRate    = 0.06f;    // exp(-0.06 r^2), report p. 67
static const float kOpqInnerR2    = 4.0f;     // zero inside 2 Re, report p. 66
static const float kOpqTaperR2    = 6.25f;    // linear taper to 2.5 Re
static const float kOpqTemplateR2 = 225.0f;   // the zero template beyond 15 Re

/// The OP-77 external field at an SM point, in nT, from a folded parameter buffer.
float3 opq_eval_sm(StructuredBuffer<float> par, float3 p_sm) {
    const float x = p_sm.x;
    const float y = p_sm.y;
    const float z = p_sm.z;
    const float y2 = y * y;
    const float r2 = x * x + y2 + z * z;
    if (!(r2 <= kOpqTemplateR2) || r2 < kOpqInnerR2) return float3(0.0f, 0.0f, 0.0f);
    const float con = r2 < kOpqTaperR2 ? (r2 - kOpqInnerR2) / (kOpqTaperR2 - kOpqInnerR2) : 1.0f;
    const float expr = exp(-kOpqExpRate * r2);

    float bx = 0.0f;
    float by = 0.0f;
    float bz = 0.0f;
    int ix = 0;        // the next B_x / B_z term: A at par[2+ix], B at par[34+ix], E at
                       // par[110+ix], F at par[142+ix]
    int iy = 0;        // the next B_y term: C at par[66+iy], D at par[88+iy]
    float xb = 1.0f;
    for (int i = 1; i <= 5; ++i) {
        float yexb = xb;
        for (int j = 1; j <= 3; ++j) {
            if (i + 2 * j > 8) break;
            int ijk = i + 2 * j + 1;
            int k = 0;
            float zeyexb = yexb;
            for (;;) {
                bz += (par[110 + ix] + par[142 + ix] * expr) * zeyexb;
                bx += (par[2 + ix] + par[34 + ix] * expr) * zeyexb;
                ++ix;
                if (ijk > 8) break;
                by += (par[66 + iy] + par[88 + iy] * expr) * zeyexb * y;
                ++iy;
                zeyexb *= z;
                ++ijk;
                ++k;
                if (!(ijk <= 9 && k <= 4)) break;
            }
            yexb *= y2;
        }
        xb *= x;
    }
    return float3(bx * con, by * con, bz * con);
}

/// The OP-77 external field at a GSM point, in nT: GSM -> SM, the series, SM -> GSM.
float3 opq_eval(StructuredBuffer<float> par, float3 p_gsm) {
    const float sps = par[0];
    const float cps = par[1];
    const float xs = p_gsm.x * cps - p_gsm.z * sps;
    const float zs = p_gsm.x * sps + p_gsm.z * cps;
    const float3 b = opq_eval_sm(par, float3(xs, p_gsm.y, zs));
    return float3(b.x * cps + b.z * sps, b.y, -b.x * sps + b.z * cps);
}

#ifdef KERNEL_OPQ

[[vk::binding(0, 0)]] StructuredBuffer<float>   opq_pos;
[[vk::binding(1, 0)]] RWStructuredBuffer<float> opq_out;
[[vk::binding(2, 0)]] StructuredBuffer<float>   opq_params;
[[vk::binding(3, 0)]] StructuredBuffer<uint>    opq_dims;

[shader("compute")]
[numthreads(256, 1, 1)]
void irbem_opq_f32(uint3 tid : SV_DispatchThreadID) {
    const uint n = opq_dims[0];
    const uint i = tid.x;
    if (i >= n) return;   // the grid is rounded up to a multiple of 256; the tail does nothing
    const float3 b = opq_eval(opq_params,
                              float3(opq_pos[3 * i + 0], opq_pos[3 * i + 1], opq_pos[3 * i + 2]));
    opq_out[3 * i + 0] = b.x;
    opq_out[3 * i + 1] = b.y;
    opq_out[3 * i + 2] = b.z;
}

#endif  // KERNEL_OPQ

// =============================================================================================
// irbem_om97_f32 — Ostapenko & Maltsev (1997), one thread per point.
//
// The external field as eq. (2) of the paper: 17 polynomial harmonics (Table 1) in solar
// magnetic coordinates normalized by 10 Re, weighted by 17 amplitudes that the HOST has already
// regressed from the drivers (eq. (3), Tables 2 and 4). The harmonics are the Cartesian forms of
// space/irbem/ext_ostapenko.hpp's om97_basis — the printed cylindrical rows with cos phi = x/rho,
// sin phi = y/rho substituted and rho cancelled, so nothing is undefined on the dipole axis and
// there is no sqrt, no division and no branch anywhere. Row 17 carries the Schmidt-normalized
// factor sqrt(15)/8 the paper's text prescribes (the header's brief has the account).
//
// This mirrors om97_components<float> expression for expression, so a disagreement between the
// two lanes is attributable to the DEVICE and not to the arithmetic having been written
// differently on the two sides.
//
// bindings: 0 pos (3N floats, GSM, Re) | 1 out (3N floats, nT, GSM)
//           2 params (19 floats: sin psi, cos psi, A1..A17)
//           3 dims {N}
// =============================================================================================

/// The 17 harmonic amplitudes' weighted sum at a normalized SM point: the paper's Table 1 in
/// Cartesian form, weighted by par[2..18]. Shared so a total-field tracer can call it.
float3 om97_eval_sm(StructuredBuffer<float> par, float sps, float3 p) {
    const float x = p.x, y = p.y, z = p.z;
    const float x2 = x * x, y2 = y * y, z2 = z * z;
    const float rho2 = x2 + y2;
    const float rho4 = rho2 * rho2;
    const float z3 = z2 * z, z4 = z2 * z2;
    const float xy = x * y, xz = x * z, yz = y * z;
    // The Schmidt-normalized leading factors: sqrt(3), sqrt(10), sqrt(3/2), sqrt(15)/8,
    // sqrt(15)/2, sqrt(15).
    const float s3 = 1.7320508f, s10 = 3.1622777f, s32 = 1.2247449f;
    const float s15_8 = 0.48412292f, s15_2 = 1.9364917f, s15 = 3.8729833f;

    float bx = 0.0f, by = 0.0f, bz = 0.0f;
    // 1: uniform B_z
    bz += par[2];
    // 2: grad(r^3 P_3)
    bx += par[3] * (-3.0f * xz);  by += par[3] * (-3.0f * yz);  bz += par[3] * (3.0f * (z2 - 0.5f * rho2));
    // 3: grad(r^5 P_5)
    const float s3v = 2.0f * z2 - 1.5f * rho2;
    bx += par[4] * (-5.0f * xz * s3v);  by += par[4] * (-5.0f * yz * s3v);
    bz += par[4] * (5.0f * (z4 - 3.0f * rho2 * z2 + 0.375f * rho4));
    // 4, 5: ring-current rows
    bz += par[5] * rho2;
    bz += par[6] * rho4;
    // 6: (-2 rho z^3, 0, z^4)
    bx += par[7] * (-2.0f * x * z3);  by += par[7] * (-2.0f * y * z3);  bz += par[7] * z4;
    // 7: sqrt(3) (z, 0, x)
    bx += par[8] * (s3 * z);  bz += par[8] * (s3 * x);
    // 8: grad(sqrt(10) rho z (z^2 - 3 rho^2/4) cos phi)
    bx += par[9] * (s10 * z * (z2 - 0.25f * (9.0f * x2 + 3.0f * y2)));
    by += par[9] * (-1.5f * s10 * xy * z);
    bz += par[9] * (3.0f * s10 * x * (z2 - 0.25f * rho2));
    // 9: B_x = z
    bx += par[10] * z;
    // 10: B_z = rho^2 x
    bz += par[11] * (rho2 * x);
    // 11: B_x = z^3
    bx += par[12] * z3;
    // 12: (-y^2 z, x y z, -x z^2 / 2)
    bx += par[13] * (-y2 * z);  by += par[13] * (xy * z);  bz += par[13] * (-0.5f * x * z2);
    // 13: sin psi grad(r^2 P_2)
    bx += par[14] * (-x * sps);  by += par[14] * (-y * sps);  bz += par[14] * (2.0f * z * sps);
    // 14: sin psi grad(r^4 P_4)
    const float s14 = 0.5f * rho2 - 2.0f * z2;
    bx += par[15] * (3.0f * x * s14 * sps);  by += par[15] * (3.0f * y * s14 * sps);
    bz += par[15] * (2.0f * z * (2.0f * z2 - 3.0f * rho2) * sps);
    // 15: sin psi uniform B_x
    bx += par[16] * sps;
    // 16: sin psi grad(sqrt(3/2) rho (2 z^2 - rho^2/2) cos phi)
    bx += par[17] * (s32 * (2.0f * z2 - 0.5f * (3.0f * x2 + y2)) * sps);
    by += par[17] * (-s32 * xy * sps);
    bz += par[17] * (4.0f * s32 * xz * sps);
    // 17: sin psi grad((sqrt(15)/8) rho (8 z^4 - 12 rho^2 z^2 + rho^4) cos phi)
    bx += par[18] * (s15_8 * (8.0f * z4 - 12.0f * rho2 * z2 + rho4 - 24.0f * x2 * z2 + 4.0f * rho2 * x2) * sps);
    by += par[18] * (s15_2 * xy * (rho2 - 6.0f * z2) * sps);
    bz += par[18] * (s15 * xz * (4.0f * z2 - 3.0f * rho2) * sps);
    return float3(bx, by, bz);
}

/// The OM97 external field at a GSM point, in nT: GSM -> SM about y by the tilt, the paper's
/// 10 Re normalization, the harmonics, and the rotation back.
float3 om97_eval(StructuredBuffer<float> par, float3 p_gsm) {
    const float sps = par[0];
    const float cps = par[1];
    const float k = 0.1f;
    const float xs = (p_gsm.x * cps - p_gsm.z * sps) * k;
    const float ys = p_gsm.y * k;
    const float zs = (p_gsm.x * sps + p_gsm.z * cps) * k;
    const float3 b = om97_eval_sm(par, sps, float3(xs, ys, zs));
    return float3(b.x * cps + b.z * sps, b.y, -b.x * sps + b.z * cps);
}

#ifdef KERNEL_OSTAPENKO

[[vk::binding(0, 0)]] StructuredBuffer<float>   om97_pos;
[[vk::binding(1, 0)]] RWStructuredBuffer<float> om97_out;
[[vk::binding(2, 0)]] StructuredBuffer<float>   om97_params;
[[vk::binding(3, 0)]] StructuredBuffer<uint>    om97_dims;

[shader("compute")]
[numthreads(256, 1, 1)]
void irbem_om97_f32(uint3 tid : SV_DispatchThreadID) {
    const uint n = om97_dims[0];
    const uint i = tid.x;
    if (i >= n) return;   // the grid is rounded up to a multiple of 256; the tail does nothing
    const float3 b = om97_eval(om97_params,
                               float3(om97_pos[3 * i + 0], om97_pos[3 * i + 1], om97_pos[3 * i + 2]));
    om97_out[3 * i + 0] = b.x;
    om97_out[3 * i + 1] = b.y;
    om97_out[3 * i + 2] = b.z;
}

#endif  // KERNEL_OSTAPENKO
