# Timing comparison: rebuilding the load -> strain operator on every call
# (`strains_and_stresses`) versus building it once (`compliance_and_operators`) and
# evaluating a whole load history as a single matrix multiply (`strains_stresses_from_B`).
#
# Run with:  julia --project=test examples/operator_benchmark.jl

using GXBeamCS, LinearAlgebra, Random, Printf

# ---------------------------------------------------------------------------
# cross section: multi-layer closed composite pipe
# ---------------------------------------------------------------------------
function pipe_clt()
    mat1 = MaterialPlane(20.59e6, 1.42e6, 0.87e6, 0.42, 1.0)
    t = [0.1, 0.1]
    laminate1 = Layer.(Ref(mat1), t, [0.0, 90.0]*pi/180)
    laminate2 = Layer.(Ref(mat1), t, [-45.0, 45.0]*pi/180)

    s = 2.0
    r = 0.4
    b1 = BeamSection(laminate1, [-s/2, s/2], [-r, -r])
    theta = range(-pi/2, pi/2, length=20)
    b2 = BeamSection(laminate2, r*cos.(theta) .+ s/2, r*sin.(theta))
    b3 = BeamSection(laminate1, [s/2, -s/2], [r, r])
    theta = range(pi/2, 3*pi/2, length=20)
    b4 = BeamSection(laminate2, r*cos.(theta) .- s/2, r*sin.(theta))

    return CLT([b1, b2, b3, b4], true)
end

"minimum elapsed time over `n` runs, after one warm-up run"
function best(f, n)
    f()
    t = Inf
    for _ = 1:n
        t = min(t, @elapsed f())
    end
    return t
end

clt = pipe_clt()
compliance_matrix(clt)                       # `strains_and_stresses` needs the cache populated

nloc = num_strain_locs(clt)
faces = collect(1:2:nloc)                    # first face of every ply: the fatigue usage
nload = 1801                                 # one fatigue station's time history

Random.seed!(1)
Fhist = randn(nload, 3)
Mhist = randn(nload, 3)

@printf("cross section: %d strain locations, %d sections\n", nloc, length(clt.sections))
@printf("load history:  %d steps\n\n", nload)

# ---------------------------------------------------------------------------
# 1. one-time setup cost
#
# The operator build replaces `compliance_matrix`, so the fair comparison is against
# `compliance_matrix` alone, not against zero.
# ---------------------------------------------------------------------------
t_K = best(() -> compliance_matrix(pipe_clt()), 20)
t_KB_all = best(() -> compliance_and_operators(pipe_clt(), 1:nloc), 20)
t_KB_sub = best(() -> compliance_and_operators(pipe_clt(), faces; outputs=(:strain_p,)), 20)

println("SETUP (once per cross section)")
@printf("  compliance_matrix                                  %8.3f ms\n", 1e3*t_K)
@printf("  compliance_and_operators, all outputs, all locs     %8.3f ms   %.2fx\n", 1e3*t_KB_all, t_KB_all/t_K)
@printf("  compliance_and_operators, strain_p at %3d locs      %8.3f ms   %.2fx\n",
        length(faces), 1e3*t_KB_sub, t_KB_sub/t_K)

# ---------------------------------------------------------------------------
# 2. evaluating the load history
#
# NOTE ON A FAIR COMPARISON: `strains_and_stresses` always computes all four output
# fields at all `nloc` locations -- it cannot be asked for a subset. So the only
# like-for-like comparison is against a `compliance_and_operators` build that also
# covers all four fields and all locations. The narrower build below is faster partly
# because it reuses the operator and partly because it simply computes 8x less; the two
# effects are reported separately so they are not conflated.
# ---------------------------------------------------------------------------
_, _, _, Ball = compliance_and_operators(pipe_clt(), 1:nloc)
_, _, _, Bsub = compliance_and_operators(pipe_clt(), faces; outputs=(:strain_p,))

function percall()
    for j = 1:nload
        strains_and_stresses(view(Fhist, j, :), view(Mhist, j, :), clt)
    end
end

t_loop = best(percall, 3)
t_all = best(() -> strains_stresses_from_B(Fhist, Mhist, clt, Ball), 10)
t_sub = best(() -> strains_stresses_from_B(Fhist, Mhist, clt, Bsub), 20)

println("\nEVALUATION (per load history)")
@printf("  per-call loop, all outputs, all %3d locs            %8.2f ms\n", nloc, 1e3*t_loop)
@printf("  batched,       all outputs, all %3d locs            %8.3f ms\n", nloc, 1e3*t_all)
@printf("  batched,       strain_p only at %3d locs            %8.3f ms\n", length(faces), 1e3*t_sub)

println("\nSPEEDUP, decomposed")
@printf("  same work (all outputs, all locs):                 %6.0fx   <- operator reuse alone\n", t_loop/t_all)
@printf("  plus computing only strain_p at %3d locs:           %6.1fx   <- narrower scope\n",
        length(faces), t_all/t_sub)
@printf("  combined, vs what the fatigue loop does today:      %6.0fx\n", t_loop/t_sub)

# ---------------------------------------------------------------------------
# 3. total cost including setup, which is what a constraint evaluation actually pays
# ---------------------------------------------------------------------------
println("\nTOTAL per cross section (setup + one history)")
@printf("  today:      %8.2f ms\n", 1e3*(t_K + t_loop))
@printf("  operators:  %8.2f ms   %.0fx\n", 1e3*(t_KB_sub + t_sub), (t_K + t_loop)/(t_KB_sub + t_sub))

# ---------------------------------------------------------------------------
# 4. agreement check, so a timing win is never mistaken for a correct one
# ---------------------------------------------------------------------------
function worst_disagreement(clt, Bs, Fhist, Mhist)
    out = strains_stresses_from_B(Fhist, Mhist, clt, Bs)
    worst = 0.0
    for j in (1, div(size(Fhist, 1), 2), size(Fhist, 1))
        strain_b, stress_b, strain_p, stress_p =
            strains_and_stresses(view(Fhist, j, :), view(Mhist, j, :), clt)
        for (got, want) in [(out.strain_b[j, :, :], strain_b), (out.stress_b[j, :, :], stress_b),
                            (out.strain_p[j, :, :], strain_p), (out.stress_p[j, :, :], stress_p)]
            worst = max(worst, maximum(abs, got .- permutedims(want)) / maximum(abs, want))
        end
    end
    return worst
end

@printf("\nmax relative disagreement with strains_and_stresses: %.2e\n",
        worst_disagreement(clt, Ball, Fhist, Mhist))
