# Unit tests for the fused compliance + load->strain/stress operators in clt.jl
using GXBeamCS, LinearAlgebra, Random, Test
using ForwardDiff, FiniteDiff
import ReverseDiff

# ---------------- fixtures ----------------
# Geometry and material properties are taken from the reference cases in test/clt.jl so the
# cross sections are known-good.

"multi-layer closed composite pipe: 4 sections, two different layups, curved segments"
function operator_pipe()
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

"single open laminated plate: exercises the closed_section=false operator branch"
function operator_plate()
    mat1 = MaterialPlane(41.5e9, 7.83e9, 3.15e9, 0.3, 1.0)
    t = fill(0.25, 8)
    theta = [25.0, 25.0, 50.0, 0.0, 50.0, 0.0, 25.0, 25.0]*pi/180
    laminate1 = Layer.(Ref(mat1), t, theta)
    return CLT([BeamSection(laminate1, [-2.0, 2.0], [0.0, 0.0])], false)
end

"closed single-laminate ring: the analogue of a turbine cylinder root, where det(cc.F) is tiny"
function operator_cylinder()
    mat1 = MaterialPlane(30.0e9, 30.0e9, 11.5e9, 0.3, 1.0)
    laminate = Layer.(Ref(mat1), fill(0.01, 4), [0.0, 45.0, -45.0, 90.0]*pi/180)
    theta = range(0, 2*pi, length=41)
    return CLT([BeamSection(laminate, cos.(theta), sin.(theta))], true)
end

const OPERATOR_FIXTURES = [("pipe", operator_pipe),
                           ("plate", operator_plate),
                           ("cylinder", operator_cylinder)]

"maximum error of `got` relative to the largest magnitude in `want`"
function relerr(got, want)
    scale = maximum(abs, want)
    scale = scale == 0 ? one(scale) : scale
    return maximum(abs, got .- want) / scale
end


@testset "clt operators" begin

# ---------------------------------------------------------------------------
@testset "K unchanged" begin
    # guards the _compliance_core! refactor: compliance_and_operators must reproduce
    # compliance_matrix exactly, not just approximately
    for (name, fixture) in OPERATOR_FIXTURES
        for shear_center in (true, false)
            Sref, scref, tcref = compliance_matrix(fixture(), shear_center)

            clt = fixture()
            S, sc, tc, _ = compliance_and_operators(clt, 1:num_strain_locs(clt);
                                                    shear_center=shear_center)

            @test S == Sref
            @test sc == scref
            @test tc == tcref

            # the cache must be populated exactly as compliance_matrix leaves it
            cltref = fixture()
            compliance_matrix(cltref, shear_center)
            @test clt.cache.F == cltref.cache.F
            @test clt.cache.L == cltref.cache.L
            @test clt.cache.Wbar == cltref.cache.Wbar
            @test clt.cache.S == cltref.cache.S
        end
    end
end

# ---------------------------------------------------------------------------
@testset "operator correctness" begin
    # the original strains_and_stresses is the oracle
    Random.seed!(1234)

    for (name, fixture) in OPERATOR_FIXTURES
        clt = fixture()
        nloc = num_strain_locs(clt)
        locs = collect(1:nloc)

        compliance_matrix(clt)
        F = randn(3)
        M = randn(3)
        strain_b, stress_b, strain_p, stress_p = strains_and_stresses(F, M, clt)

        _, _, _, Bs = compliance_and_operators(fixture(), locs)
        out = strains_stresses_from_B(F, M, clt, Bs)

        # single load case comes back as (nloc, 6), the transpose of the original layout
        @test size(out.strain_p) == (nloc, 6)
        @test relerr(out.strain_b, permutedims(strain_b)) <= 1e-12
        @test relerr(out.stress_b, permutedims(stress_b)) <= 1e-12
        @test relerr(out.strain_p, permutedims(strain_p)) <= 1e-12
        @test relerr(out.stress_p, permutedims(stress_p)) <= 1e-12

        # in ply coordinates this method only populates 11, 22 and 12; the 33, 13 and 23
        # components stay zero. (In beam coordinates they do not: the clt_Talpha rotation
        # spreads the in-plane components across all six rows for an angled element.)
        for field in (:strain_p, :stress_p)
            @test all(getfield(out, field)[:, [3, 5, 6]] .== 0)
        end

        # per-location operators: drive the four FMvec basis vectors through the *original*
        # function and check each column of B. FMvec = [F[1], M[2], -M[3], M[1]], so the
        # basis vectors in natural (F, M) coordinates are:
        basis = [([1.0, 0, 0], [0.0, 0, 0]),     # Nxbar
                 ([0.0, 0, 0], [0.0, 1, 0]),     # Mybar
                 ([0.0, 0, 0], [0.0, 0, -1]),    # Mzbar (M[3] enters negated)
                 ([0.0, 0, 0], [1.0, 0, 0])]     # Txbar

        columns = [strains_and_stresses(Fb, Mb, clt) for (Fb, Mb) in basis]
        for (fieldidx, field) in enumerate((:strain_b, :stress_b, :strain_p, :stress_p))
            # assemble what B should be, column by column, from the oracle
            for i in (1, 2, div(nloc, 2), nloc)
                want = hcat((columns[c][fieldidx][:, i] for c = 1:4)...)   # 6x4
                @test relerr(node_operator(Bs, i, field), want) <= 1e-12
            end
        end
    end
end

# ---------------------------------------------------------------------------
@testset "linearity" begin
    Random.seed!(5678)

    for (name, fixture) in OPERATOR_FIXTURES
        clt = fixture()
        _, _, _, Bs = compliance_and_operators(clt, 1:num_strain_locs(clt))

        # zero load gives exactly zero
        zero_out = strains_stresses_from_B(zeros(3), zeros(3), clt, Bs)
        for field in GXBeamCS.OPERATOR_FIELDS
            @test all(getfield(zero_out, field) .== 0)
        end

        F1 = randn(3); M1 = randn(3)
        F2 = randn(3); M2 = randn(3)

        base = strains_stresses_from_B(F1, M1, clt, Bs)
        doubled = strains_stresses_from_B(2*F1, 2*M1, clt, Bs)
        other = strains_stresses_from_B(F2, M2, clt, Bs)
        summed = strains_stresses_from_B(F1 + F2, M1 + M2, clt, Bs)

        for field in GXBeamCS.OPERATOR_FIELDS
            b = getfield(base, field)
            @test relerr(getfield(doubled, field), 2*b) <= 1e-14
            @test relerr(getfield(summed, field), b + getfield(other, field)) <= 1e-14
        end
    end
end

# ---------------------------------------------------------------------------
@testset "batched equivalence" begin
    # a load history evaluated in one shot must match the per-call loop
    Random.seed!(9012)
    nload = 200

    for (name, fixture) in OPERATOR_FIXTURES
        clt = fixture()
        nloc = num_strain_locs(clt)
        compliance_matrix(clt)

        Fhist = randn(nload, 3)
        Mhist = randn(nload, 3)

        _, _, _, Bs = compliance_and_operators(fixture(), 1:nloc)
        out = strains_stresses_from_B(Fhist, Mhist, clt, Bs)

        @test size(out.strain_p) == (nload, nloc, 6)

        for j in (1, 2, div(nload, 3), nload)
            strain_b, stress_b, strain_p, stress_p = strains_and_stresses(Fhist[j, :], Mhist[j, :], clt)
            @test relerr(out.strain_b[j, :, :], permutedims(strain_b)) <= 1e-14
            @test relerr(out.stress_b[j, :, :], permutedims(stress_b)) <= 1e-14
            @test relerr(out.strain_p[j, :, :], permutedims(strain_p)) <= 1e-14
            @test relerr(out.stress_p[j, :, :], permutedims(stress_p)) <= 1e-14
        end

        # a single case must agree with the corresponding row of the batch
        single = strains_stresses_from_B(Fhist[1, :], Mhist[1, :], clt, Bs)
        @test relerr(single.strain_p, out.strain_p[1, :, :]) <= 1e-14
    end
end

# ---------------------------------------------------------------------------
@testset "location subsets" begin
    # requesting a subset must give the same operators as the corresponding slice of the
    # full set, in the requested order
    clt = operator_pipe()
    nloc = num_strain_locs(clt)

    _, _, _, Ball = compliance_and_operators(operator_pipe(), 1:nloc)

    subset = [nloc, 3, 17, 1]
    _, _, _, Bsub = compliance_and_operators(operator_pipe(), subset)

    @test Bsub.nodes == subset
    @test length(Bsub) == 4
    for (i, idx) in enumerate(subset)
        for field in GXBeamCS.OPERATOR_FIELDS
            @test node_operator(Bsub, i, field) == node_operator(Ball, idx, field)
        end
    end

    # tuple addressing must resolve to the same operators as the equivalent indices
    tuples = [strain_loc_tuple(clt, idx) for idx in subset]
    _, _, _, Btup = compliance_and_operators(operator_pipe(), tuples)
    @test Btup.nodes == subset
    @test node_operator(Btup, 1, :strain_p) == node_operator(Bsub, 1, :strain_p)
end

# ---------------------------------------------------------------------------
@testset "outputs keyword" begin
    clt = operator_pipe()
    locs = 1:2:num_strain_locs(clt)

    _, _, _, Bfull = compliance_and_operators(operator_pipe(), locs)
    _, _, _, Bpart = compliance_and_operators(operator_pipe(), locs; outputs=(:strain_p,))

    @test operator_fields(Bfull) == GXBeamCS.OPERATOR_FIELDS
    @test operator_fields(Bpart) == (:strain_p,)
    @test Bpart.strain_b === nothing
    @test Bpart.stress_b === nothing
    @test Bpart.stress_p === nothing
    @test Bpart.strain_p == Bfull.strain_p

    out = strains_stresses_from_B(randn(3), randn(3), clt, Bpart)
    @test out.strain_p !== nothing
    @test out.strain_b === nothing
    @test out.stress_b === nothing
    @test out.stress_p === nothing
end

# ---------------------------------------------------------------------------
@testset "location addressing" begin
    for (name, fixture) in OPERATOR_FIXTURES
        clt = fixture()
        nloc = num_strain_locs(clt)

        # index <-> tuple round trip over every location
        for idx = 1:nloc
            s, e, p, f = strain_loc_tuple(clt, idx)
            @test strain_loc_index(clt, s, e, p, f) == idx
        end

        # per-section ranges tile 1:nloc, and each section's face==1 locations are its 1:2:N
        covered = Int[]
        for i in eachindex(clt.sections)
            rng = strain_loc_indices(clt, i)
            append!(covered, rng)

            faces = [strain_loc_tuple(clt, idx)[4] for idx in rng]
            @test faces[1:2:end] == fill(1, length(1:2:length(rng)))
            @test faces[2:2:end] == fill(2, length(2:2:length(rng)))
        end
        @test covered == collect(1:nloc)
    end

    clt = operator_pipe()
    nloc = num_strain_locs(clt)
    @test_throws ArgumentError strain_loc_tuple(clt, 0)
    @test_throws ArgumentError strain_loc_tuple(clt, nloc + 1)
    @test_throws ArgumentError strain_loc_index(clt, length(clt.sections) + 1, 1, 1, 1)
    @test_throws ArgumentError strain_loc_index(clt, 1, 1, 1, 3)
    @test_throws ArgumentError strain_loc_index(clt, 1, 1, 99, 1)
end

# ---------------------------------------------------------------------------
@testset "errors" begin
    clt = operator_pipe()
    nloc = num_strain_locs(clt)

    @test_throws ArgumentError compliance_and_operators(operator_pipe(), [nloc + 1])
    @test_throws ArgumentError compliance_and_operators(operator_pipe(), [1]; outputs=(:bogus,))
    @test_throws ArgumentError compliance_and_operators(operator_pipe(), [1]; outputs=())

    _, _, _, Bs = compliance_and_operators(operator_pipe(), 1:8; outputs=(:strain_p,))

    # asking for a field that was not built
    @test_throws ArgumentError node_operator(Bs, 1, :stress_b)
    @test_throws ArgumentError node_operator(Bs, 1, :bogus)
    @test_throws ArgumentError node_operator(Bs, 0, :strain_p)
    @test_throws ArgumentError node_operator(Bs, 9, :strain_p)

    # mismatched or wrongly shaped load batches
    @test_throws DimensionMismatch strains_stresses_from_B(randn(4, 3), randn(5, 3), clt, Bs)
    @test_throws DimensionMismatch strains_stresses_from_B(randn(4, 2), randn(4, 2), clt, Bs)
    @test_throws DimensionMismatch strains_stresses_from_B(randn(3), randn(4, 3), clt, Bs)
    @test_throws DimensionMismatch strains_stresses_from_B(randn(2), randn(2), clt, Bs)

    # operators built for a bigger cross section than the one being evaluated
    _, _, _, Bbig = compliance_and_operators(operator_pipe(), [nloc])
    @test_throws ArgumentError strains_stresses_from_B(randn(3), randn(3), operator_plate(), Bbig)
end

# ---------------------------------------------------------------------------
@testset "AD compatibility" begin
    # design variables: ply thickness, ply angle, E1, and a section coordinate
    xvec = [0.1, -45.0*pi/180, 20.59e6, 0.4]

    function operatorwrapper(x)
        TF = eltype(x)
        mat1 = MaterialPlane(x[3], TF(1.42e6), TF(0.87e6), TF(0.42), one(TF))

        t = TF[x[1], 0.1]
        laminate1 = Layer.(Ref(mat1), t, TF[0.0, pi/2])
        laminate2 = Layer.(Ref(mat1), t, TF[x[2], -x[2]])

        s = TF(2.0)
        r = x[4]
        b1 = BeamSection(laminate1, TF[-s/2, s/2], TF[-r, -r])
        theta = range(-pi/2, pi/2, length=8)
        b2 = BeamSection(laminate2, r*cos.(theta) .+ s/2, r*sin.(theta))
        b3 = BeamSection(laminate1, TF[s/2, -s/2], TF[r, r])
        theta = range(pi/2, 3*pi/2, length=8)
        b4 = BeamSection(laminate2, r*cos.(theta) .- s/2, r*sin.(theta))

        clt = CLT([b1, b2, b3, b4], true)
        locs = 1:2:num_strain_locs(clt)
        S, _, _, Bs = compliance_and_operators(clt, locs; outputs=(:strain_p,))

        F = TF[1.0e3, 0, 0]
        M = TF[0.0, -1.0e3, 5.0e2]
        out = strains_stresses_from_B(F, M, clt, Bs)

        return vcat(vec(Matrix(S)), vec(out.strain_p))
    end

    # the wrapper must be usable, and its plain-Float64 result must agree with the
    # per-call strains_and_stresses path
    y = operatorwrapper(xvec)
    @test all(isfinite, y)

    Jfd = FiniteDiff.finite_difference_jacobian(operatorwrapper, xvec)
    Jfor = ForwardDiff.jacobian(operatorwrapper, xvec)
    @test size(Jfor) == size(Jfd)
    @test relerr(Jfor, Jfd) <= 1e-6

    Jrev = ReverseDiff.jacobian(operatorwrapper, xvec)
    @test relerr(Jrev, Jfor) <= 1e-10
end

end # clt operators
