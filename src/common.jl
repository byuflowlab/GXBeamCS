# ------ material properties ------

"""
    Material(E1, E2, E3, G12, G13, G23, nu12, nu13, nu23, rho,
        S1t=1.0, S1c=1.0, S2t=1.0, S2c=1.0, S3t=1.0, S3c=1.0, S12=1.0, S13=1.0, S23=1.0)

Orthotropic material properties. Axis 1 is the main ply axis, axis 2 is the transverse
ply axis, and axis 3 is normal to the ply.  For a fiber orientation of zero, axis 1 is
along the beam axis.

# Arguments
- `Ei::float`: Young's modulus along ith axis.
- `Gij::float`: Shear moduli
- `nuij::float`: Poisson's ratio.  ``nu_ij E_j = nu_ji E_i``
- `rho::float`: Density
- `Sit::float`: Tensile strength in ith direction
- `Sic::float`: Compressive strength in ith direction
- `Sij::float`: Strength in ij direction
"""
struct Material{TF}
    E1::TF
    E2::TF
    E3::TF
    G12::TF
    G13::TF
    G23::TF
    nu12::TF
    nu13::TF
    nu23::TF
    rho::TF
    S1t::TF
    S1c::TF
    S2t::TF
    S2c::TF
    S3t::TF
    S3c::TF
    S12::TF
    S13::TF
    S23::TF
end

# floating point type
Base.eltype(::Material{TF}) where TF = TF
Base.eltype(::Type{Material{TF}}) where TF = TF

# floating point type conversions
Material{TF}(m::Material) where {TF} = Material{TF}(
    m.E1, m.E2, m.E3, m.G12, m.G13, m.G23, m.nu12, m.nu13, m.nu23, m.rho,
    m.S1t, m.S1c, m.S2t, m.S2c, m.S3t, m.S3c, m.S12, m.S13, m.S23)
Base.convert(::Type{Material{TF}}, m::Material) where {TF} = Material{TF}(m)

# promote all inputs to the same type
function Material(E1, E2, E3, G12, G13, G23, nu12, nu13, nu23, rho,
    S1t, S1c, S2t, S2c, S3t, S3c, S12, S13, S23)

    if S1c < 0 || S2c < 0 || S3c < 0
        @warn "the compressive strengths should be input as positive numbers"
    end

    return Material(promote(E1, E2, E3, G12, G13, G23, nu12, nu13, nu23, rho,
    S1t, S1c, S2t, S2c, S3t, S3c, S12, S13, S23)...)
end

# make strength properties optional
Material(E1, E2, E3, G12, G13, G23, nu12, nu13, nu23, rho) = Material(promote(E1, E2, E3, G12, G13, G23, nu12, nu13, nu23, rho, ones(9)...)...)
Material{TF}(E1, E2, E3, G12, G13, G23, nu12, nu13, nu23, rho) where TF = Material{TF}(E1, E2, E3, G12, G13, G23, nu12, nu13, nu23, rho, ones(TF, 9)...)

# for a plane stress case we have fewer constants

MaterialPlane(E1, E2, G12, nu12, rho) = MaterialPlane(E1, E2, G12, nu12, rho, 1.0, 1.0, 1.0, 1.0, 1.0)

"""
    MaterialPlane(E1, E2, G12, nu12, rho,
        S1t=1.0, S1c=1.0, S2t=1.0, S2c=1.0, S12=1.0)

Orthotropic material properties for a plane stress case (which assumes no stress along axis 3).
Axis 1 is the main ply axis, axis 2 is the transverse
ply axis, and axis 3 is normal to the ply.  For a fiber orientation of zero, axis 1 is
along the beam axis.

# Arguments
- `Ei::float`: Young's modulus along ith axis.
- `Gij::float`: Shear moduli
- `nuij::float`: Poisson's ratio.  ``nu_ij E_j = nu_ji E_i``
- `rho::float`: Density
- `Sit::float`: Tensile strength in ith direction
- `Sic::float`: Compressive strength in ith direction
- `Sij::float`: Strength in ij direction
"""
MaterialPlane(E1, E2, G12, nu12, rho, S1t, S1c, S2t, S2c, S12) = Material(E1, E2, 0.0, G12, 0.0, 0.0, nu12, 0.0, 0.0, rho,
                                                                            S1t, S1c, S2t, S2c, 0.0, 0.0, S12, 0.0, 0.0)

# --------------------------

# --------- composite layer -------------


"""
    Layer(material, t, theta)

A layer (could be one ply or many plys of same material).
A layup is a vector of layers.

**Arguments**
- `material::Material`: material of layer
- `t::float`: thickness of layer
- `theta::float`: fiber orientation (rad)
"""
struct Layer{TF}
    material::Material{TF}
    t::TF
    theta::TF
end
Base.eltype(::Layer{TF}) where {TF} = TF
Base.eltype(::Type{Layer{TF}}) where {TF} = TF

Layer{TF}(l::Layer) where {TF} = Layer{TF}(l.material, l.t, l.theta)
Base.convert(::Type{Layer{TF}}, l::Layer) where {TF} = Layer{TF}(l)

function thickness(layup::Vector{<:Layer})
    return sum(l.t for l in layup)
end

function get_laminate_type(laminate::Vector{<:Layer})
    return typeof(laminate[1].t)
end

# -----------------------------------------------------------


abstract type CompositeSectionAnalysis end

"""
    compliance_matrix(::CompositeSectionAnalysis, shear_center=true)

Compute compliance matrix for the section.

# Arguments
- `shear_center::Bool`: Indicates whether the compliance matrix should be provided about the
    shear center

# Returns
- `S::Matrix`: compliance matrix in the order expected by GXBeam.
- `sc::Vector{float}`: x, y location of shear center (location where a transverse/shear
    force will not produce any torsion, i.e., beam will not twist)
- `tc::Vector{float}`: x, y location of tension center, aka elastic center, aka centroid
    (location where an axial force will not produce any bending, i.e., beam will remain
    straight)
"""
function compliance_matrix(::CompositeSectionAnalysis, shear_center=true)
    S = zeros(6, 6)
    sc = [0.0, 0.0]
    tc = [0.0, 0.0]
    return S, sc, tc
end

"""
    mass_matrix(::CompositeSectionAnalysis)

Compute mass matrix for the section

# Returns
- `M::Matrix`: mass matrix in the order expected by GXBeam.
- `mc::Vector{float}`: x, y location of mass center
"""
function mass_matrix(::CompositeSectionAnalysis)
    M = zeros(6, 6)
    mc = [0.0, 0.0]
    return M, mc
end

"""
    strains_and_stresses(F, M, ::CompositeSectionAnalysis)

Compute stresses and strains across the cross section (locations depend on the method).
Uses GXBeam local beam axis (e.g. axial stresses would correspond to the `xx` direction, or first index).

# Arguments
- `F::Vector(3)`: force at this cross section in x, y, z directions.  Using GXBeam local beam axis.
- `M::Vector(3)`: moment at this cross section in x, y, z directions. Using GXBeam local beam axis.

# Returns
- `strain_b::Vector(6, nloc)`: strains in beam coordinate system for each element. order: xx, yy, zz, xy, xz, yz
    the meaning of the locations varies with the different methods (e.g., elements vs edges of plys)
- `stress_b::Vector(6, nloc)`: stresses in beam coordinate system for each element. order: xx, yy, zz, xy, xz, yz
- `strain_p::Vector(6, nloc)`: strains in ply coordinate system for each element. order: 11, 22, 33, 12, 13, 23
- `stress_p::Vector(6, nloc)`: stresses in ply coordinate system for each element. order: 11, 22, 33, 12, 13, 23
"""
function strains_and_stresses(F, M, ::CompositeSectionAnalysis)
    return zeros(6, 1), zeros(6, 1), zeros(6, 1), zeros(6, 1)
end


"""
    tsai_wu(stress_p, ::CompositeSectionAnalysis)

Tsai Wu failure criteria

# Arguments
- `stress_p::vector(6, nloc)`: stresses in ply coordinate system

# Returns
- `failure::vector(nloc)`: tsai-wu failure criteria at each location.  fails if >= 1
"""
function tsai_wu(stress_p, ::CompositeSectionAnalysis)
    return zeros(1)
end



"""
    plotgeometry(::CompositeSectionAnalysis, pyplot; plotnumbers=false)

plot geometry for a quick visualization.
Need to pass in a PyPlot object as PyPlot is not loaded by this package.
"""
function plotgeometry(::CompositeSectionAnalysis, pyplot; plotnumbers=false) end

"""
    plotsoln(::CompositeSectionAnalysis, soln, pyplot)

plot stress/strain on mesh
soln could be any vector that is of appropriate length for the analysis method, e.g., sigma_b[3, :]
Need to pass in a PyPlot object as PyPlot is not loaded by this package.
"""
function plotsoln(::CompositeSectionAnalysis, soln, pyplot) end