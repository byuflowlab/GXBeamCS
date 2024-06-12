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


# -----------------------------------------------------------