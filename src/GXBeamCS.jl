module GXBeamCS

import ForwardDiff
using UnPack
using SparseArrays
using StaticArrays
using FLOWMath
using LinearAlgebra

# from common
export Material, MaterialPlane, Layer

# from fem
export Node, MeshElement
export initialize_cache, compliance_matrix, mass_matrix, plotmesh
export strain_recovery, plotsoln, tsai_wu

# from clt
export MaterialPlane, BeamSection
export beamstiffness

# from afmesh
export afmesh

# shared definitions and functions
include("common.jl")

# section properties using finite element based approach
include("fem.jl")

# section properties using classical laminate theory
include("clt.jl")

# airfoil meshing
include("afmesh.jl")

end