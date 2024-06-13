module GXBeamCS

import ForwardDiff
using SparseArrays
using StaticArrays
using FLOWMath
using LinearAlgebra

# from common
export Material, MaterialPlane, Layer
export compliance_matrix, mass_matrix, strains_and_stresses, tsai_wu
export plotmesh, plotsoln

# from fem
export Node, MeshElement
export initialize_cache

# from afmesh
export afmesh

# from clt
export BeamSection

# shared definitions and functions
include("common.jl")

# section properties using finite element based approach
include("fem.jl")

# airfoil meshing
include("afmesh.jl")

# section properties using classical laminate theory
include("clt.jl")

end