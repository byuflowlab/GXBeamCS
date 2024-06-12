module GXBeamCS

import ForwardDiff
using UnPack
using SparseArrays
using StaticArrays
using FLOWMath
using LinearAlgebra

# from common
export Material, MaterialPlane, Layer
export compliance_matrix, mass_matrix, strains_and_stresses, tsai_wu

# from fem
export Node, MeshElement
export initialize_cache, plotmesh
export plotsoln

# from clt
export BeamSection

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